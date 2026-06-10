"""Binary-head predictor for SaProt + atom-graph checkpoints, on novel PDBs.

Same role and output as ``predict_binary`` (a per-residue ``Binding site
probability`` CSV that proseek reads), but it additionally computes the two
cached-at-training features INLINE so a SaProt/atomgraph checkpoint can run with
no precomputed cache:

  * SaProt 1280-d PLM embedding  (foldseek 3Di -> SaProt_650M_AF2)
  * heavy-atom graph             (gemmi)

Which features are computed is auto-detected from the checkpoint's
``model_cfg`` hparams (``esm_dim``, ``atomgraph``) — no flags needed; a plain
(no-PLM, no-atomgraph) checkpoint runs exactly like ``predict_binary``.

Inference regime: these checkpoints are trained with ``drop_partner_chain=ON``
(the atom-graph ``resid`` requires it), i.e. an **apo / query-chain-only**
predictor. We therefore parse the complex, then slice to the query chain (remap
it to chain 0) and build all features over the query residues only.

Usage:
    gvpbind-predict <pdb>_<chain> \
        [--checkpoint CKPT]  (default: bundled GVP+SaProt+atom-graph model) \
        --name RUN --predictions_folder OUT [--temperature T | --rank-normalize] \
        [--saprot-model-dir DIR] [--foldseek-bin PATH]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import torch

from ..data.dataset import build_knn_graph
from ..data.dockground import parse_pdb_complex
from ..infer.atomgraph_infer import compute_atomgraph
from ..infer.saprot_embed import DEFAULT_FOLDSEEK, DEFAULT_MODEL_DIR, compute_saprot_embedding
from ..task.multitask import MultiTask
from ._annotate import write_annotated_pdb, write_chimerax_script

_ONE_LETTER = "ARNDCQEGHILKMFPSTWYV"

# Bundled best model (GVP + SaProt + atom-graph, DCL contrastive). Resolved
# relative to the repo root so an editable install (`pip install -e .`) finds it
# with no flags. Override with --checkpoint.
DEFAULT_CHECKPOINT = (
    Path(__file__).resolve().parents[3] / "checkpoints" / "gvpbind_saprot_atg_dcl.ckpt"
)


def _split_target(target: str) -> tuple[Path, str]:
    if "_" not in target:
        raise SystemExit(f"target must look like 'path/to.pdb_C', got {target!r}")
    pdb_str, chain = target.rsplit("_", 1)
    return Path(pdb_str), chain


def _model_cfg(module) -> dict:
    hp = getattr(module, "hparams", {}) or {}
    cfg = hp.get("model_cfg", hp)
    return dict(cfg) if cfg is not None else {}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("target", help="<pdb_path>_<query_chain>")
    p.add_argument("--mode", default="interface", help="(ignored, compat)")
    p.add_argument("--name", default="scannetpp")
    p.add_argument("--predictions_folder", default="./predictions")
    p.add_argument("--pdb", action="store_true", help="(ignored, compat)")
    p.add_argument("--noMSA", action="store_true", help="(ignored, compat)")
    p.add_argument("--checkpoint", default=str(DEFAULT_CHECKPOINT),
                   help="model checkpoint (default: bundled gvpbind_saprot_atg_dcl.ckpt)")
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--temperature", type=float, default=1.0)
    p.add_argument("--rank-normalize", action="store_true")
    p.add_argument("--saprot-model-dir", default=DEFAULT_MODEL_DIR)
    p.add_argument("--foldseek-bin", default=DEFAULT_FOLDSEEK)
    args = p.parse_args()
    if args.temperature <= 0:
        raise SystemExit("--temperature must be > 0")
    if not Path(args.checkpoint).exists():
        raise SystemExit(
            f"checkpoint not found: {args.checkpoint}\n"
            "Pass --checkpoint, or install from the repo root so the bundled "
            "checkpoint under checkpoints/ is on disk.")

    pdb_path, query_chain = _split_target(args.target)
    out_dir = Path(args.predictions_folder)
    out_dir.mkdir(parents=True, exist_ok=True)

    parsed = parse_pdb_complex(pdb_path, query_chain=query_chain)
    if parsed is None:
        raise SystemExit(f"Could not parse {pdb_path} (or chain {query_chain} missing).")

    # --- slice to the query chain (apo / drop_partner regime) -----------------
    q = np.asarray(parsed["is_query"], dtype=bool)
    coords = torch.from_numpy(parsed["coords"][q])
    res_types = torch.from_numpy(parsed["res_types"][q].astype("int64"))
    seq_idx = torch.from_numpy(parsed["seq_idx"][q].astype("int64"))
    q_resnum = parsed["pdb_resnum"][q]
    nq = int(q.sum())
    chain_ids = torch.zeros(nq, dtype=torch.long)        # query remapped to chain 0
    is_query = torch.ones(nq, dtype=torch.float32)
    edge_index = build_knn_graph(coords[:, 1])
    aa_seq = "".join(_ONE_LETTER[i] if i < 20 else "X" for i in res_types.tolist())

    module = MultiTask.load_from_checkpoint(args.checkpoint, map_location=args.device, strict=False)
    module.eval().to(args.device)
    cfg = _model_cfg(module)
    need_saprot = int(cfg.get("esm_dim", 0)) > 0
    need_atom = bool(cfg.get("atomgraph", False))

    batch = {
        "coords": coords.to(args.device),
        "res_types": res_types.to(args.device),
        "chain_ids": chain_ids.to(args.device),
        "is_query": is_query.to(args.device),
        "edge_index": edge_index.to(args.device),
        "seq_idx": seq_idx.to(args.device),
    }

    if need_saprot:
        esm = compute_saprot_embedding(
            pdb_path, query_chain, expected_aa_seq=aa_seq,
            model_dir=args.saprot_model_dir, foldseek_bin=args.foldseek_bin,
            device=args.device)
        if esm.shape[0] != nq:
            raise SystemExit(f"SaProt L={esm.shape[0]} != query residues {nq}")
        batch["esm_features"] = esm.to(args.device)

    if need_atom:
        atom = compute_atomgraph(pdb_path, query_chain, q_resnum)
        for k, v in atom.items():
            batch[k] = v.to(args.device)

    with torch.no_grad():
        out = module(batch)
    binary_logit = out["binary_logit"].cpu()
    prob = torch.sigmoid(binary_logit / args.temperature)

    if args.rank_normalize:
        ranks = torch.empty_like(binary_logit)
        order = binary_logit.argsort()
        n = max(1, nq - 1)
        ranks[order] = torch.arange(nq, dtype=ranks.dtype) / n
        prob = ranks

    chain_letters = parsed["meta"]["chains"]
    qchar = query_chain
    out_path = out_dir / f"predictions_{args.name}.csv"
    prob_by_res = {}
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Model", "Chain", "Residue Index", "Sequence", "Binding site probability"])
        for i in range(nq):
            w.writerow(["0", qchar, int(seq_idx[i]) + 1,
                        aa_seq[i], f"{float(prob[i]):.4f}"])
            prob_by_res[(qchar, int(q_resnum[i]))] = float(prob[i])
    print(f"wrote {out_path}  (saprot={need_saprot} atomgraph={need_atom} L={nq})")

    annotated_pdb = out_dir / f"annotated_{args.name}.pdb"
    write_annotated_pdb(pdb_path, annotated_pdb, prob_by_res)
    write_chimerax_script(out_dir / f"annotated_{args.name}.cxc", annotated_pdb.name)
    print(f"wrote {annotated_pdb}")


if __name__ == "__main__":
    main()
