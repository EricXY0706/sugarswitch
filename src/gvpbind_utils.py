"""
GVP-Bind: per-residue protein binding-site (interface) prediction.
Credit to: https://github.com/LAJ-THU/GVP-Bind
"""
from __future__ import annotations
import csv

import numpy as np
import torch

from GVP_Bind.src.gvpbind.data.dataset import build_knn_graph
from GVP_Bind.src.gvpbind.data.dockground import parse_pdb_complex
from GVP_Bind.src.gvpbind.infer.atomgraph_infer import compute_atomgraph
from GVP_Bind.src.gvpbind.infer.saprot_embed import compute_saprot_embedding
from GVP_Bind.src.gvpbind.task.multitask import MultiTask

_ONE_LETTER = "ARNDCQEGHILKMFPSTWYV"
DEFAULT_CHECKPOINT = "./GVP_Bind/checkpoints/gvpbind_saprot_atg_dcl.ckpt"

def _model_cfg(module) -> dict:
    hp = getattr(module, "hparams", {}) or {}
    cfg = hp.get("model_cfg", hp)
    return dict(cfg) if cfg is not None else {}

def gvpbind_predict(
    structure_file: str,
    chain_id: str,
    out_dir: str,
    filename: str,
):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    parsed = parse_pdb_complex(structure_file, query_chain=chain_id)

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

    module = MultiTask.load_from_checkpoint(DEFAULT_CHECKPOINT, map_location=device, strict=False)
    module.eval().to(device)
    cfg = _model_cfg(module)
    need_saprot = int(cfg.get("esm_dim", 0)) > 0
    need_atom = bool(cfg.get("atomgraph", False))

    batch = {
        "coords": coords.to(device),
        "res_types": res_types.to(device),
        "chain_ids": chain_ids.to(device),
        "is_query": is_query.to(device),
        "edge_index": edge_index.to(device),
        "seq_idx": seq_idx.to(device),
    }

    if need_saprot:
        esm = compute_saprot_embedding(
            structure_file, chain_id, expected_aa_seq=aa_seq,
            device=device)
        if esm.shape[0] != nq:
            raise SystemExit(f"SaProt L={esm.shape[0]} != query residues {nq}")
        batch["esm_features"] = esm.to(device)

    if need_atom:
        atom = compute_atomgraph(structure_file, chain_id, q_resnum)
        for k, v in atom.items():
            batch[k] = v.to(device)

    with torch.no_grad():
        out = module(batch)
    binary_logit = out["binary_logit"].cpu()
    prob = torch.sigmoid(binary_logit)

    ranks = torch.empty_like(binary_logit)
    order = binary_logit.argsort()
    n = max(1, nq - 1)
    ranks[order] = torch.arange(nq, dtype=ranks.dtype) / n
    prob = ranks

    out_path = f"{out_dir}/{filename}_gvpbind_results.csv"
    prob_by_res = {}
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Model", "Chain", "Residue Index", "Sequence", "Binding site probability"])
        for i in range(nq):
            w.writerow(["0", chain_id, int(seq_idx[i]) + 1,
                        aa_seq[i], f"{float(prob[i]):.4f}"])
            prob_by_res[int(q_resnum[i])] = float(prob[i])
    
    return prob_by_res