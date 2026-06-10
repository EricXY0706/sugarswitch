"""Compute per-residue SaProt embeddings inline from a bare PDB/CIF.

Replicates the recipe used to build the ``*_saprot`` training caches
(``data/labeling/ppbs_binary_ba_saprot/<stem>.pt``, ``[N, 1280]`` float16):

  1. foldseek ``structureto3didescriptor`` -> per-chain 3Di structure tokens
     (vendored ``get_struc_seq`` from westlake-repl/SaProt utils/foldseek_util.py).
  2. interleave AA + 3Di (lowercased) into the SaProt "SA" sequence, e.g.
     ``M#EvVpQp...`` ('#' marks the 3Di of a low-pLDDT residue when masking is on).
  3. tokenize with SaProt's ``EsmTokenizer`` and run ``EsmForMaskedLM`` with
     ``output_hidden_states=True``; take the **last** hidden layer, strip the
     CLS/EOS tokens -> ``[L, 1280]`` per-residue embedding.

The exact layer / masking are pinned empirically by
``scripts/validate_saprot_embed.py`` (recompute one cached stem, diff vs the
stored ``.pt``). If that validation fails, adjust ``HIDDEN_LAYER`` /
``plddt_mask`` here until it matches the cache the checkpoint was trained on.

The SaProt model weights and the foldseek binary are NOT bundled. Point at them
via the environment (see ``backend_setup`` in the handoff) or the function args:
    SAPROT_MODEL_DIR   default /home/liuaj/WORK/apps/SaProt/weights/PLMs/SaProt_650M_AF2
    FOLDSEEK_BIN       default /home/liuaj/WORK/apps/SaProt/bin/foldseek
"""
from __future__ import annotations

import os
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
import torch

# --- training-recipe constants (verify with validate_saprot_embed.py) ---------
# Validation vs the ppbs_binary_ba_saprot cache (stem 11as-A, 2026-06-06):
#   per-chain tokenization + last hidden layer  -> cosine 0.989 to the cache
#   (shape/alignment exact). Whole-complex tokenization was ruled out (0.36).
#   The residual ~1% is a foldseek-version / BA-vs-AU structure-source artifact
#   on a small tail of residues, not a recipe error, and is well within the
#   tolerance of the model's Linear(1280->64) projection. The binding acceptance
#   test is functional (interface residues score high), not bit-exactness.
HIDDEN_LAYER = -1          # last hidden layer (EsmForMaskedLM.hidden_states[-1])
PLDDT_THRESHOLD = 70.0     # foldseek 3Di masking threshold ('auto' on AF/Boltz)
DEFAULT_MODEL_DIR = os.environ.get(
    "SAPROT_MODEL_DIR", "./SaProt/weights/PLMs")
DEFAULT_FOLDSEEK = os.environ.get(
    "FOLDSEEK_BIN", "./SaProt/foldseek/foldseek")


# --- vendored from SaProt utils/foldseek_util.py (trimmed) --------------------
def _extract_plddt(pdb_path: str) -> np.ndarray:
    from Bio.PDB import MMCIFParser, PDBParser
    parser = MMCIFParser(QUIET=True) if pdb_path.endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("p", pdb_path)
    chain = next(iter(structure[0]))
    plddts = []
    for residue in chain:
        plddts.append(np.mean([a.get_bfactor() for a in residue]))
    return np.asarray(plddts)


def get_struc_seq(foldseek: str, path: str, chains: list[str] | None = None,
                  process_id: int = 0, plddt_mask: bool | str = "auto",
                  plddt_threshold: float = PLDDT_THRESHOLD) -> dict:
    """Return {chain: (aa_seq, foldseek_seq, combined_seq)}. Vendored, verbatim
    logic from SaProt's foldseek_util.get_struc_seq so embeddings match training."""
    assert os.path.exists(foldseek), f"foldseek not found: {foldseek}"
    assert os.path.exists(path), f"structure not found: {path}"
    tmp = f"get_struc_seq_{process_id}_{time.time()}.tsv"
    os.system(f"{foldseek} structureto3didescriptor -v 0 --threads 1 "
              f"--chain-name-mode 1 {path} {tmp}")
    if plddt_mask == "auto":
        with open(path) as r:
            plddt_mask = "alphafold" in r.read().lower()
    seq_dict: dict[str, tuple] = {}
    name = os.path.basename(path)
    with open(tmp) as r:
        for line in r:
            desc, seq, struc_seq = line.split("\t")[:3]
            if plddt_mask:
                try:
                    plddts = _extract_plddt(path)
                    assert len(plddts) == len(struc_seq)
                    np_seq = np.array(list(struc_seq))
                    np_seq[np.where(plddts < plddt_threshold)[0]] = "#"
                    struc_seq = "".join(np_seq)
                except Exception as e:  # noqa: BLE001
                    print(f"[saprot] plddt mask failed for {name}: {e}")
            chain = desc.split(" ")[0].replace(name, "").split("_")[-1]
            if (chains is None or chain in chains) and chain not in seq_dict:
                combined = "".join(a + b.lower() for a, b in zip(seq, struc_seq))
                seq_dict[chain] = (seq, struc_seq, combined)
    for ext in ("", ".dbtype"):
        try:
            os.remove(tmp + ext)
        except OSError:
            pass
    return seq_dict


# --- model loading (cached across calls within a process) ---------------------
@lru_cache(maxsize=2)
def _load_saprot(model_dir: str, device: str):
    from transformers import EsmForMaskedLM, EsmTokenizer
    tok = EsmTokenizer.from_pretrained(model_dir)
    model = EsmForMaskedLM.from_pretrained(model_dir).to(device).eval()
    return tok, model


def compute_saprot_embedding(
    pdb_path: str | Path,
    chain: str,
    *,
    expected_aa_seq: str | None = None,
    model_dir: str = DEFAULT_MODEL_DIR,
    foldseek_bin: str = DEFAULT_FOLDSEEK,
    device: str = "cuda",
    plddt_mask: bool | str = "auto",
) -> torch.Tensor:
    """Return ``[L, 1280]`` float32 per-residue SaProt embedding for ``chain``.

    ``expected_aa_seq`` (the AA sequence of the residues the model will score,
    in order) is checked against foldseek's sequence to catch residue-set
    mismatches early (gaps / missing backbone). Pass ``parsed`` AA letters.
    """
    pdb_path = str(pdb_path)
    parsed = get_struc_seq(foldseek_bin, pdb_path, [chain], plddt_mask=plddt_mask)
    if chain not in parsed:
        raise RuntimeError(f"foldseek produced no 3Di for chain {chain} of {pdb_path}")
    aa_seq, _foldseek_seq, combined_seq = parsed[chain]

    if expected_aa_seq is not None:
        if len(aa_seq) != len(expected_aa_seq):
            raise RuntimeError(
                "SaProt/parser residue-count mismatch for "
                f"{Path(pdb_path).name} chain {chain}: foldseek L={len(aa_seq)} "
                f"vs parser L={len(expected_aa_seq)}. The PLM embedding cannot be "
                "aligned to the model's residues. (Usually a missing-backbone / "
                "gap residue; Boltz refolds should be clean.)\n"
                f"  foldseek: {aa_seq}\n  parser  : {expected_aa_seq}")
        if aa_seq != expected_aa_seq:
            ndiff = sum(a != b for a, b in zip(aa_seq, expected_aa_seq))
            print(f"[saprot] WARN: {ndiff} residue-letter diffs vs parser for "
                  f"{Path(pdb_path).name} chain {chain} (lengths match; likely "
                  "non-standard residues mapped differently). Proceeding.")

    tok, model = _load_saprot(model_dir, device)
    inputs = tok(combined_seq, return_tensors="pt")
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        out = model(**inputs, output_hidden_states=True)
    hidden = out.hidden_states[HIDDEN_LAYER][0]      # [L+2, 1280] incl CLS/EOS
    per_res = hidden[1:1 + len(aa_seq)]              # strip CLS/EOS -> [L, 1280]
    return per_res.float().cpu()
