"""Compute the heavy-atom point cloud inline from a bare PDB/CIF.

Replicates ``scripts/compute_atomgraph.py`` (which reads a cache pickle for the
residue map) but instead takes the query residue order directly from the live
parse, so it works on novel structures with no cache. Emits the exact arrays the
AtomEncoder expects, aligned to the QUERY residues (``resid`` in ``[0, nq)``).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import torch

from ..data.dataset import build_knn_graph

# identical to scripts/compute_atomgraph.py
_BB = {"N", "CA", "C", "O"}
_ELEM = {"C": 0, "N": 1, "O": 2, "S": 3}
_ATOM_K = 16


def compute_atomgraph(
    pdb_path: str | Path,
    chain: str,
    q_resnums: np.ndarray,
    *,
    atom_k: int = _ATOM_K,
) -> dict:
    """Return atom-graph tensors for ``chain``'s heavy atoms.

    ``q_resnums`` is the array of PDB residue numbers of the QUERY residues in
    model order (``parsed['pdb_resnum'][is_query]``). Atoms are mapped back to
    residue index via this order, so ``atom_resid`` aligns 1:1 with the model's
    query residues — exactly the contract the training cache guarantees under
    ``drop_partner_chain=ON``.

    Keys: ``atom_xyz``[A,3] f32, ``atom_elem``[A] i64, ``atom_scflag``[A] f32,
    ``atom_resid``[A] i64, ``atom_edge_index``[2,E] i64.
    """
    import gemmi

    rn_to_idx = {int(rn): i for i, rn in enumerate(q_resnums)}
    nq = len(q_resnums)

    st = gemmi.read_structure(str(pdb_path))
    if len(st) == 0:
        raise RuntimeError(f"empty structure: {pdb_path}")
    ch = next((c for c in st[0] if c.name == chain), None)
    if ch is None:
        raise RuntimeError(f"chain {chain} not found in {pdb_path}")

    xyz, elem, scflag, resid = [], [], [], []
    for res in ch:
        if res.het_flag == "H":
            continue
        ri = rn_to_idx.get(res.seqid.num)
        if ri is None:
            continue
        for a in res:
            if a.element.is_hydrogen:
                continue
            xyz.append([a.pos.x, a.pos.y, a.pos.z])
            elem.append(_ELEM.get(a.element.name, 4))
            scflag.append(0 if a.name.strip() in _BB else 1)
            resid.append(ri)
    if not xyz:
        raise RuntimeError(f"no heavy atoms mapped for chain {chain} of {pdb_path}")

    axyz = torch.tensor(np.asarray(xyz, np.float32))
    resid_t = torch.tensor(np.asarray(resid, np.int64))
    if int(resid_t.max()) >= nq:
        raise RuntimeError(
            f"atom resid out of range ({int(resid_t.max())} >= nq={nq}); "
            "query residue map is inconsistent with the structure")
    return {
        "atom_xyz": axyz,
        "atom_elem": torch.tensor(np.asarray(elem, np.int64)),
        "atom_scflag": torch.tensor(np.asarray(scflag, np.float32)),
        "atom_resid": resid_t,
        "atom_edge_index": build_knn_graph(axyz, k=atom_k),
    }
