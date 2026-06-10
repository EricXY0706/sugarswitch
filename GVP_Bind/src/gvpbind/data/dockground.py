"""Dockground complex parsing.

Dockground (https://dockground.compbio.ku.edu/) provides curated
non-redundant protein-protein complexes. We use the Bound benchmark v4.

This module assumes the raw PDB files have already been downloaded into
`{data_root}/dockground/raw/` (the actual download is a one-shot manual
step — see README) and provides a parser to convert each PDB into the
per-complex .pt schema used by `dataset.ComplexDataset`.

We deliberately keep the download out of the Python path (curl/wget in a
shell script) because the Dockground server prefers session cookies and
is best driven from a one-off bash script the user runs interactively.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

# 3-letter -> 1-letter mapping; X for unknown.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA_TO_IDX = {a: i for i, a in enumerate("ARNDCQEGHILKMFPSTWYV")}
UNK_IDX = 20


def parse_pdb_complex(pdb_path: str | Path, query_chain: str | None = None) -> dict | None:
    """Parse a PDB file and return the per-complex dict (without ddg labels).

    Uses Biopython for parsing. Skips residues missing any of N/CA/C/O.
    Returns None if the structure has fewer than 10 residues or only one chain.
    """
    try:
        from Bio.PDB import PDBParser
    except ImportError as e:
        raise RuntimeError("biopython is required to parse Dockground PDBs") from e

    parser = PDBParser(QUIET=True)
    pdb_path = Path(pdb_path)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(structure.get_models())

    coords_list, res_types, chain_ids, seq_idx_list, pdb_resnums = [], [], [], [], []
    chain_id_to_int: dict[str, int] = {}
    chain_seq_counter: dict[str, int] = {}

    for chain in model:
        cid = chain.id
        if cid == " " or cid == "":
            continue
        # Register chain lazily — skip chains that contain only HETATMs
        # (e.g. chain C in 1nam is all NAG sugars). Otherwise they'd appear
        # in meta["chains"] but have zero residues, breaking downstream
        # FoldX --analyseComplexChains arg construction.
        cint = None

        for res in chain:
            if res.id[0] != " ":  # skip HETATM/water
                continue
            try:
                n = res["N"].coord
                ca = res["CA"].coord
                c = res["C"].coord
                o = res["O"].coord
            except KeyError:
                continue
            if cint is None:
                if cid not in chain_id_to_int:
                    chain_id_to_int[cid] = len(chain_id_to_int)
                    chain_seq_counter[cid] = 0
                cint = chain_id_to_int[cid]
            three = res.get_resname().upper()
            one = THREE_TO_ONE.get(three, "X")
            res_types.append(AA_TO_IDX.get(one, UNK_IDX))
            coords_list.append(np.stack([n, ca, c, o], axis=0))
            chain_ids.append(cint)
            seq_idx_list.append(chain_seq_counter[cid])
            # Preserve the actual PDB residue number — needed by Rosetta's
            # pdb2pose() and FoldX's mutation strings. Insertion codes are
            # ignored (rare in PPBS); res.id is (hetflag, resnum, icode).
            pdb_resnums.append(int(res.id[1]))
            chain_seq_counter[cid] += 1

    # Inference is apo / query-chain-only, so a single-chain (monomer) input is
    # valid here — unlike the training-data parser this is derived from, we do
    # NOT require >=2 chains. Only reject empty / near-empty structures.
    if len(coords_list) < 10 or len(chain_id_to_int) < 1:
        return None

    coords = np.stack(coords_list, axis=0).astype(np.float32)
    res_types_arr = np.array(res_types, dtype=np.int8)
    chain_ids_arr = np.array(chain_ids, dtype=np.int8)
    seq_idx_arr = np.array(seq_idx_list, dtype=np.int16)
    pdb_resnum_arr = np.array(pdb_resnums, dtype=np.int32)

    if query_chain is None:
        # Default: largest chain is the query.
        unique, counts = np.unique(chain_ids_arr, return_counts=True)
        query_int = int(unique[counts.argmax()])
    else:
        if query_chain not in chain_id_to_int:
            return None
        query_int = chain_id_to_int[query_chain]
    is_query = (chain_ids_arr == query_int).astype(bool)

    return {
        "coords": coords,
        "res_types": res_types_arr,
        "chain_ids": chain_ids_arr,
        "seq_idx": seq_idx_arr,
        "pdb_resnum": pdb_resnum_arr,
        "is_query": is_query,
        "meta": {
            "pdb_id": pdb_path.stem,
            "chains": list(chain_id_to_int.keys()),
            "query_chain": [c for c, i in chain_id_to_int.items() if i == query_int][0],
        },
    }


def iter_dockground_pdbs(raw_dir: str | Path) -> Iterable[Path]:
    raw_dir = Path(raw_dir)
    yield from sorted(raw_dir.glob("*.pdb"))
    yield from sorted(raw_dir.glob("*.ent"))
