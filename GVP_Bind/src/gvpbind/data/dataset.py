"""Dataset over per-complex .pt cache files.

Each cached file is a dict with keys:
    coords      : float32 [N, 4, 3]   N, CA, C, O
    res_types   : int8    [N]         0..20 (X = 20)
    chain_ids   : int8    [N]         0..n_chains-1
    is_query    : bool    [N]         True for query-chain residues to score
    is_interface: bool    [N]         True for interface + 6 A shell residues
    ddg         : float32 [N]         Rosetta/FoldX per-residue ddG; 0 outside ddg_mask
    ddg_mask    : bool    [N]         True iff ddg is a real measurement (not a fill)
    seq_idx     : int16   [N]         residue position within chain (0-indexed)
    meta        : dict                {pdb_id, chains, engine, ...}

The loader builds a k-NN graph from CA coordinates on the fly (cheap, ~1 ms
per complex for k=30), so we don't bake graph structure into the cache —
this lets us re-run with a different cutoff without regenerating labels.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import torch
from torch.utils.data import Dataset


def build_knn_graph(ca_coords: torch.Tensor, k: int = 30, max_radius: float = 20.0) -> torch.Tensor:
    """Returns edge_index [2, E] (src, dst) for k nearest neighbors of each residue
    within max_radius. Self-loops are excluded. Edges are directed (i->j and j->i
    both included if both qualify) — caller can dedupe if symmetric is desired.
    """
    n = ca_coords.shape[0]
    if n <= 1:
        return torch.zeros((2, 0), dtype=torch.long, device=ca_coords.device)
    # Pairwise distances [N, N]
    diff = ca_coords.unsqueeze(0) - ca_coords.unsqueeze(1)
    dists = torch.linalg.vector_norm(diff, dim=-1)
    dists.fill_diagonal_(float("inf"))
    k_eff = min(k, n - 1)
    topk = torch.topk(dists, k=k_eff, dim=-1, largest=False)
    nbr_idx = topk.indices  # [N, k]
    nbr_d = topk.values     # [N, k]
    valid = nbr_d <= max_radius
    src = torch.arange(n, device=ca_coords.device).unsqueeze(-1).expand(-1, k_eff)[valid]
    dst = nbr_idx[valid]
    return torch.stack([src, dst], dim=0).long()


class ComplexDataset(Dataset):
    HEADER = "_index.json"  # {stem: {"split": str, "n": int}}

    def __init__(self, cache_dir: str | Path, pdb_ids: list[str] | None = None, k_neighbors: int = 30,
                 max_radius: float = 20.0, splits: list[str] | None = None,
                 ddg_clip: float | None = None,
                 esm_cache_dir: str | Path | None = None,
                 pssm_cache_dir: str | Path | None = None,
                 sidechain_cache_dir: str | Path | None = None,
                 atomgraph_cache_dir: str | Path | None = None,
                 atom_k: int = 16,
                 drop_partner_chain: bool = False,
                 pairwise_mode: bool = False,
                 pairwise_contact_cutoff: float = 8.0,
                 max_complex_residues: int | None = None,
                 force_ddg_mask_false: bool = False):
        """If `splits` is provided, keeps only caches whose `ppbs_split` field is in
        that list. Requires the cache to have been enriched by add_ppbs_binary_labels.

        `ddg_clip`: if set, any residue with |ddg| > ddg_clip is masked out by
        forcing `ddg_mask = False` (and the ddg value is zeroed). Needed because
        ppbs_rosetta_v2 contains Rosetta failure sentinels up to ~10^7 kcal/mol —
        see reports/2026-05-15_ddg_target_audit.html. Typical value: 15.0.

        `esm_cache_dir`: if set, expects `<esm_cache_dir>/<stem>.pt` to contain a
        per-residue ESM-2 embedding tensor of shape [N, esm_dim] (float16 OK). The
        loader will yield an `esm_features` key with the tensor cast to float32.

        `drop_partner_chain`: if True, at load time keep only residues with
        `is_query=True` and recompute the kNN graph from the remaining residues.
        Required for PINDER pre-training: the cache stores both chains, but the
        binary interface label is fully determined by chain L's geometry plus
        `chain_ids`, so leaving both chains in the input lets the model trivially
        solve the task (val_aucpr→1.0) without learning transferable structure.
        Hiding the partner makes the task analogous to PPBS: predict binding
        sites on an isolated chain.

        `pairwise_mode`: if True, keep both chains but FILTER cross-chain edges
        from edge_index (so the backbone encodes each chain independently — no
        distance-ruler shortcut). Additionally compute pair_index [2, P] (all
        query x partner residue pairs, local indices) and pair_label [P] (CA-CA
        distance < `pairwise_contact_cutoff` Å). Used for the PairwisePPI task.
        Mutually exclusive with drop_partner_chain.
        """
        self.cache_dir = Path(cache_dir)
        self.ddg_clip = ddg_clip
        self.esm_cache_dir = Path(esm_cache_dir) if esm_cache_dir is not None else None
        self.pssm_cache_dir = Path(pssm_cache_dir) if pssm_cache_dir is not None else None
        self.sidechain_cache_dir = Path(sidechain_cache_dir) if sidechain_cache_dir is not None else None
        self.atomgraph_cache_dir = Path(atomgraph_cache_dir) if atomgraph_cache_dir is not None else None
        self.atom_k = atom_k
        self.drop_partner_chain = drop_partner_chain
        self.pairwise_mode = pairwise_mode
        self.pairwise_contact_cutoff = float(pairwise_contact_cutoff)
        self.force_ddg_mask_false = bool(force_ddg_mask_false)
        if pairwise_mode and drop_partner_chain:
            raise ValueError("pairwise_mode and drop_partner_chain are mutually exclusive — pairwise needs both chains.")
        if pdb_ids is None:
            self.files = sorted(p for p in self.cache_dir.glob("*.pt") if not p.stem.startswith("_"))
        else:
            self.files = [self.cache_dir / f"{p}.pt" for p in pdb_ids]
            self.files = [f for f in self.files if f.exists()]
        self.k = k_neighbors
        self.max_radius = max_radius
        index = self._load_or_build_index(self.files)
        if splits is not None:
            keep = set(splits)
            self.files = [f for f in self.files if index.get(f.stem, {}).get("split", "") in keep]
        if max_complex_residues is not None:
            before = len(self.files)
            self.files = [f for f in self.files if int(index[f.stem]["n"]) <= max_complex_residues]
            dropped = before - len(self.files)
            if dropped > 0:
                print(f"  filtered {dropped}/{before} complexes with > {max_complex_residues} residues "
                      f"(kept {len(self.files)})")
        self._lengths = [int(index[f.stem]["n"]) for f in self.files]

    def _load_or_build_index(self, files: list[Path]) -> dict[str, dict]:
        """Load {stem: {split, n}} from `_index.json`, building it on first use.
        Eliminates the per-process pickle scan that otherwise dominates startup."""
        import json
        import pickle
        header = self.cache_dir / self.HEADER
        index: dict[str, dict] = {}
        if header.exists():
            try:
                index = json.loads(header.read_text())
            except Exception:
                index = {}
        missing = [f for f in files if f.stem not in index]
        if not missing:
            return index
        for f in missing:
            with open(f, "rb") as fh:
                d = pickle.load(fh)
            index[f.stem] = {
                "split": d.get("ppbs_split", ""),
                "n": int(d["coords"].shape[0]),
            }
        # Best-effort write; tolerate races between concurrent jobs.
        try:
            header.write_text(json.dumps(index))
        except Exception:
            pass
        return index

    def lengths(self) -> list[int]:
        return list(self._lengths)

    def __len__(self) -> int:
        return len(self.files)

    def __getitem__(self, idx: int) -> dict[str, Any]:
        import pickle
        with open(self.files[idx], "rb") as f:
            d = pickle.load(f)
        coords = torch.as_tensor(d["coords"], dtype=torch.float32)
        res_types = torch.as_tensor(d["res_types"], dtype=torch.long)
        chain_ids = torch.as_tensor(d["chain_ids"], dtype=torch.long)
        is_query = torch.as_tensor(d["is_query"], dtype=torch.bool)
        is_interface = torch.as_tensor(d["is_interface"], dtype=torch.bool)
        ddg_mask = torch.as_tensor(d["ddg_mask"], dtype=torch.bool) if "ddg_mask" in d else is_interface
        ddg = torch.as_tensor(d["ddg"], dtype=torch.float32)
        n_full = coords.shape[0]
        is_interface_ppbs = (
            torch.as_tensor(d["is_interface_ppbs"], dtype=torch.bool)
            if "is_interface_ppbs" in d else torch.zeros(n_full, dtype=torch.bool)
        )
        ppbs_mask = (
            torch.as_tensor(d["ppbs_mask"], dtype=torch.bool)
            if "ppbs_mask" in d else torch.zeros(n_full, dtype=torch.bool)
        )
        seq_idx = torch.as_tensor(d["seq_idx"], dtype=torch.long)
        esm_features = None
        if self.esm_cache_dir is not None:
            esm_path = self.esm_cache_dir / f"{self.files[idx].stem}.pt"
            esm_features = torch.load(esm_path, map_location="cpu", weights_only=True).to(torch.float32)
            assert esm_features.shape[0] == n_full, \
                f"ESM shape mismatch for {self.files[idx].stem}: got {tuple(esm_features.shape)} vs N={n_full}"
        pssm_features = None
        if self.pssm_cache_dir is not None:
            import numpy as np
            pssm_path = self.pssm_cache_dir / f"{self.files[idx].stem}.npy"
            pssm_features = torch.from_numpy(np.load(pssm_path)).to(torch.float32)
            assert pssm_features.shape[0] == n_full, \
                f"PSSM shape mismatch for {self.files[idx].stem}: got {tuple(pssm_features.shape)} vs N={n_full}"
        sidechain_features = None
        if self.sidechain_cache_dir is not None:
            import numpy as np
            sc_path = self.sidechain_cache_dir / f"{self.files[idx].stem}.npy"
            sidechain_features = torch.from_numpy(np.load(sc_path)).to(torch.float32)
            assert sidechain_features.shape[0] == n_full, \
                f"sidechain shape mismatch for {self.files[idx].stem}: got {tuple(sidechain_features.shape)} vs N={n_full}"

        if self.drop_partner_chain:
            keep = is_query
            coords = coords[keep]
            res_types = res_types[keep]
            chain_ids = chain_ids[keep]
            is_query = is_query[keep]
            is_interface = is_interface[keep]
            ddg = ddg[keep]
            ddg_mask = ddg_mask[keep]
            is_interface_ppbs = is_interface_ppbs[keep]
            ppbs_mask = ppbs_mask[keep]
            seq_idx = seq_idx[keep]
            if esm_features is not None:
                esm_features = esm_features[keep]
            if pssm_features is not None:
                pssm_features = pssm_features[keep]
            if sidechain_features is not None:
                sidechain_features = sidechain_features[keep]

        ca = coords[:, 1]
        edge_index = build_knn_graph(ca, k=self.k, max_radius=self.max_radius)
        n = coords.shape[0]

        pair_index = None
        pair_label = None
        if self.pairwise_mode:
            # Drop cross-chain edges so the backbone has no relative-pose signal.
            same_chain_edge = chain_ids[edge_index[0]] == chain_ids[edge_index[1]]
            edge_index = edge_index[:, same_chain_edge]
            q_idx = torch.nonzero(is_query, as_tuple=False).squeeze(-1)
            p_idx = torch.nonzero(~is_query, as_tuple=False).squeeze(-1)
            if q_idx.numel() > 0 and p_idx.numel() > 0:
                # All query x partner pairs; local indices into this complex.
                ii, jj = torch.meshgrid(q_idx, p_idx, indexing="ij")
                pair_index = torch.stack([ii.reshape(-1), jj.reshape(-1)], dim=0)  # [2, P]
                pair_dist = torch.linalg.vector_norm(ca[pair_index[0]] - ca[pair_index[1]], dim=-1)
                pair_label = (pair_dist < self.pairwise_contact_cutoff)
            else:
                pair_index = torch.zeros((2, 0), dtype=torch.long)
                pair_label = torch.zeros((0,), dtype=torch.bool)

        if self.force_ddg_mask_false:
            ddg = torch.zeros_like(ddg)
            ddg_mask = torch.zeros_like(ddg_mask)
        elif self.ddg_clip is not None:
            bad = ddg.abs() > self.ddg_clip
            if bad.any():
                ddg_mask = ddg_mask & ~bad
                ddg = ddg.masked_fill(bad, 0.0)
        ppbs_weight = float(d.get("ppbs_weight", 1.0))
        out = {
            "coords": coords,
            "res_types": res_types,
            "chain_ids": chain_ids,
            "is_query": is_query,
            "is_interface": is_interface,
            "ddg": ddg,
            "ddg_mask": ddg_mask,
            # multi-AA compatibility: PPBS/synthetic items have no mutant target,
            # ddg_mask gates whether this is read. All-zero is the safe default.
            "ddg_target_aa": torch.zeros(ddg.shape[0], dtype=torch.long),
            "is_interface_ppbs": is_interface_ppbs,
            "ppbs_mask": ppbs_mask,
            "ppbs_weight": torch.full((n,), ppbs_weight, dtype=torch.float32),
            "seq_idx": seq_idx,
            "edge_index": edge_index,
            "meta": d.get("meta", {}),
        }
        if esm_features is not None:
            out["esm_features"] = esm_features
        if pssm_features is not None:
            out["pssm_features"] = pssm_features
        if sidechain_features is not None:
            out["sidechain_features"] = sidechain_features
        if self.atomgraph_cache_dir is not None:
            import numpy as np
            z = np.load(self.atomgraph_cache_dir / f"{self.files[idx].stem}.npz")
            axyz = torch.from_numpy(z["xyz"].astype("float32"))
            assert int(z["resid"].max(initial=-1)) < n, \
                f"atomgraph resid out of range for {self.files[idx].stem}: needs drop_partner_chain=ON"
            out["atom_xyz"] = axyz
            out["atom_elem"] = torch.from_numpy(z["elem"].astype("int64"))
            out["atom_scflag"] = torch.from_numpy(z["scflag"].astype("float32"))
            out["atom_resid"] = torch.from_numpy(z["resid"].astype("int64"))
            out["atom_edge_index"] = build_knn_graph(axyz, k=self.atom_k)
        if pair_index is not None:
            out["pair_index"] = pair_index
            out["pair_label"] = pair_label
        return out


class SyntheticDataset(Dataset):
    """Random-but-shape-correct complexes for smoke training without real data."""

    def __init__(self, n_complexes: int = 32, seed: int = 0, n_chains: int = 2):
        self.n_complexes = n_complexes
        self.seed = seed
        self.n_chains = n_chains

    def __len__(self) -> int:
        return self.n_complexes

    def __getitem__(self, idx: int) -> dict[str, Any]:
        g = torch.Generator().manual_seed(self.seed + idx)
        n_per_chain = torch.randint(60, 150, (self.n_chains,), generator=g).tolist()
        n_total = sum(n_per_chain)

        # Random backbone-ish coords; CAs roughly on a chain, partner chain offset.
        coords = torch.randn(n_total, 4, 3, generator=g) * 0.5
        offset = 0
        for c, length in enumerate(n_per_chain):
            base = torch.tensor([c * 15.0, 0.0, 0.0])
            for i in range(length):
                coords[offset + i, 1] += base + torch.tensor([float(i) * 3.8, 0.0, 0.0])
            coords[offset:offset + length, 0] = coords[offset:offset + length, 1] + torch.tensor([-1.0, 0.5, 0.0])
            coords[offset:offset + length, 2] = coords[offset:offset + length, 1] + torch.tensor([1.0, 0.5, 0.0])
            coords[offset:offset + length, 3] = coords[offset:offset + length, 2] + torch.tensor([0.0, 1.0, 0.0])
            offset += length

        chain_ids = torch.cat([torch.full((n,), c, dtype=torch.long) for c, n in enumerate(n_per_chain)])
        seq_idx = torch.cat([torch.arange(n) for n in n_per_chain])
        res_types = torch.randint(0, 20, (n_total,), generator=g)
        is_query = chain_ids == 0  # query chain = 0
        # Pretend the closest 30% of query residues to chain 1 are "interface".
        ca = coords[:, 1]
        partner_mask = chain_ids != 0
        if partner_mask.any():
            min_d_to_partner = torch.cdist(ca[is_query], ca[partner_mask]).min(dim=1).values
            thresh = torch.quantile(min_d_to_partner, 0.3)
            is_interface_q = min_d_to_partner <= thresh
        else:
            is_interface_q = torch.zeros(is_query.sum(), dtype=torch.bool)
        is_interface = torch.zeros(n_total, dtype=torch.bool)
        is_interface[is_query] = is_interface_q

        # Synthetic ddG: positively correlated with depth into the partner cluster.
        ddg = torch.zeros(n_total)
        if is_interface.any():
            d = min_d_to_partner[is_interface_q]
            ddg[is_interface] = (thresh - d) * 0.5 + 0.3 * torch.randn(is_interface.sum(), generator=g)

        # Synthetic engine "scans" every interface residue.
        ddg_mask = is_interface.clone()
        # Synthetic binary supervision: PPBS label tracks is_interface on the query chain.
        is_interface_ppbs = is_interface & is_query
        ppbs_mask = is_query.clone()
        edge_index = build_knn_graph(ca, k=30)
        return {
            "coords": coords,
            "res_types": res_types,
            "chain_ids": chain_ids,
            "is_query": is_query,
            "is_interface": is_interface,
            "ddg": ddg,
            "ddg_mask": ddg_mask,
            # multi-AA compatibility: PPBS/synthetic items have no mutant target,
            # ddg_mask gates whether this is read. All-zero is the safe default.
            "ddg_target_aa": torch.zeros(ddg.shape[0], dtype=torch.long),
            "is_interface_ppbs": is_interface_ppbs,
            "ppbs_mask": ppbs_mask,
            "ppbs_weight": torch.ones(n_total, dtype=torch.float32),
            "seq_idx": seq_idx,
            "edge_index": edge_index,
            "meta": {"pdb_id": f"synth_{idx:04d}"},
        }
