import numpy as np
import pandas as pd
import random
import torch
from src.util import *

GLY_MOTIFS = ["NXTXX", "NXSXX", "XNXTX", "XNXSX", "XXNXT", "XXNXS",
              "NXTXXX", "NXSXXX", "XNXTXX", "XNXSXX", "XXNXTX", "XXNXSX", "XXXNXT", "XXXNXS",
              "NXTXXXX", "NXSXXXX", "XNXTXXX", "XNXSXXX", "XXNXTXX", "XXNXSXX", "XXXNXTX", "XXXNXSX", "XXXXNXT", "XXXXNXS"]

def set_seed(
    seed: int = 42
):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

def find_cliques(contact_map, editable_sites, n):
    N = contact_map.shape[0]
    results = set()

    def backtrack(start, current):

        if len(current) == n:
            sites = sorted(editable_sites[i] for i in current)
            results.add(tuple(sites))
            return

        for next_idx in range(start, N):

            if all(contact_map[next_idx, i] and contact_map[i, next_idx]
                   for i in current):
                backtrack(next_idx + 1, current + [next_idx])

    backtrack(0, [])
    return results

def sample_sites(
    structure_file: str,
    scoring_df: str, 
    chain_id: str,
    num_sites_per_comb: int = 3,
):
    editable_sites = np.array(sorted(pd.read_csv(scoring_df)["Site"].tolist()))
    ca_coords = StructureLoader.get_ca_coords(structure_file=structure_file, chain_id=chain_id)[editable_sites - 1]
    contact_1d = np.abs(editable_sites[:, None] - editable_sites[None, :])
    contact_3d = np.linalg.norm(ca_coords[:, None, :] - ca_coords[None, :, :], axis=-1)

    dist_1d = np.max(contact_1d) / num_sites_per_comb
    dist_3d = np.max(contact_3d) / num_sites_per_comb

    contact_map = (contact_1d >= dist_1d) & (contact_3d >= dist_3d)
    sites_combs = find_cliques(contact_map, editable_sites, num_sites_per_comb)

    df = pd.read_csv(scoring_df)
    sites_scores = {sites: round(np.sum(np.array([df.loc[df["Site"] == s]["SASA"].values[0] for s in sites])), 4) for sites in sites_combs}
    sites_scores = sorted(sites_scores.items(), key=lambda x: x[1], reverse=True)
    sampled_sites = sites_scores[0][0]
    
    return sampled_sites