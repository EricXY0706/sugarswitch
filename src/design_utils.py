import numpy as np
import pandas as pd
import random
import torch
from src.util import *

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
    sites_scores = {sites: round(np.sum(np.array([df.loc[df["Site"] == s]["rSASA"].values[0] for s in sites])), 4) for sites in sites_combs}
    sites_scores = sorted(sites_scores.items(), key=lambda x: x[1], reverse=True)
    sampled_sites = sites_scores[0][0]
    
    return sampled_sites

def report_site_status(
    site: int,
    wt_seq: str,
    structure_unfav_sites: set,
    sequence_unfav_sites: set,
    low_rsasa_sites: set,
    fav_sites: set,
):
    if site < 1 or site > len(wt_seq):
        return "V"
    elif site in fav_sites:
        return "Y"
    elif site in structure_unfav_sites:
        return "N"
    else:
        return "Y"
    
def add_single_mutation(
    site: int,
    wt_seq: str,
    structure_unfav_sites: set,
    sequence_unfav_sites: set,
    low_rsasa_sites: set,
    fav_sites: set,
    N_minus_1_aa: str = "X",
    N_plus_1_aa: str = "X",
    N_plus_2_aa: str = "X",
    trial_times: int = 0,
):
    
    
    GLY_MOTIFS = {
        "YYYY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXNXT", f"PXNXS",
                 f"PXNXT{N_plus_2_aa}", f"PXNXS{N_plus_2_aa}",
                 f"PXNXT{N_plus_1_aa}{N_plus_2_aa}", f"PXNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "NYYY": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}NXT", f"{N_minus_1_aa}NXS",
                 f"{N_minus_1_aa}NXT{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_2_aa}",
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "YYYN": [f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"PXNXT{N_plus_2_aa}", f"PXNXS{N_plus_2_aa}",
                 f"PXNXT{N_plus_1_aa}{N_plus_2_aa}", f"PXNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "NYYN": [f"{N_minus_1_aa}NXT{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_2_aa}",
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "YYNY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXN{N_plus_1_aa}T{N_plus_2_aa}", f"PXN{N_plus_1_aa}S{N_plus_2_aa}"
                 ],
        "NYNY": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}N{N_plus_1_aa}T{N_plus_2_aa}", f"{N_minus_1_aa}N{N_plus_1_aa}S{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}"
                 ],
        "YYNN": [f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"PXNXT{N_plus_1_aa}{N_plus_2_aa}", f"PXNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "NYNN": [f"{N_minus_1_aa}NXT{N_plus_1_aa}{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "VYYY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 ],
        "VYYN": [f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 ],
        "VYNY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",
                 ],
        "VYNN": [f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 ],
        "YYVV": [f"XNXT", f"XNXS",
                 f"XNXTX", f"XNXSX",
                 f"XNXTXX", f"XNXSXX",
                 f"PXNXT", f"PXNXS",
                 f"PXNXTX", f"PXNXSX",
                 f"PXNXTXX", f"PXNXSXX"
                 ],
        "NYVV": [f"{N_minus_1_aa}NXT", f"{N_minus_1_aa}NXS",
                 f"{N_minus_1_aa}NXTX", f"{N_minus_1_aa}NXSX",
                 f"{N_minus_1_aa}NXTXX", f"{N_minus_1_aa}NXSXX",
                 f"XNXT", f"XNXS",
                 f"XNXTX", f"XNXSX",
                 f"XNXTXX", f"XNXSXX"
                 ],
        "YYYV": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXTX", f"XNXSX",
                 f"XNXT{N_plus_1_aa}X", f"XNXS{N_plus_1_aa}X",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXNXT", f"PXNXS",
                 f"PXNXTX", f"PXNXSX",
                 f"PXNXT{N_plus_1_aa}X", f"PXNXS{N_plus_1_aa}X"
                 ],
        "NYYV": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}NXT", f"{N_minus_1_aa}NXS",
                 f"{N_minus_1_aa}NXTX", f"{N_minus_1_aa}NXSX",
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}X", f"{N_minus_1_aa}NXS{N_plus_1_aa}X",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXTX", f"XNXSX",
                 f"XNXT{N_plus_1_aa}X", f"XNXS{N_plus_1_aa}X"
                 ],
        "YYNV": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}TX", f"XN{N_plus_1_aa}SX",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXN{N_plus_1_aa}TX", f"PXN{N_plus_1_aa}SX"
                 ],
        "NYNV": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}N{N_plus_1_aa}TX", f"{N_minus_1_aa}N{N_plus_1_aa}SX",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}TX", f"XN{N_plus_1_aa}SX"
                 ],
    }
    
    status = ""
    for s in [site - 1, site, site + 1, site + 2]:
        status += report_site_status(site=s, wt_seq=wt_seq, structure_unfav_sites=structure_unfav_sites, sequence_unfav_sites=sequence_unfav_sites, low_rsasa_sites=low_rsasa_sites, fav_sites=fav_sites)
    
    wt_seq = list(wt_seq)
    mut_seq = "".join(wt_seq[:site - 2] + list(GLY_MOTIFS[status][min(len(GLY_MOTIFS[status]), trial_times)]) + wt_seq[site + 1:])
    
    return mut_seq