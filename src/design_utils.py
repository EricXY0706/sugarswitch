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
    conservation_df: pd.DataFrame,
    coupling_stength: np.array,
    interaction_dict: dict,
    rsasa_index_dict: dict,
    chain_id: str,
    wt_seq: str,
    num_sites_per_comb: int = 3,
    combination_id: int = 0,
    return_num_combinations: bool = False,
):
    df = pd.read_csv(scoring_df)

    sites_states = {
        s: report_four_sites_state(
            chain_id=chain_id,
            center_site=s,
            wt_seq=wt_seq,
            conservation_df=conservation_df,
            coupling_stength=coupling_stength,
            interaction_dict=interaction_dict,
            rsasa_index_dict=rsasa_index_dict,
        )
        for s in df["Site"].tolist()
    }

    GLY_MOTIFS_RANK = ["*YYY", "*YNY", "*YYN", "*YNN", "*Y*V", "*YVV"]

    def _match_state_pattern(state: str, pattern: str) -> bool:
        if len(state) != len(pattern):
            return False
        return all(p == "*" or p == c for c, p in zip(state, pattern))

    def _state_rank(state: str) -> int:
        for rank, pattern in enumerate(GLY_MOTIFS_RANK):
            if _match_state_pattern(state, pattern):
                return rank
        return len(GLY_MOTIFS_RANK)

    borda_scores = {
        int(row["Site"]): float(row["Borda_score"])
        for _, row in df.iterrows()
    }

    # Unit-site priority:
    # 1. GLY_MOTIFS_RANK order
    # 2. Higher Borda_score within the same state class
    # 3. Smaller site index as a stable tie-breaker
    ranked_sites = sorted(
        df["Site"].astype(int).tolist(),
        key=lambda s: (
            _state_rank(sites_states[s]),
            -borda_scores[s],
            s,
        ),
    )
    site_priority = {site: rank for rank, site in enumerate(ranked_sites)}

    editable_sites = np.array(sorted(df["Site"].astype(int).tolist()))
    ca_coords = StructureLoader.get_ca_coords(
        structure_file=structure_file,
        chain_id=chain_id,
    )[editable_sites - 1]

    contact_1d = np.abs(editable_sites[:, None] - editable_sites[None, :])
    contact_3d = np.linalg.norm(
        ca_coords[:, None, :] - ca_coords[None, :, :],
        axis=-1,
    )

    # dist_1d = np.max(contact_1d) / num_sites_per_comb
    dist_1d = 10.
    dist_3d = np.max(contact_3d) / num_sites_per_comb

    contact_map = (contact_1d >= dist_1d) & (contact_3d >= dist_3d)
    sites_combs = find_cliques(contact_map, editable_sites, num_sites_per_comb)

    if not sites_combs:
        raise ValueError("No valid site combinations found under distance constraints.")

    def _comb_priority(comb):
        priorities = sorted(site_priority[int(s)] for s in comb)
        sites = tuple(sorted(int(s) for s in comb))
        return (
            tuple(priorities),
            sum(priorities),
            sites,
        )

    sites_combs = sorted(sites_combs, key=_comb_priority)

    sampled_sites = sites_combs[min(combination_id, len(sites_combs) - 1)]

    if return_num_combinations:
        return sampled_sites, len(sites_combs)
    return sampled_sites

def report_site_state(
    chain_id: str,
    site: int,
    wt_seq: str,
    conservation_df: pd.DataFrame,
    coupling_stength: np.array,
    interaction_dict: dict,
    rsasa_index_dict: dict,
):
    if site < 1 or site > len(wt_seq):
        return "V"
    else:
        if conservation_df is not None:
            aa = conservation_df.loc[conservation_df["i"] == site, "A_i"].values[0]
            freq = conservation_df.loc[conservation_df["i"] == site, aa].values[0] <= 0.8
        else:
            freq = True
        if coupling_stength is not None:
            co_evo = coupling_stength[site-1] <= 1.0
        else:
            co_evo = True
        
        interaction_num = len(interaction_dict.get(f"{chain_id}:{site}:{ONE_TO_THREE[list(wt_seq)[site-1]]}", [])) <= 1
        rsasa = rsasa_index_dict[site] >= 0.2
        
        if freq and co_evo and interaction_num and rsasa:
            return "Y"
        else:
            return "N"

def report_four_sites_state(
    chain_id: str,
    center_site: int,
    wt_seq: str,
    conservation_df: pd.DataFrame,
    coupling_stength: np.array,
    interaction_dict: dict,
    rsasa_index_dict: dict,
):
    state = ""
    for s in [center_site - 1, center_site, center_site + 1, center_site + 2]:
        if s == center_site - 1:
            if wt_seq[s-1] == "P":
                state += "Y"
            else:
                state += "N"
        else:
            state += report_site_state(
                chain_id=chain_id,
                site=s, 
                wt_seq=wt_seq, 
                conservation_df=conservation_df,
                coupling_stength=coupling_stength,
                interaction_dict=interaction_dict,
                rsasa_index_dict=rsasa_index_dict,
            )
    return state
    
def add_single_mutation(
    chain_id: str,
    site: int,
    wt_seq: str,
    conservation_df: pd.DataFrame,
    coupling_stength: np.array,
    interaction_dict: dict,
    rsasa_index_dict: dict,
    N_minus_1_aa: str = "X",
    N_plus_1_aa: str = "X",
    N_plus_2_aa: str = "X",
    trial_times: list = [0],
    gly_site_seqid: int = 0,
):
    
    GLY_MOTIFS = {
        "YYYY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",              
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXNXT", f"PXNXS",
                 f"PXNXT{N_plus_2_aa}", f"PXNXS{N_plus_2_aa}",
                 f"PXN{N_plus_1_aa}T{N_plus_2_aa}", f"PXN{N_plus_1_aa}S{N_plus_2_aa}",
                 f"PXNXT{N_plus_1_aa}{N_plus_2_aa}", f"PXNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "NYYY": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}NXT", f"{N_minus_1_aa}NXS",
                 f"{N_minus_1_aa}NXT{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_2_aa}",
                 f"{N_minus_1_aa}N{N_plus_1_aa}T{N_plus_2_aa}", f"{N_minus_1_aa}N{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "YYYN": [f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"PXNXT{N_plus_2_aa}", f"PXNXS{N_plus_2_aa}",
                 f"PXN{N_plus_1_aa}T{N_plus_2_aa}", f"PXN{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"PXNXT{N_plus_1_aa}{N_plus_2_aa}", f"PXNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "NYYN": [f"{N_minus_1_aa}NXT{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_2_aa}",
                 f"{N_minus_1_aa}N{N_plus_1_aa}T{N_plus_2_aa}", f"{N_minus_1_aa}N{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "YYNY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXN{N_plus_1_aa}T{N_plus_2_aa}", f"PXN{N_plus_1_aa}S{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"PXNXT{N_plus_1_aa}{N_plus_2_aa}", f"PXNXS{N_plus_1_aa}{N_plus_2_aa}"
                 ],
        "NYNY": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}N{N_plus_1_aa}T{N_plus_2_aa}", f"{N_minus_1_aa}N{N_plus_1_aa}S{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}{N_plus_2_aa}", f"{N_minus_1_aa}NXS{N_plus_1_aa}{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}", 
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
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",              
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 ],
        "VYYN": [f"XNXT{N_plus_2_aa}", f"XNXS{N_plus_2_aa}",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}", 
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
                 ],
        "VYNY": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}T{N_plus_2_aa}", f"XN{N_plus_1_aa}S{N_plus_2_aa}",
                 f"XNXT{N_plus_1_aa}{N_plus_2_aa}", f"XNXS{N_plus_1_aa}{N_plus_2_aa}",
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
                 f"XN{N_plus_1_aa}TX", f"XN{N_plus_1_aa}SX",              
                 f"XNXT{N_plus_1_aa}X", f"XNXS{N_plus_1_aa}X",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXNXT", f"PXNXS",
                 f"PXNXTX", f"PXNXSX",
                 f"PXN{N_plus_1_aa}TX", f"PXN{N_plus_1_aa}SX",
                 f"PXNXT{N_plus_1_aa}X", f"PXNXS{N_plus_1_aa}X"
                 ],
        "NYYV": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}NXT", f"{N_minus_1_aa}NXS",
                 f"{N_minus_1_aa}NXTX", f"{N_minus_1_aa}NXSX",
                 f"{N_minus_1_aa}N{N_plus_1_aa}TX", f"{N_minus_1_aa}N{N_plus_1_aa}SX", 
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}X", f"{N_minus_1_aa}NXS{N_plus_1_aa}X",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XNXT", f"XNXS",
                 f"XNXTX", f"XNXSX",
                 f"XN{N_plus_1_aa}TX", f"XN{N_plus_1_aa}SX", 
                 f"XNXT{N_plus_1_aa}X", f"XNXS{N_plus_1_aa}X"
                 ],
        "YYNV": [f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}TX", f"XN{N_plus_1_aa}SX",
                 f"PXN{N_plus_1_aa}T", f"PXN{N_plus_1_aa}S",
                 f"PXN{N_plus_1_aa}TX", f"PXN{N_plus_1_aa}SX",
                 f"XNXT{N_plus_1_aa}X", f"XNXS{N_plus_1_aa}X",
                 f"PXNXT{N_plus_1_aa}X", f"PXNXS{N_plus_1_aa}X"
                 ],
        "NYNV": [f"{N_minus_1_aa}N{N_plus_1_aa}T", f"{N_minus_1_aa}N{N_plus_1_aa}S",
                 f"{N_minus_1_aa}N{N_plus_1_aa}TX", f"{N_minus_1_aa}N{N_plus_1_aa}SX",
                 f"XN{N_plus_1_aa}T", f"XN{N_plus_1_aa}S",
                 f"XN{N_plus_1_aa}TX", f"XN{N_plus_1_aa}SX",
                 f"{N_minus_1_aa}NXT{N_plus_1_aa}X", f"{N_minus_1_aa}NXS{N_plus_1_aa}X",
                 f"XNXT{N_plus_1_aa}X", f"XNXS{N_plus_1_aa}X", 
                 ],
    }
    
    state = report_four_sites_state(
        chain_id=chain_id,
        center_site=site,
        wt_seq=wt_seq,
        conservation_df=conservation_df,
        coupling_stength=coupling_stength,
        interaction_dict=interaction_dict,
        rsasa_index_dict=rsasa_index_dict,
    )

    motifs = GLY_MOTIFS[state]
    motif = motifs[min(len(motifs) - 1, trial_times[gly_site_seqid])]
    n_offset = (site - 1 - max(site - 2, 0)) + (1 if motif.startswith("PX") else 0)

    wt_seq = list(wt_seq)
    mut_seq = "".join(
        wt_seq[:max(site - 2, 0)]
        + list(motif)
        + wt_seq[site + 2:]
    )

    return mut_seq, "".join(wt_seq[max(site - 2, 0):site + 2]), motif, state, n_offset, len(motifs)