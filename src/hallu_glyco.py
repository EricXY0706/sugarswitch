'''
找到位点--> DONE
定位空间位置--> DONE
根据输入个数在一维和三维中均匀采样位点--> DONE
结合上一步的得分表，找出每个loc得分最高的site--> DONE

以采样位点为中心，前后各2个氨基酸共5个位点重设计，装入5-7个氨基酸（随机）-->
NXS/T sequon滑窗放置，其他位置为X，结构预测-->
ligandmpnn设计序列，设置--redesigned_residues，设置NXS/T中的X为P-negative偏置，剩余区域无偏置-->
结构预测，和WT结构对比，设置filter
恢复糖链结构
'''
import numpy as np
import pandas as pd
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.util import StructureLoader

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
    editable_sites: set[int] = None,
    num_sites_per_comb: int = 3,

):
    editable_sites = np.array(sorted(list(editable_sites)))
    ca_coords = StructureLoader.get_ca_coords(structure_file=structure_file, chain_id=chain_id)[editable_sites - 1]
    contact_1d = np.abs(editable_sites[:, None] - editable_sites[None, :])
    contact_3d = np.linalg.norm(ca_coords[:, None, :] - ca_coords[None, :, :], axis=-1)

    dist_1d = np.max(contact_1d) / num_sites_per_comb
    dist_3d = np.max(contact_3d) / num_sites_per_comb

    contact_map = (contact_1d >= dist_1d) & (contact_3d >= dist_3d)
    sites_combs = find_cliques(contact_map, editable_sites, num_sites_per_comb)
    sites_per_loc = {f"loc{i+1}": set() for i in range(num_sites_per_comb)}
    for sites_comb in sites_combs:
        for i, s in enumerate(sites_comb):
            sites_per_loc[f"loc{i+1}"].add(s)
    
    df = pd.read_csv(scoring_df)
    sampled_sites = []
    for l in sites_per_loc:
        sites = sorted(list(sites_per_loc[l]))
        sampled_sites.append(
            df.loc[df["Site"].isin(sites)].loc[
                lambda x: x["Borda_score"].idxmax()
            ]["Site"]
        )
    
    return sampled_sites


def boltz_predict(
    sequence,  
):
    pass

def ligandmpnn_redesign(
    
):
    pass



structure_file = '/sdata2/WORK/Xuyi/NB1/NB1.pdb'
scoring_df = "/sdata2/WORK/Xuyi/NB1/NB1_single_points.csv"
chain_id = 'A'
editable_sites = {111, 27, 112, 43, 110, 61, 62, 77, 113, 1, 13, 41, 42, 127, 128, 87, 44, 26, 114, 88, 65, 85, 117, 120, 10, 5, 63, 123, 23, 17, 84, 14, 19, 66, 9, 78, 93, 15, 8}
num_sites = 4
sampled_sites = sample_sites(structure_file, scoring_df, chain_id, editable_sites, num_sites)
print(sampled_sites)
