'''
找到位点-->
定位空间位置-->
根据输入个数在一维和三维中均匀采样位点-->
以采样位点为中心，前后各2个氨基酸共5个位点重设计，装入5-7个氨基酸（随机）-->
NXS/T sequon滑窗放置，其他位置为X，结构预测-->
ligandmpnn设计序列，设置--redesigned_residues，设置NXS/T中的X为P-negative偏置，剩余区域无偏置-->
结构预测，和WT结构对比，设置filter
恢复糖链结构
'''
import numpy as np

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.util import StructureLoader
def sample_sites(
    structure_file, 
    chain_id: str,
    editable_sites: set[int] = None,
    num_sites: int = 3,

):
    editable_sites = np.array(sorted(list(editable_sites)))
    ca_coords = StructureLoader.get_ca_coords(structure_file=structure_file, chain_id=chain_id)[editable_sites]
    dist_1d = np.abs(editable_sites[:, None] - editable_sites[None, :])
    dist_3d = np.linalg.norm(ca_coords[:, None, :] - ca_coords[None, :, :], axis=-1)

    print(dist_1d)
    print(dist_3d)

    

def boltz_predict(
    sequence,  
):
    pass

def ligandmpnn_redesign(
    
):
    pass



structure_file = '/sdata2/WORK/Xuyi/NB1/NB1.pdb'
chain_id = 'A'
editable_sites = {111, 27, 112, 43, 110, 61, 62, 77, 113, 1, 13, 41, 42, 127, 87, 44, 26, 114, 88, 65, 85, 117, 120, 128, 10, 5, 63, 123, 23, 17, 84, 14, 19, 66, 9, 78, 93, 15, 8}
num_sites = 3
sample_sites(structure_file, chain_id, editable_sites, num_sites)
