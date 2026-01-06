'''
找到位点-->
定位空间位置-->
根据输入个数均匀采样位点-->
以采样位点为中心，前后各2/3个氨基酸共5/7个位点重设计-->
NXS/T sequon滑窗放置，其他位置为X，结构预测-->
ligandmpnn设计序列，设置--redesigned_residues，设置NXS/T中的X为P-negative偏置，剩余区域无偏置-->
结构预测，和WT结构对比，设置filter
'''

