"""
GVP-Bind: per-residue protein binding-site (interface) prediction.
Credit to: https://github.com/LAJ-THU/GVP-Bind
"""
from __future__ import annotations
# import sys
# from pathlib import Path
# # 把 sugarswitch/ 加入 sys.path，使 GVP_Bind 可被导入
# sys.path.insert(0, str(Path(__file__).parent.parent))
import torch

from GVP_Bind.src.gvpbind.data.dataset import build_knn_graph
from GVP_Bind.src.gvpbind.data.parse import parse_pdb_complex, parse_single_chain
from GVP_Bind.src.gvpbind.task.multitask import MultiTask

def gvpbind_predict(
    structure_file: str,
    chain_id: str,
    mode: str = "Complex",
):
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if mode == "Complex":
        parsed = parse_pdb_complex(structure_file, query_chain=chain_id)
        batch = {
            "coords": torch.from_numpy(parsed["coords"]).to(device),
            "res_types": torch.from_numpy(parsed["res_types"]).long().to(device),
            "chain_ids": torch.from_numpy(parsed["chain_ids"]).long().to(device),
            "is_query": torch.from_numpy(parsed["is_query"]).to(device),
            "edge_index": build_knn_graph(torch.from_numpy(parsed["coords"]).to(device)[:, 1]).to(device),
            "seq_idx": torch.from_numpy(parsed["seq_idx"]).long().to(device),
    }
        pdb_resnum = parsed["pdb_resnum"]
    elif mode == "Monomer":
        coords_np, rt_np, si_np, pdb_resnum = parse_single_chain(structure_file, chain_id)
        batch = {
            "coords": torch.from_numpy(coords_np).to(device),
            "res_types": torch.from_numpy(rt_np).long().to(device),
            "chain_ids": torch.zeros(len(torch.from_numpy(rt_np).long()), dtype=torch.long, device=device),
            "is_query": torch.ones(len(torch.from_numpy(rt_np).long()), dtype=torch.bool, device=device),
            "edge_index": build_knn_graph(torch.from_numpy(coords_np)[:, 1]).to(device),
            "seq_idx": torch.from_numpy(si_np).long().to(device),
    }
    
    model_type = {"Complex": "ba", "Monomer": "au"}
    module = MultiTask.load_from_checkpoint(f"./GVP_Bind/checkpoints/gvpbind_au.ckpt", map_location=device, strict=False)
    module.eval().to(device)
    
    with torch.no_grad():
        out = module(batch)
    binary_logit = out["binary_logit"].cpu()
    prob = torch.sigmoid(binary_logit)

    q = batch["is_query"].cpu().bool()
    q_logit = binary_logit[q]
    ranks = torch.empty_like(q_logit)
    order = q_logit.argsort()
    n = max(1, len(q_logit) - 1)
    ranks[order] = torch.arange(len(q_logit), dtype=ranks.dtype) / n
    # scatter back into a full-length vector aligned with prob's indexing
    prob = prob.clone()
    prob[q] = ranks
    
    prob_by_res = {}
    for i in range(len(batch["res_types"])):
        prob_by_res[int(pdb_resnum[i])] = float(prob[i])
        
    return prob_by_res