"""Atom-level encoder (option 2 / realization 2): a small GVP message-passing net
over the all-heavy-atom point cloud, mean-pooled atoms -> residue. Produces a
per-residue (scalar, vector) embedding that the residue-level ScanNetPP appends to
its node features — giving the model cross-residue atomic packing / side-chain
shape (the signal PeSTo's all-atom geometry captures, e.g. on test_topology).

Atom inputs are query-chain heavy atoms (see scripts/compute_atomgraph.py); the
atom kNN graph is built per-complex in the dataset and batch-offset in collate.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from .gvp_layers import GVP, GVPConv

N_ELEM = 5  # C, N, O, S, other


def _rbf(d: torch.Tensor, d_max: float = 8.0, n: int = 16) -> torch.Tensor:
    centers = torch.linspace(0.0, d_max, n, device=d.device)
    spacing = centers[1] - centers[0]
    return torch.exp(-(((d.unsqueeze(-1) - centers) / spacing) ** 2))


class AtomEncoder(nn.Module):
    def __init__(self, elem_dim: int = 16, node_dims=(32, 4), edge_dims=(16, 1),
                 n_layers: int = 2, drop_rate: float = 0.1) -> None:
        super().__init__()
        self.out_dims = node_dims
        self.elem_emb = nn.Embedding(N_ELEM, elem_dim)
        # atoms have no intrinsic vector feature: project scalars, init node vectors
        # to zero; GVPConv lifts them to vectors via the (vector-carrying) edges.
        self.node_s_proj = nn.Linear(elem_dim + 1, node_dims[0])
        self.edge_in = GVP((16, 1), edge_dims)                   # rbf(16) + unit disp(1)
        self.layers = nn.ModuleList(
            [GVPConv(node_dims, edge_dims, drop_rate=drop_rate) for _ in range(n_layers)]
        )

    def forward(self, xyz, elem, scflag, edge_index, resid, n_res):
        """xyz[A,3], elem[A], scflag[A], edge_index[2,E] (atom-local, batch-offset),
        resid[A] (global residue index), n_res = total residues in batch.
        Returns (res_s[n_res, hs], res_v[n_res, hv, 3])."""
        s = torch.cat([self.elem_emb(elem.long()), scflag.float().unsqueeze(-1)], dim=-1)
        s0 = self.node_s_proj(s)
        v0 = torch.zeros(s0.shape[0], self.out_dims[1], 3, device=s0.device, dtype=s0.dtype)
        node = (s0, v0)                                          # atoms start with zero vectors
        src, dst = edge_index[0], edge_index[1]
        disp = xyz[dst] - xyz[src]
        d = torch.linalg.vector_norm(disp, dim=-1)
        ev = (disp / d.clamp(min=1e-6).unsqueeze(-1)).unsqueeze(-2)   # [E,1,3]
        edge = self.edge_in((_rbf(d), ev))
        for layer in self.layers:
            node = layer(node, edge_index, edge)
        s_a, v_a = node                                          # [A,hs], [A,hv,3]
        hs, hv = self.out_dims
        resid = resid.long()
        res_s = torch.zeros(n_res, hs, device=s_a.device, dtype=s_a.dtype).index_add_(0, resid, s_a)
        res_v = torch.zeros(n_res, hv, 3, device=v_a.device, dtype=v_a.dtype).index_add_(0, resid, v_a)
        cnt = torch.zeros(n_res, device=s_a.device, dtype=s_a.dtype).index_add_(
            0, resid, torch.ones_like(resid, dtype=s_a.dtype)).clamp(min=1.0)
        return res_s / cnt.unsqueeze(-1), res_v / cnt.unsqueeze(-1).unsqueeze(-1)
