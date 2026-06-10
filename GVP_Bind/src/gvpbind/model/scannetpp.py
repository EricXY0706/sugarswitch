"""ScanNet++ — per-residue ddG regressor over a multi-chain residue graph.

Architecture:
    residue features (s, V) ── 6× GVPConv ── per-residue scalar head ── ddG

Node scalars: residue-type embedding [20+1] + chain embedding [8] + is_query flag
Node vectors: 3 (N->Ca, Ca->C, sidechain direction proxy)
Edge scalars: distance RBF (16) + sequence-separation embedding (8) + same-chain flag
Edge vectors: 1 unit displacement vector

The model takes a *batched flat* node tensor (all complexes concatenated)
plus a `batch` index. Edge_index is global (already offset by collate_fn).
"""
from __future__ import annotations

from typing import Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from .gvp_layers import GVP, GVPConv, GVPLayerNorm
from .atom_encoder import AtomEncoder

ScalarVector = Tuple[torch.Tensor, torch.Tensor]

# 20 standard AAs + UNK
N_AA = 21
N_CHAIN_EMB = 8
SEQ_SEP_BUCKETS = 8


def _rbf(d: torch.Tensor, d_min: float = 0.0, d_max: float = 20.0, n: int = 16) -> torch.Tensor:
    """Gaussian radial basis features for distances [..., 1] -> [..., n]."""
    centers = torch.linspace(d_min, d_max, n, device=d.device)
    width = (d_max - d_min) / n
    return torch.exp(-((d.unsqueeze(-1) - centers) / width) ** 2)


def _seq_sep_bucket(diff: torch.Tensor, n_buckets: int = SEQ_SEP_BUCKETS) -> torch.Tensor:
    """Bucket signed sequence separation into log-spaced bins; -1 for inter-chain."""
    abs_d = torch.abs(diff).clamp(min=1).float()
    # log2(1)=0 .. log2(>=64)=6 -> 0..6 buckets, last bucket = inter-chain
    bucket = torch.floor(torch.log2(abs_d)).clamp(max=n_buckets - 2).long()
    return bucket


class ScanNetPP(nn.Module):
    def __init__(
        self,
        node_h_dims: Tuple[int, int] = (128, 16),
        edge_h_dims: Tuple[int, int] = (32, 1),
        n_layers: int = 6,
        drop_rate: float = 0.1,
        n_chain_emb: int = N_CHAIN_EMB,
        ddg_out_dim: int = 1,
        n_task_layers: int = 0,
        esm_dim: int = 0,
        esm_proj_dim: int = 64,
        pssm_dim: int = 0,
        pssm_proj_dim: int = 16,
        sidechain: bool = False,
        atomgraph: bool = False,
        atom_layers: int = 2,
        contrastive: bool = False,
        contrast_dim: int = 64,
    ) -> None:
        super().__init__()
        self.node_h_dims = node_h_dims
        self.edge_h_dims = edge_h_dims
        self.ddg_out_dim = ddg_out_dim
        self.n_task_layers = n_task_layers
        self.esm_dim = esm_dim
        self.esm_proj_dim = esm_proj_dim if esm_dim > 0 else 0
        self.pssm_dim = pssm_dim
        self.pssm_proj_dim = pssm_proj_dim if pssm_dim > 0 else 0
        self.sidechain = sidechain
        self.atomgraph = atomgraph
        if atomgraph:
            self.atom_enc = AtomEncoder(n_layers=atom_layers)
        atom_s = self.atom_enc.out_dims[0] if atomgraph else 0
        atom_v = self.atom_enc.out_dims[1] if atomgraph else 0

        # Contrastive (EGCPPIS-style) alignment of the atom-tower per-residue
        # representation with the residue-tower trunk output. Both are projected
        # into a shared unit-sphere and pulled together (same residue) / pushed
        # apart (other residues) via InfoNCE in the LightningModule. Requires an
        # atom tower to align against.
        self.contrastive = contrastive
        self.contrast_dim = contrast_dim if contrastive else 0
        if contrastive:
            assert atomgraph, "contrastive=True requires atomgraph=True (need an atom tower to align)"
            self.proj_atom = nn.Linear(atom_s, contrast_dim)
            self.proj_res = nn.Linear(node_h_dims[0], contrast_dim)

        # Node input scalar features:
        #   res_type_emb (16) + chain_emb (8) + is_query (1) = 25
        # Node input vector features: 3 (N->Ca, Ca->C, backbone "out") [+2 side-chain] [+atom_v]
        self.res_emb = nn.Embedding(N_AA, 16)
        self.chain_emb = nn.Embedding(n_chain_emb, 8)
        # side-chain adds 3 scalars (count, extent, burial) + 2 vector channels;
        # atom encoder adds (atom_s scalars, atom_v vector channels) per residue.
        node_in_s = 16 + 8 + 1 + self.esm_proj_dim + self.pssm_proj_dim + (3 if sidechain else 0) + atom_s
        node_in_v = 3 + (2 if sidechain else 0) + atom_v
        self.node_in = GVP((node_in_s, node_in_v), node_h_dims)
        if esm_dim > 0:
            self.esm_proj = nn.Linear(esm_dim, esm_proj_dim)
        if pssm_dim > 0:
            self.pssm_proj = nn.Linear(pssm_dim, pssm_proj_dim)

        # Edge input scalar features:
        #   rbf (16) + seq_sep (8) + same_chain (1) = 25
        # Edge input vector features: 1 unit displacement
        self.seq_sep_emb = nn.Embedding(SEQ_SEP_BUCKETS, 8)
        edge_in_s = 16 + 8 + 1
        edge_in_v = 1
        self.edge_in = GVP((edge_in_s, edge_in_v), edge_h_dims)

        # n_layers includes the per-task layers: the shared trunk gets
        # (n_layers - n_task_layers) layers, and each task gets n_task_layers
        # branched layers on top. This keeps total per-token compute roughly
        # constant when n_task_layers > 0 vs the baseline 6-layer shared trunk.
        n_shared = max(1, n_layers - n_task_layers)
        self.layers = nn.ModuleList(
            [GVPConv(node_h_dims, edge_h_dims, drop_rate=drop_rate) for _ in range(n_shared)]
        )
        # Per-task branched GVPConv blocks (empty ModuleList when n_task_layers == 0).
        self.task_layers_ddg = nn.ModuleList(
            [GVPConv(node_h_dims, edge_h_dims, drop_rate=drop_rate) for _ in range(n_task_layers)]
        )
        self.task_layers_binary = nn.ModuleList(
            [GVPConv(node_h_dims, edge_h_dims, drop_rate=drop_rate) for _ in range(n_task_layers)]
        )

        self.head = nn.Sequential(
            GVP(node_h_dims, (node_h_dims[0], 0)),
            nn.Linear(node_h_dims[0], ddg_out_dim),
        )
        # Parallel head for binary interface (PPBS) classification.
        self.head_binary = nn.Sequential(
            GVP(node_h_dims, (node_h_dims[0], 0)),
            nn.Linear(node_h_dims[0], 1),
        )
        # Pairwise cross-chain contact head: scores (i, j) residue pairs from
        # per-residue scalar embeddings. Used only by forward_pairwise(); the
        # backbone is fed an intra-chain-only edge_index by the dataset, so
        # node_x already has no relative-pose signal between chains.
        self.head_pairwise = nn.Sequential(
            nn.Linear(3 * node_h_dims[0], node_h_dims[0]),
            nn.ReLU(),
            nn.Dropout(drop_rate),
            nn.Linear(node_h_dims[0], 1),
        )

    @staticmethod
    def build_node_vectors(coords: torch.Tensor) -> torch.Tensor:
        """coords: [N, 4, 3] for (N, CA, C, O). Returns [N, 3, 3] vector channels:
        N->CA, CA->C, sidechain direction proxy (CA - centroid(N,C))."""
        N_atom, CA, C = coords[:, 0], coords[:, 1], coords[:, 2]
        v1 = N_atom - CA
        v2 = C - CA
        v3 = CA - 0.5 * (N_atom + C)  # local backbone "out" vector; proxy for sidechain orientation
        # Normalize so they're roughly unit-length; preserves direction equivariance.
        def _u(v):
            n = torch.linalg.vector_norm(v, dim=-1, keepdim=True).clamp(min=1e-6)
            return v / n
        return torch.stack([_u(v1), _u(v2), _u(v3)], dim=-2)  # [N, 3, 3]

    def featurize(
        self,
        coords: torch.Tensor,         # [N, 4, 3]
        res_types: torch.Tensor,      # [N] long, in [0, 21)
        chain_ids: torch.Tensor,      # [N] long, in [0, n_chain_emb)
        is_query: torch.Tensor,       # [N] bool/float
        edge_index: torch.Tensor,     # [2, E] long (global, post-batch offset)
        seq_idx: torch.Tensor,        # [N] long, residue sequence position within its chain
        esm_features: torch.Tensor | None = None,  # [N, esm_dim] optional pre-computed ESM-2 per-residue
        pssm_features: torch.Tensor | None = None,  # [N, pssm_dim] optional MSA PWM per-residue
        sidechain_features: torch.Tensor | None = None,  # [N, 9]: d_cb(3) d_centroid(3) count/extent/burial(3)
        atom: dict | None = None,  # atomgraph inputs: xyz/elem/scflag/edge_index/resid
    ) -> Tuple[ScalarVector, ScalarVector, torch.Tensor | None]:
        # Bio-assemblies may have >n_chain_emb chains; collapse the overflow
        # into the last bucket. is_query already carries the query/partner
        # distinction we actually need.
        chain_ids = chain_ids.clamp(min=0, max=self.chain_emb.num_embeddings - 1)
        # Node features
        parts = [self.res_emb(res_types), self.chain_emb(chain_ids), is_query.float().unsqueeze(-1)]
        if self.esm_dim > 0:
            assert esm_features is not None, "model built with esm_dim>0 but no esm_features passed in"
            parts.append(self.esm_proj(esm_features.to(parts[0].dtype)))
        if self.pssm_dim > 0:
            assert pssm_features is not None, "model built with pssm_dim>0 but no pssm_features passed in"
            parts.append(self.pssm_proj(pssm_features.to(parts[0].dtype)))
        if self.sidechain:
            assert sidechain_features is not None, "model built with sidechain=True but no sidechain_features passed in"
            parts.append(sidechain_features[:, 6:9].to(parts[0].dtype))   # count, extent, burial
        atom_s = atom_v = None
        if self.atomgraph:
            assert atom is not None, "model built with atomgraph=True but no atom inputs passed in"
            atom_s, atom_v = self.atom_enc(atom["xyz"], atom["elem"], atom["scflag"],
                                           atom["edge_index"], atom["resid"], coords.shape[0])
            parts.append(atom_s.to(parts[0].dtype))
        s_node = torch.cat(parts, dim=-1)
        V_node = self.build_node_vectors(coords)
        if self.sidechain:
            def _u(v):
                return v / torch.linalg.vector_norm(v, dim=-1, keepdim=True).clamp(min=1e-6)
            sc_v = torch.stack([_u(sidechain_features[:, 0:3]), _u(sidechain_features[:, 3:6])], dim=-2)  # [N,2,3]
            V_node = torch.cat([V_node, sc_v.to(V_node.dtype)], dim=-2)
        if self.atomgraph:
            V_node = torch.cat([V_node, atom_v.to(V_node.dtype)], dim=-2)
        node_x = self.node_in((s_node, V_node))

        # Edge features
        ca = coords[:, 1]  # [N, 3]
        src, dst = edge_index[0], edge_index[1]
        disp = ca[dst] - ca[src]  # [E, 3]
        d = torch.linalg.vector_norm(disp, dim=-1)  # [E]
        rbf = _rbf(d)  # [E, 16]
        same_chain = (chain_ids[src] == chain_ids[dst]).float().unsqueeze(-1)  # [E, 1]
        seq_diff = seq_idx[dst] - seq_idx[src]
        # Inter-chain edges get the last bucket regardless of seq_diff.
        ss_bucket = _seq_sep_bucket(seq_diff)
        ss_bucket = torch.where(same_chain.squeeze(-1).bool(), ss_bucket, torch.full_like(ss_bucket, SEQ_SEP_BUCKETS - 1))
        s_edge = torch.cat([rbf, self.seq_sep_emb(ss_bucket), same_chain], dim=-1)
        V_edge = (disp / d.clamp(min=1e-6).unsqueeze(-1)).unsqueeze(-2)  # [E, 1, 3]
        edge_x = self.edge_in((s_edge, V_edge))

        return node_x, edge_x, atom_s

    def _apply_head(self, head: nn.Sequential, node_x: ScalarVector) -> torch.Tensor:
        head_out = node_x
        for m in head:
            if isinstance(m, GVP):
                head_out = m(head_out)
            else:
                head_out = m(head_out if not isinstance(head_out, tuple) else head_out[0])
        # Squeeze trailing dim only when it is 1 (scalar regression head).
        # Keep [N, K] for multi-class ddg bin heads.
        return head_out.squeeze(-1) if head_out.shape[-1] == 1 else head_out

    def forward(
        self,
        coords: torch.Tensor,
        res_types: torch.Tensor,
        chain_ids: torch.Tensor,
        is_query: torch.Tensor,
        edge_index: torch.Tensor,
        seq_idx: torch.Tensor,
        return_dict: bool = False,
        esm_features: torch.Tensor | None = None,
        pssm_features: torch.Tensor | None = None,
        sidechain_features: torch.Tensor | None = None,
        atom: dict | None = None,
    ) -> torch.Tensor | dict:
        """Returns per-residue ddG predictions [N] (default), or
        {"ddg": [N], "binary_logit": [N]} if return_dict=True."""
        node_x, edge_x, atom_s = self.featurize(coords, res_types, chain_ids, is_query, edge_index, seq_idx,
                                                esm_features=esm_features, pssm_features=pssm_features,
                                                sidechain_features=sidechain_features, atom=atom)
        for layer in self.layers:
            node_x = layer(node_x, edge_index, edge_x)
        # Per-task branched layers (no-ops when n_task_layers == 0).
        node_x_ddg = node_x
        for layer in self.task_layers_ddg:
            node_x_ddg = layer(node_x_ddg, edge_index, edge_x)
        ddg = self._apply_head(self.head, node_x_ddg)
        if not return_dict:
            return ddg
        node_x_bin = node_x
        for layer in self.task_layers_binary:
            node_x_bin = layer(node_x_bin, edge_index, edge_x)
        binary_logit = self._apply_head(self.head_binary, node_x_bin)
        out = {"ddg": ddg, "binary_logit": binary_logit}
        # Contrastive alignment projections (atom tower vs residue trunk), aligned
        # per-residue. node_x[0] is the shared-trunk scalar embedding before task
        # branches; atom_s is the pre-fusion atom-pooled scalar embedding.
        if self.contrastive and atom_s is not None:
            out["z_atom"] = self.proj_atom(atom_s.to(node_x[0].dtype))
            out["z_res"] = self.proj_res(node_x[0])
        return out

    def forward_pairwise(
        self,
        coords: torch.Tensor,
        res_types: torch.Tensor,
        chain_ids: torch.Tensor,
        is_query: torch.Tensor,
        edge_index: torch.Tensor,
        seq_idx: torch.Tensor,
        pair_index: torch.Tensor,    # [2, P] long, batched node indices for (i, j) pairs
        esm_features: torch.Tensor | None = None,
        pssm_features: torch.Tensor | None = None,
        sidechain_features: torch.Tensor | None = None,
        atom: dict | None = None,
    ) -> torch.Tensor:
        """Returns per-pair logits [P] for cross-chain contact prediction.

        `edge_index` is expected to contain intra-chain edges only — the dataset
        filters cross-chain edges in pairwise_mode so the backbone never sees
        relative L↔R pose.
        """
        node_x, edge_x, _ = self.featurize(coords, res_types, chain_ids, is_query, edge_index, seq_idx,
                                           esm_features=esm_features, pssm_features=pssm_features,
                                           sidechain_features=sidechain_features, atom=atom)
        for layer in self.layers:
            node_x = layer(node_x, edge_index, edge_x)
        h = node_x[0]  # scalar channel: [N, s_dim]
        h_i = h[pair_index[0]]
        h_j = h[pair_index[1]]
        feat = torch.cat([h_i, h_j, h_i * h_j], dim=-1)  # [P, 3*s_dim]
        return self.head_pairwise(feat).squeeze(-1)
