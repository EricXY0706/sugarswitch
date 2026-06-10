"""GVP — Geometric Vector Perceptron (Jing et al. 2021, ICLR).

Ported from the reference implementation at
https://github.com/drorlab/gvp-pytorch (MIT). Trimmed to what we need for
ScanNet++ — node-only message passing on a residue graph with scalar +
vector features at each node and edge.

Key idea: every layer takes (s, V) where s is scalar features [N, ds] and
V is vector features [N, dv, 3]. Linear maps on V are applied as right
multiplication of a learned [dv_in, dv_out] matrix on the channel
dimension, leaving the spatial 3-vector axis untouched — this preserves
SO(3) equivariance by construction. Non-linearities act on scalars and on
vector norms only (vectors are then rescaled by gated scalars).
"""
from __future__ import annotations

from typing import Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

ScalarVector = Tuple[torch.Tensor, torch.Tensor]  # (s [..., ds], V [..., dv, 3])


def _norm_no_nan(x: torch.Tensor, axis: int = -1, keepdims: bool = False, eps: float = 1e-8) -> torch.Tensor:
    """Safe vector norm; returns 0 instead of NaN at the origin."""
    return torch.sqrt(torch.sum(x * x, dim=axis, keepdim=keepdims) + eps)


class GVP(nn.Module):
    """Geometric Vector Perceptron.

    in_dims = (ds_in, dv_in), out_dims = (ds_out, dv_out).
    `vector_gate=True` enables the vector-gating variant (Jing 2021 §3.4),
    which gates each output vector channel by a sigmoid of a learned
    scalar — strictly more expressive than the original norm-only variant
    and what we use throughout.
    """

    def __init__(
        self,
        in_dims: Tuple[int, int],
        out_dims: Tuple[int, int],
        h_dim: int | None = None,
        activations: Tuple = (F.relu, torch.sigmoid),
        vector_gate: bool = True,
    ) -> None:
        super().__init__()
        self.si, self.vi = in_dims
        self.so, self.vo = out_dims
        self.vector_gate = vector_gate
        self.scalar_act, self.vector_act = activations

        if self.vi:
            self.h_dim = h_dim if h_dim is not None else max(self.vi, self.vo)
            self.wh = nn.Linear(self.vi, self.h_dim, bias=False)
            self.ws = nn.Linear(self.h_dim + self.si, self.so)
            if self.vo:
                self.wv = nn.Linear(self.h_dim, self.vo, bias=False)
                if self.vector_gate:
                    self.wsv = nn.Linear(self.so, self.vo)
        else:
            self.ws = nn.Linear(self.si, self.so)

    def forward(self, x: ScalarVector | torch.Tensor) -> ScalarVector | torch.Tensor:
        if self.vi:
            s, V = x
            # vector linear: [..., vi, 3] -> [..., h_dim, 3]
            Vh = self.wh(V.transpose(-1, -2)).transpose(-1, -2)
            vn = _norm_no_nan(Vh, axis=-1)  # [..., h_dim]
            s_out = self.ws(torch.cat([s, vn], dim=-1))
            if self.vo:
                Vo = self.wv(Vh.transpose(-1, -2)).transpose(-1, -2)  # [..., vo, 3]
                if self.vector_gate:
                    gate = torch.sigmoid(self.wsv(self.scalar_act(s_out)))
                    Vo = Vo * gate.unsqueeze(-1)
                else:
                    vno = _norm_no_nan(Vo, axis=-1, keepdims=True)
                    Vo = Vo * self.vector_act(vno)
                s_out = self.scalar_act(s_out)
                return s_out, Vo
            return self.scalar_act(s_out)
        else:
            s_out = self.ws(x)
            if self.vo:
                # No input vectors -> output vectors are zero. Caller responsibility.
                raise ValueError("GVP with vi=0 cannot output vector features.")
            return self.scalar_act(s_out)


class GVPDropout(nn.Module):
    """Dropout for (s, V): scalar dropout + per-channel vector dropout."""

    def __init__(self, p: float) -> None:
        super().__init__()
        self.scalar_dropout = nn.Dropout(p)
        self.p = p

    def forward(self, x: ScalarVector) -> ScalarVector:
        s, V = x
        if not self.training or self.p == 0:
            return s, V
        s_out = self.scalar_dropout(s)
        # Same mask across the spatial 3-axis so a dropped channel stays equivariant.
        mask = torch.bernoulli((1 - self.p) * torch.ones(V.shape[:-1], device=V.device))
        V_out = V * mask.unsqueeze(-1) / (1 - self.p)
        return s_out, V_out


class GVPLayerNorm(nn.Module):
    """LayerNorm over scalars + per-channel vector-norm normalization."""

    def __init__(self, dims: Tuple[int, int]) -> None:
        super().__init__()
        self.s, self.v = dims
        self.scalar_norm = nn.LayerNorm(self.s)

    def forward(self, x: ScalarVector) -> ScalarVector:
        s, V = x
        s_out = self.scalar_norm(s)
        if self.v == 0:
            return s_out, V
        # Normalize each vector channel by the RMS of its 3-vector magnitudes
        # over the batch/node dim — but for residue-graph use we just RMS over
        # the spatial axis to keep things local and translation-equivariant.
        vn = _norm_no_nan(V, axis=-1, keepdims=True)
        rms = torch.sqrt(torch.mean(vn ** 2, dim=-2, keepdim=True) + 1e-8)
        return s_out, V / rms


class GVPConv(nn.Module):
    """Message-passing layer over a residue graph.

    Inputs: per-node (s_n [N, sn], V_n [N, vn, 3]) and per-edge
    (s_e [E, se], V_e [E, ve, 3]) plus edge_index [2, E].

    Each edge (i, j) computes a message via a stack of GVPs on the
    concatenation [s_j, V_j, s_e, V_e]. Messages are summed at the
    destination node. A skip MLP then mixes message + node state.
    """

    def __init__(
        self,
        node_dims: Tuple[int, int],
        edge_dims: Tuple[int, int],
        n_message_layers: int = 3,
        drop_rate: float = 0.1,
    ) -> None:
        super().__init__()
        sn, vn = node_dims
        se, ve = edge_dims
        in_dims = (sn + se, vn + ve)
        msg_layers = []
        for i in range(n_message_layers):
            d_in = in_dims if i == 0 else node_dims
            d_out = node_dims
            msg_layers.append(GVP(d_in, d_out))
        self.message_fn = nn.Sequential(*msg_layers)

        self.norm = nn.ModuleList([GVPLayerNorm(node_dims), GVPLayerNorm(node_dims)])
        self.dropout = GVPDropout(drop_rate)
        self.update_fn = nn.Sequential(
            GVP(node_dims, node_dims),
            GVP(node_dims, node_dims, activations=(F.relu, None), vector_gate=True),
        )

    def _message(self, s_src: torch.Tensor, V_src: torch.Tensor, s_e: torch.Tensor, V_e: torch.Tensor) -> ScalarVector:
        s = torch.cat([s_src, s_e], dim=-1)
        V = torch.cat([V_src, V_e], dim=-2)
        # Sequential expects a single arg; we feed (s, V).
        out = (s, V)
        for layer in self.message_fn:
            out = layer(out)
        return out

    def forward(
        self,
        x: ScalarVector,
        edge_index: torch.Tensor,
        edge_attr: ScalarVector,
    ) -> ScalarVector:
        s, V = x
        s_e, V_e = edge_attr
        src, dst = edge_index[0], edge_index[1]

        msg_s, msg_V = self._message(s[src], V[src], s_e, V_e)

        # Scatter-sum aggregation at destination nodes.
        agg_s = torch.zeros_like(s)
        agg_V = torch.zeros_like(V)
        agg_s.index_add_(0, dst, msg_s)
        agg_V.index_add_(0, dst, msg_V)

        # Residual + norm + update
        s_res, V_res = self.norm[0]((s + self.dropout((agg_s, agg_V))[0], V + self.dropout((agg_s, agg_V))[1]))
        upd_s, upd_V = self.update_fn((s_res, V_res))
        s_out, V_out = self.norm[1]((s_res + upd_s, V_res + upd_V))
        return s_out, V_out
