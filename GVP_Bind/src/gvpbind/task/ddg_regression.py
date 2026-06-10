"""Lightning module for per-residue ddG regression."""
from __future__ import annotations

from typing import Any

import pytorch_lightning as pl
import torch
import torch.nn as nn
import torch.nn.functional as F

from ..model.scannetpp import ScanNetPP


def _pearson(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    if x.numel() < 2:
        return torch.tensor(0.0, device=x.device)
    xm, ym = x - x.mean(), y - y.mean()
    denom = (xm.norm() * ym.norm()).clamp(min=1e-8)
    return (xm * ym).sum() / denom


def _spearman(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    if x.numel() < 2:
        return torch.tensor(0.0, device=x.device)
    rx = x.argsort().argsort().float()
    ry = y.argsort().argsort().float()
    return _pearson(rx, ry)


def _aucpr_at_threshold(pred: torch.Tensor, target: torch.Tensor, thresh: float = 1.0) -> torch.Tensor:
    """AUCPR with binary label `target > thresh`. Cheap trapezoidal approximation."""
    if target.numel() < 2:
        return torch.tensor(0.0, device=pred.device)
    binary = (target > thresh).float()
    if binary.sum() == 0 or binary.sum() == binary.numel():
        return torch.tensor(0.0, device=pred.device)
    # Sort predictions descending; sweep precision-recall curve.
    order = pred.argsort(descending=True)
    sorted_b = binary[order]
    tp = torch.cumsum(sorted_b, dim=0)
    fp = torch.cumsum(1 - sorted_b, dim=0)
    precision = tp / (tp + fp).clamp(min=1)
    recall = tp / binary.sum()
    # Trapezoidal area under (recall, precision).
    recall = torch.cat([torch.zeros(1, device=recall.device), recall])
    precision = torch.cat([torch.tensor([1.0], device=precision.device), precision])
    return torch.trapz(precision, recall)


class DDGRegression(pl.LightningModule):
    def __init__(
        self,
        model_cfg: dict | None = None,
        lr: float = 3e-4,
        weight_decay: float = 1e-2,
        warmup_steps: int = 500,
        max_steps: int = 50_000,
        huber_delta: float = 1.0,
        non_interface_weight: float = 0.0,
        binarize_thresh: float = 1.0,
    ) -> None:
        super().__init__()
        self.save_hyperparameters()
        cfg = model_cfg or {}
        self.model = ScanNetPP(**cfg)
        self.lr = lr
        self.weight_decay = weight_decay
        self.warmup_steps = warmup_steps
        self.max_steps_ = max_steps
        self.huber_delta = huber_delta
        self.non_interface_weight = non_interface_weight
        self.binarize_thresh = binarize_thresh

    def forward(self, batch: dict) -> torch.Tensor:
        return self.model(
            coords=batch["coords"],
            res_types=batch["res_types"],
            chain_ids=batch["chain_ids"],
            is_query=batch["is_query"],
            edge_index=batch["edge_index"],
            seq_idx=batch["seq_idx"],
        )

    def _loss(self, pred: torch.Tensor, batch: dict) -> tuple[torch.Tensor, dict]:
        target = batch["ddg"]
        is_query = batch["is_query"]
        is_interface = batch["is_interface"]
        # ddg_mask: True iff target ddg is a real engine measurement.
        # Falls back to is_interface for old caches / synthetic data lacking the field.
        ddg_mask = batch.get("ddg_mask", is_interface)

        # Loss is computed only on query residues.
        q = is_query
        if q.sum() == 0:
            zero = pred.sum() * 0.0
            return zero, {"loss": zero.detach()}

        pred_q = pred[q]
        target_q = target[q]
        interface_q = is_interface[q]
        # Huber is restricted to residues that were actually scanned;
        # A/G/P interface residues with fill ddg=0 contribute no gradient.
        scored_q = ddg_mask[q] & interface_q

        loss_interface = F.huber_loss(pred_q[scored_q], target_q[scored_q],
                                      delta=self.huber_delta, reduction="mean") if scored_q.any() else 0.0 * pred_q.sum()
        # Optional non-interface anchor (default weight 0). When enabled, restrict to
        # truly non-interface residues (their ddg=0 fill is a deliberate prior).
        non_q = ~interface_q
        loss_non = F.mse_loss(pred_q[non_q], target_q[non_q],
                              reduction="mean") if (self.non_interface_weight > 0 and non_q.any()) else 0.0 * pred_q.sum()
        loss = loss_interface + self.non_interface_weight * loss_non

        with torch.no_grad():
            pearson = _pearson(pred_q[scored_q], target_q[scored_q]) if scored_q.any() else torch.tensor(0.0, device=pred.device)
            spearman = _spearman(pred_q[scored_q], target_q[scored_q]) if scored_q.any() else torch.tensor(0.0, device=pred.device)
            aucpr = _aucpr_at_threshold(pred_q[scored_q], target_q[scored_q], self.binarize_thresh) if scored_q.any() else torch.tensor(0.0, device=pred.device)

        metrics = {
            "loss": loss.detach(),
            "loss_interface": (loss_interface.detach() if torch.is_tensor(loss_interface) else torch.tensor(0.0, device=pred.device)),
            "loss_non_interface": (loss_non.detach() if torch.is_tensor(loss_non) else torch.tensor(0.0, device=pred.device)),
            "pearson": pearson,
            "spearman": spearman,
            "aucpr_at1": aucpr,
        }
        return loss, metrics

    def training_step(self, batch: dict, batch_idx: int) -> torch.Tensor:
        pred = self(batch)
        loss, metrics = self._loss(pred, batch)
        for k, v in metrics.items():
            self.log(f"train/{k}", v, on_step=True, on_epoch=True, batch_size=int(batch["is_query"].sum()))
        return loss

    def validation_step(self, batch: dict, batch_idx: int) -> None:
        pred = self(batch)
        _, metrics = self._loss(pred, batch)
        for k, v in metrics.items():
            self.log(f"val/{k}", v, on_step=False, on_epoch=True, batch_size=int(batch["is_query"].sum()))

    def configure_optimizers(self) -> Any:
        opt = torch.optim.AdamW(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        # Cosine after linear warmup.
        def lr_lambda(step: int) -> float:
            if step < self.warmup_steps:
                return step / max(1, self.warmup_steps)
            progress = (step - self.warmup_steps) / max(1, self.max_steps_ - self.warmup_steps)
            progress = min(max(progress, 0.0), 1.0)
            return 0.5 * (1 + torch.cos(torch.tensor(progress * 3.14159265)).item())
        sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_lambda)
        return {"optimizer": opt, "lr_scheduler": {"scheduler": sched, "interval": "step"}}
