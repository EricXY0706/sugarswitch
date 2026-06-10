"""Multi-task Lightning module: binary PPBS interface + continuous ddG regression.

Loss = BCE_logits(p_interface, y_PPBS) [over PPBS-labeled residues]
     + lambda_ddg * Huber(ddg_pred, ddg_true) [over ddg-scanned interface residues]

Selecting `task=binary_only` or `task=ddg_only` zeroes one of the two terms, so a
single module covers all three ablations.
"""
from __future__ import annotations

from typing import Any

import pytorch_lightning as pl
import torch
import torch.nn as nn
import torch.nn.functional as F

from torchmetrics.classification import BinaryAUROC, BinaryAveragePrecision

from ..model.scannetpp import ScanNetPP
from .ddg_regression import _pearson, _spearman


def _aucpr_binary(score: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    """AUCPR when target is already 0/1. Trapezoidal."""
    if target.numel() < 2 or target.sum() == 0 or target.sum() == target.numel():
        return torch.tensor(0.0, device=score.device)
    order = score.argsort(descending=True)
    sorted_b = target[order].float()
    tp = torch.cumsum(sorted_b, dim=0)
    fp = torch.cumsum(1 - sorted_b, dim=0)
    precision = tp / (tp + fp).clamp(min=1)
    recall = tp / sorted_b.sum()
    recall = torch.cat([torch.zeros(1, device=recall.device), recall])
    precision = torch.cat([torch.tensor([1.0], device=precision.device), precision])
    return torch.trapz(precision, recall)


def _auc_roc_binary(score: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    """Mann-Whitney U-statistic form of AUC."""
    if target.numel() < 2 or target.sum() == 0 or target.sum() == target.numel():
        return torch.tensor(0.0, device=score.device)
    pos = score[target.bool()]
    neg = score[~target.bool()]
    # Average rank of positives minus pos-only contribution.
    all_scores = torch.cat([pos, neg])
    ranks = all_scores.argsort().argsort().float() + 1
    pos_rank_sum = ranks[: pos.numel()].sum()
    n_pos = float(pos.numel())
    n_neg = float(neg.numel())
    auc = (pos_rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return auc


class MultiTask(pl.LightningModule):
    def __init__(
        self,
        model_cfg: dict | None = None,
        lr: float = 3e-4,
        weight_decay: float = 1e-2,
        warmup_steps: int = 500,
        max_steps: int = 50_000,
        huber_delta: float = 1.0,
        lambda_ddg: float = 1.0,
        pos_weight: float = 5.0,
        task: str = "multitask",          # one of: multitask, binary_only, ddg_only
        use_sample_weight: bool = True,
        normalize_ddg: bool = False,      # divide Huber by mean(|target|) so loss_ddg is O(1)
        ddg_transform: str = "none",      # "none" or "asinh" (asinh(x/3), tames heavy tail)
        ddg_bins: int = 0,                # if >0, train ddG as K-way classification (T2.d)
        ddg_bin_edges: list | None = None,  # K-1 quantile-based edges; required when ddg_bins>0
        multi_aa_ddg: bool = False,       # if True, ddG head outputs 20 scalars/residue (one per mutant AA)
        contrastive_weight: float = 0.0,  # EGCPPIS-style atom/residue cross-scale aux-loss weight
        contrast_temp: float = 0.1,       # InfoNCE/DCL/SupCon temperature
        contrast_mode: str = "infonce",   # infonce | dcl | vicreg | decur
        contrast_hard_beta: float = 0.0,  # hard-negative reweighting (Robinson) for infonce
        contrast_kc: int = 48,            # DeCUR common-block size (rest of contrast_dim = unique)
        supcon_weight: float = 0.0,       # supervised contrastive over binding labels (on z_res)
    ) -> None:
        super().__init__()
        self.save_hyperparameters()
        cfg = model_cfg or {}
        # When bin-classification is on the model needs to emit K logits per residue
        # instead of a single scalar. Capacity is otherwise identical to E/G.
        if ddg_bins > 0:
            cfg = {**cfg, "ddg_out_dim": ddg_bins}
        if multi_aa_ddg:
            assert ddg_bins == 0, "multi_aa_ddg is mutually exclusive with ddg_bins"
            cfg = {**cfg, "ddg_out_dim": 20}
        self.multi_aa_ddg = bool(multi_aa_ddg)
        self.contrastive_weight = float(contrastive_weight)
        self.contrast_temp = float(contrast_temp)
        self.contrast_mode = str(contrast_mode)
        self.contrast_hard_beta = float(contrast_hard_beta)
        self.contrast_kc = int(contrast_kc)
        self.supcon_weight = float(supcon_weight)
        assert self.contrast_mode in {"infonce", "dcl", "vicreg", "decur"}, self.contrast_mode
        # The model must emit the z_atom/z_res projections whenever any cross-scale
        # or supervised-contrastive term is active.
        if self.contrastive_weight > 0 or self.supcon_weight > 0:
            cfg = {**cfg, "contrastive": True}
            if "contrast_dim" not in cfg:
                cfg["contrast_dim"] = 64
        self.model = ScanNetPP(**cfg)
        self.lr = lr
        self.weight_decay = weight_decay
        self.warmup_steps = warmup_steps
        self.max_steps_ = max_steps
        self.huber_delta = huber_delta
        self.lambda_ddg = lambda_ddg
        self.register_buffer("pos_weight", torch.tensor(float(pos_weight)))
        assert task in {"multitask", "binary_only", "ddg_only"}
        self.task = task
        self.use_sample_weight = use_sample_weight
        self.normalize_ddg = normalize_ddg
        assert ddg_transform in {"none", "asinh"}, ddg_transform
        self.ddg_transform = ddg_transform
        self.ddg_bins = int(ddg_bins)
        if self.ddg_bins > 0:
            assert ddg_bin_edges is not None and len(ddg_bin_edges) == self.ddg_bins - 1, \
                f"need {self.ddg_bins - 1} edges for {self.ddg_bins} bins, got {ddg_bin_edges}"
            edges = torch.tensor(list(ddg_bin_edges), dtype=torch.float32)
            self.register_buffer("ddg_bin_edges", edges)
            # Bin midpoints for Pearson reconstruction. Open bins on the ends —
            # use the boundary +/- half the median interior width as a representative point.
            interior_widths = edges[1:] - edges[:-1]
            w = float(interior_widths.median()) if interior_widths.numel() else 1.0
            mids = [float(edges[0]) - w / 2]
            for i in range(len(edges) - 1):
                mids.append(float((edges[i] + edges[i + 1]) / 2))
            mids.append(float(edges[-1]) + w / 2)
            self.register_buffer("ddg_bin_mids", torch.tensor(mids, dtype=torch.float32))
        # DDP-correct, threshold-free ranking metrics (auto-synced by Lightning).
        self.val_auroc = BinaryAUROC()
        self.val_aucpr = BinaryAveragePrecision()

    def _cross_scale_loss(self, za: torch.Tensor, zr: torch.Tensor) -> torch.Tensor:
        """Cross-scale alignment between atom-tower (za) and residue-tower (zr)
        per-residue projections, [N, d]. Mode selects the objective family."""
        mode, temp = self.contrast_mode, self.contrast_temp
        N = za.shape[0]
        if mode in ("infonce", "dcl"):
            a = F.normalize(za, dim=-1)
            r = F.normalize(zr, dim=-1)
            sim = a @ r.t()                                # cosine [N,N]
            eye = torch.eye(N, dtype=torch.bool, device=za.device)
            if mode == "infonce":
                logits = sim / temp
                if self.contrast_hard_beta > 0:
                    # Robinson hard-negative reweighting: tilt negatives toward
                    # near-anchor (high-sim) pairs. Add a hardness bonus to the
                    # off-diagonal (negative) logits only.
                    logits = logits + self.contrast_hard_beta * sim.masked_fill(eye, 0.0)
                labels = torch.arange(N, device=za.device)
                return 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.t(), labels))
            # DCL: delete the positive term from the denominator (decoupled).
            s = sim / temp
            pos = s.diagonal()
            neg = s.masked_fill(eye, float("-inf"))
            row = (-pos + torch.logsumexp(neg, dim=1)).mean()
            col = (-pos + torch.logsumexp(neg.t(), dim=1)).mean()
            return 0.5 * (row + col)
        if mode == "vicreg":
            sim_loss = F.mse_loss(za, zr)
            def _std(x):
                std = torch.sqrt(x.var(dim=0) + 1e-4)
                return torch.mean(F.relu(1.0 - std))
            std_loss = _std(za) + _std(zr)
            def _cov(x):
                xc = x - x.mean(dim=0)
                cov = (xc.t() @ xc) / max(1, N - 1)
                off = cov - torch.diag(torch.diagonal(cov))
                return (off ** 2).sum() / x.shape[1]
            cov_loss = _cov(za) + _cov(zr)
            return 25.0 * sim_loss + 25.0 * std_loss + 1.0 * cov_loss
        # decur: Barlow-style cross-correlation; common block -> identity,
        # unique block (modality-specific dims) -> decorrelated to zero.
        za_n = (za - za.mean(0)) / (za.std(0) + 1e-4)
        zr_n = (zr - zr.mean(0)) / (zr.std(0) + 1e-4)
        C = (za_n.t() @ zr_n) / N                          # cross-correlation [d,d]
        kc = min(self.contrast_kc, C.shape[0])
        Cc = C[:kc, :kc]
        on = ((torch.diagonal(Cc) - 1.0) ** 2).sum()
        off = (Cc ** 2).sum() - (torch.diagonal(Cc) ** 2).sum()
        L_com = on + 0.005 * off
        Cu = C[kc:, kc:]
        L_uni = (Cu ** 2).sum()
        return L_com + L_uni

    def _supcon_loss(self, z: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        """Supervised contrastive (SupCon) over binding labels: pull same-class
        (interface vs non-interface) residues together. z [M,d], y [M] in {0,1}."""
        z = F.normalize(z, dim=-1)
        M = z.shape[0]
        sim = (z @ z.t()) / self.contrast_temp
        eye = torch.eye(M, dtype=torch.bool, device=z.device)
        sim = sim.masked_fill(eye, float("-inf"))
        logp = sim - torch.logsumexp(sim, dim=1, keepdim=True)
        same = (y.unsqueeze(0) == y.unsqueeze(1)) & ~eye
        denom = same.sum(1).clamp(min=1)
        # masked_fill (not multiply) so the -inf on the diagonal never hits a 0 weight
        # (-inf * 0 = NaN); non-positive entries contribute exactly 0.
        per = -logp.masked_fill(~same, 0.0).sum(1) / denom
        # only anchors that have at least one same-class positive contribute
        valid = same.sum(1) > 0
        return per[valid].mean() if valid.any() else z.sum() * 0.0

    def forward(self, batch: dict) -> dict:
        return self.model(
            coords=batch["coords"],
            res_types=batch["res_types"],
            chain_ids=batch["chain_ids"],
            is_query=batch["is_query"],
            edge_index=batch["edge_index"],
            seq_idx=batch["seq_idx"],
            return_dict=True,
            esm_features=batch.get("esm_features"),
            pssm_features=batch.get("pssm_features"),
            sidechain_features=batch.get("sidechain_features"),
            atom=({"xyz": batch["atom_xyz"], "elem": batch["atom_elem"],
                   "scflag": batch["atom_scflag"], "edge_index": batch["atom_edge_index"],
                   "resid": batch["atom_resid"]} if "atom_xyz" in batch else None),
        )

    def _loss(self, out: dict, batch: dict) -> tuple[torch.Tensor, dict]:
        ddg_pred = out["ddg"]
        bin_logit = out["binary_logit"]
        is_query = batch["is_query"]
        is_interface = batch["is_interface"]
        ddg_mask = batch.get("ddg_mask", is_interface)
        target_ddg = batch["ddg"]
        target_bin = batch["is_interface_ppbs"].float()
        ppbs_mask = batch["ppbs_mask"]
        ppbs_w = batch["ppbs_weight"] if self.use_sample_weight else torch.ones_like(batch["ppbs_weight"])

        device = ddg_pred.device
        zero = ddg_pred.sum() * 0.0

        # --- Binary BCE ---
        bce_mask = ppbs_mask & is_query
        if bce_mask.any() and self.task != "ddg_only":
            target_b = target_bin[bce_mask]
            logit_b = bin_logit[bce_mask]
            w = ppbs_w[bce_mask]
            pw = torch.where(target_b > 0.5, self.pos_weight, torch.ones_like(self.pos_weight))
            bce_per = F.binary_cross_entropy_with_logits(
                logit_b, target_b, reduction="none"
            ) * pw * w
            loss_bce = bce_per.sum() / (w * pw).sum().clamp(min=1e-6)
        else:
            loss_bce = zero

        # --- DDG loss: Huber regression OR K-way ordinal cross-entropy ---
        scored = ddg_mask & is_interface & is_query
        if scored.any() and self.task != "binary_only":
            pred_s = ddg_pred[scored]              # [N_s] or [N_s, K] or [N_s, 20]
            targ_s = target_ddg[scored]            # [N_s], raw kcal/mol
            if self.ddg_bins > 0:
                # Ordinal-bin head: pred is [N_s, K] logits; target binned by quantile edges.
                target_idx = torch.bucketize(targ_s, self.ddg_bin_edges)  # 0..K-1
                loss_ddg = F.cross_entropy(pred_s, target_idx, reduction="mean")
            else:
                if self.multi_aa_ddg:
                    # Per-residue 20-way scalar head; gather the column corresponding
                    # to this row's mutant AA so the loss is a regression against the
                    # single experimental ΔΔG value at that (residue, mutant_AA) pair.
                    target_aa = batch["ddg_target_aa"][scored]            # [N_s]
                    pred_s = pred_s.gather(1, target_aa.unsqueeze(-1)).squeeze(-1)  # [N_s]
                if self.ddg_transform == "asinh":
                    # asinh(x): near-linear for |x|<1, sublinear thereafter
                    # (asinh(15)≈3.4, asinh(86)≈5.2). Keeps gradient strength near
                    # the median while damping the right tail.
                    pred_s = torch.asinh(pred_s)
                    targ_s = torch.asinh(targ_s)
                loss_ddg = F.huber_loss(pred_s, targ_s, delta=self.huber_delta, reduction="mean")
                if self.normalize_ddg:
                    scale = targ_s.abs().mean().detach().clamp(min=1.0)
                    loss_ddg = loss_ddg / scale
        else:
            loss_ddg = zero

        if self.task == "binary_only":
            loss = loss_bce
        elif self.task == "ddg_only":
            loss = loss_ddg
        else:
            loss = loss_bce + self.lambda_ddg * loss_ddg

        # --- Cross-scale atom/residue alignment aux loss (mode-selectable) ---
        loss_contrast = zero
        if self.contrastive_weight > 0 and "z_atom" in out:
            cmask = is_query
            za = out["z_atom"][cmask]
            zr = out["z_res"][cmask]
            if za.shape[0] >= 2:
                loss_contrast = self._cross_scale_loss(za, zr)
                loss = loss + self.contrastive_weight * loss_contrast

        # --- Supervised contrastive over binding labels (optional) ---
        loss_supcon = zero
        if self.supcon_weight > 0 and "z_res" in out:
            smask = ppbs_mask & is_query
            if smask.sum() >= 4:
                loss_supcon = self._supcon_loss(out["z_res"][smask], target_bin[smask])
                loss = loss + self.supcon_weight * loss_supcon

        with torch.no_grad():
            if bce_mask.any():
                p = torch.sigmoid(bin_logit[bce_mask])
                t = target_bin[bce_mask]
                aucpr = _aucpr_binary(p, t)
                auc = _auc_roc_binary(p, t)
            else:
                aucpr = torch.tensor(0.0, device=device)
                auc = torch.tensor(0.0, device=device)
            if scored.any():
                if self.ddg_bins > 0:
                    # Continuous prediction via softmax-weighted expectation over
                    # bin midpoints — rank-preserving, uses the full distribution
                    # information instead of collapsing to argmax.
                    probs = F.softmax(ddg_pred[scored], dim=-1)
                    pred_real = (probs * self.ddg_bin_mids).sum(dim=-1)
                    pearson = _pearson(pred_real, target_ddg[scored])
                    spearman = _spearman(pred_real, target_ddg[scored])
                elif self.multi_aa_ddg:
                    target_aa = batch["ddg_target_aa"][scored]
                    pred_real = ddg_pred[scored].gather(1, target_aa.unsqueeze(-1)).squeeze(-1)
                    pearson = _pearson(pred_real, target_ddg[scored])
                    spearman = _spearman(pred_real, target_ddg[scored])
                else:
                    pearson = _pearson(ddg_pred[scored], target_ddg[scored])
                    spearman = _spearman(ddg_pred[scored], target_ddg[scored])
            else:
                pearson = torch.tensor(0.0, device=device)
                spearman = torch.tensor(0.0, device=device)

        # Skip head-specific keys when their mask is empty so epoch-mean logs aren't
        # biased by zeros from batches that don't supervise that head (matters for
        # mixed PPBS+SKEMPI val where some batches have no PPBS residues and others
        # have no ddG residues).
        out_metrics: dict = {"loss": loss.detach()}
        if self.contrastive_weight > 0 and "z_atom" in out:
            out_metrics["loss_contrast"] = loss_contrast.detach()
        if self.supcon_weight > 0 and "z_res" in out:
            out_metrics["loss_supcon"] = loss_supcon.detach()
        if bce_mask.any():
            out_metrics["loss_bce"] = loss_bce.detach()
            out_metrics["aucpr"] = aucpr
            out_metrics["auc"] = auc
        if scored.any():
            out_metrics["loss_ddg"] = loss_ddg.detach()
            out_metrics["pearson"] = pearson
            out_metrics["spearman"] = spearman
        return loss, out_metrics

    def _batch_size(self, batch: dict) -> int:
        return int(max(1, batch["is_query"].sum().item()))

    def training_step(self, batch: dict, batch_idx: int) -> torch.Tensor:
        out = self(batch)
        loss, metrics = self._loss(out, batch)
        bs = self._batch_size(batch)
        for k, v in metrics.items():
            self.log(f"train/{k}", v, on_step=True, on_epoch=True, batch_size=bs)
        return loss

    def validation_step(self, batch: dict, batch_idx: int) -> None:
        out = self(batch)
        _, metrics = self._loss(out, batch)
        bs = self._batch_size(batch)
        # Per-batch scalars (loss, pearson, etc.) go through standard logging.
        for k, v in metrics.items():
            if k in {"aucpr", "auc"}:
                continue  # superseded by torchmetrics accumulators below
            self.log(f"val/{k}", v, on_step=False, on_epoch=True, batch_size=bs)
        # Rank-correct AUROC / AP across the whole val set — but only when the
        # task uses the binary head at all. In ddg_only fine-tune (e.g. SKEMPI),
        # ppbs_mask is all-False on every batch; the torchmetrics accumulators
        # would compute on an empty state at epoch end and crash. Skip cleanly.
        if self.task == "ddg_only":
            return
        bin_logit = out["binary_logit"]
        ppbs_mask = batch["ppbs_mask"] & batch["is_query"]
        if ppbs_mask.any():
            p = torch.sigmoid(bin_logit[ppbs_mask]).detach()
            t = batch["is_interface_ppbs"][ppbs_mask].long()
            self.val_auroc.update(p, t)
            self.val_aucpr.update(p, t)
        self.log("val/auc", self.val_auroc, on_epoch=True, batch_size=bs)
        self.log("val/aucpr", self.val_aucpr, on_epoch=True, batch_size=bs)

    def configure_optimizers(self) -> Any:
        opt = torch.optim.AdamW(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        def lr_lambda(step: int) -> float:
            if step < self.warmup_steps:
                return step / max(1, self.warmup_steps)
            progress = (step - self.warmup_steps) / max(1, self.max_steps_ - self.warmup_steps)
            progress = min(max(progress, 0.0), 1.0)
            return 0.5 * (1 + torch.cos(torch.tensor(progress * 3.14159265)).item())
        sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_lambda)
        return {"optimizer": opt, "lr_scheduler": {"scheduler": sched, "interval": "step"}}
