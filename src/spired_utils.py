from copy import deepcopy
import numpy as np
import torch
from Spired.scripts.model import SPIRED_Stab

class Spired_funcs:

    def __init__(self) -> None:
        # Seven-GPU layout in this process:
        # cuda:0     ESM2-3B
        # cuda:1     ESM2-650M
        # cuda:2-5   SPIRED folding units 1-4
        # cuda:6     stability head
        self.esm3b_device = torch.device("cuda:0")
        self.esm650m_device = torch.device("cuda:1")
        self.spired_devices = ["cuda:2", "cuda:3", "cuda:4", "cuda:5"]
        self.spired_input_device = torch.device(self.spired_devices[0])
        self.stab_device = torch.device("cuda:6")

        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is not available")
        if torch.cuda.device_count() < 7:
            raise RuntimeError(
                "The seven-GPU SPIRED layout requires at least 7 visible GPUs, "
                f"but only {torch.cuda.device_count()} are visible"
            )

        (
            self.model,
            self.esm2_650M,
            self.esm2_3B,
            self.esm2_batch_converter,
        ) = self.init_spired_stab()

        self._cached_wt_sequence = None
        self._cached_wt_features = None
        self._validate_devices()

    def init_spired_stab(self):
        
        model = SPIRED_Stab(
            device_list=self.spired_devices,
            stab_device=self.stab_device,
        )
        state_dict = torch.load(
            "./Spired/model/SPIRED-Stab.pth",
            map_location="cpu",
        )
        model.load_state_dict(state_dict)
        model.eval()

        esm2_650M, _ = torch.hub.load(
            "facebookresearch/esm:main",
            "esm2_t33_650M_UR50D",
        )
        esm2_650M = esm2_650M.to(self.esm650m_device).eval()

        esm2_3B, esm2_alphabet = torch.hub.load(
            "facebookresearch/esm:main",
            "esm2_t36_3B_UR50D",
        )
        esm2_3B = esm2_3B.to(self.esm3b_device).eval()
        esm2_batch_converter = esm2_alphabet.get_batch_converter()

        return model, esm2_650M, esm2_3B, esm2_batch_converter

    def _validate_devices(self):
        
        expected_devices = {
            "ESM2-3B": self.esm3b_device,
            "ESM2-650M": self.esm650m_device,
            "SPIRED block0": torch.device("cuda:2"),
            "SPIRED block1": torch.device("cuda:3"),
            "SPIRED block2": torch.device("cuda:4"),
            "SPIRED block3": torch.device("cuda:5"),
            "Stab head": self.stab_device,
        }
        actual_devices = {
            "ESM2-3B": next(self.esm2_3B.parameters()).device,
            "ESM2-650M": next(self.esm2_650M.parameters()).device,
            "SPIRED block0": next(self.model.SPIRED.block0.parameters()).device,
            "SPIRED block1": next(self.model.SPIRED.block1.parameters()).device,
            "SPIRED block2": next(self.model.SPIRED.block2.parameters()).device,
            "SPIRED block3": next(self.model.SPIRED.block3.parameters()).device,
            "Stab head": next(self.model.Stab.parameters()).device,
        }

        for name, expected_device in expected_devices.items():
            actual_device = actual_devices[name]
            if actual_device != expected_device:
                raise RuntimeError(
                    f"{name} is on {actual_device}, expected {expected_device}"
                )

    def _extract_esm_features(self, sequence: str):
        
        batch = [("protein", sequence)]
        _, _, batch_tokens = self.esm2_batch_converter(batch)

        tokens_3b = batch_tokens.to(self.esm3b_device)
        output_3b = self.esm2_3B(
            tokens_3b,
            repr_layers=list(range(37)),
            return_contacts=False,
        )
        # [B, tokens, 37, 2560] -> remove BOS and EOS -> SPIRED unit 1.
        f1d_esm2_3B = torch.stack(
            [output_3b["representations"][layer] for layer in range(37)],
            dim=2,
        )[:, 1:-1].detach().to(self.spired_input_device)
        del output_3b, tokens_3b

        tokens_650m = batch_tokens.to(self.esm650m_device)
        output_650m = self.esm2_650M(
            tokens_650m,
            repr_layers=[33],
            return_contacts=False,
        )
        # [B, tokens, 1280] -> remove BOS and EOS -> stability head.
        f1d_esm2_650M = (
            output_650m["representations"][33][:, 1:-1]
            .detach()
            .to(self.stab_device)
        )
        del output_650m, tokens_650m

        target_tokens = batch_tokens[:, 1:-1].detach().to(
            self.spired_input_device
        )
        return f1d_esm2_3B, f1d_esm2_650M, target_tokens

    def _get_wt_features(self, wt_seq: str):
        
        if (
            self._cached_wt_sequence != wt_seq
            or self._cached_wt_features is None
        ):
            self._cached_wt_features = self._extract_esm_features(wt_seq)
            self._cached_wt_sequence = wt_seq
        return self._cached_wt_features

    def prepare_wt(self, wt_seq: str):
        
        with torch.inference_mode():
            self._get_wt_features(wt_seq)

    def clear_wt_cache(self):
        
        self._cached_wt_sequence = None
        self._cached_wt_features = None
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    def get_mutation_effect(self, wt_seq: str, mutations: dict):
        
        wt_seq_list = list(wt_seq)
        mut_seq_list = deepcopy(wt_seq_list)
        for position, amino_acid in mutations.items():
            mut_seq_list[position - 1] = amino_acid
        mut_seq = "".join(mut_seq_list)

        mutation_mask = (
            np.asarray(wt_seq_list) != np.asarray(mut_seq_list)
        ).astype(np.float32)
        mut_pos_torch_list = torch.as_tensor(
            mutation_mask,
            dtype=torch.float32,
            device=self.stab_device,
        )

        with torch.inference_mode():
            (
                wt_f1d_esm2_3B,
                wt_f1d_esm2_650M,
                wt_target_tokens,
            ) = self._get_wt_features(wt_seq)
            (
                mut_f1d_esm2_3B,
                mut_f1d_esm2_650M,
                mut_target_tokens,
            ) = self._extract_esm_features(mut_seq)

            wt_data = {
                "target_tokens": wt_target_tokens,
                "esm2-3B": wt_f1d_esm2_3B,
                "embedding": wt_f1d_esm2_650M,
            }
            mut_data = {
                "target_tokens": mut_target_tokens,
                "esm2-3B": mut_f1d_esm2_3B,
                "embedding": mut_f1d_esm2_650M,
            }
            ddG, dTm, _, _ = self.model(
                wt_data,
                mut_data,
                mut_pos_torch_list,
            )

        return round(ddG.item(), 3), round(dTm.item(), 3)