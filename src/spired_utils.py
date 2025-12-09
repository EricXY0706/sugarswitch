import numpy as np
from copy import deepcopy
import torch
from Spired.scripts.model import SPIRED_Stab
from Spired.scripts.utils_train_valid import getStabDataTest

class Spired_funcs:
    
    def __init__(self) -> None:
        self.model, self.esm2_650M, self.esm2_3B, self.esm2_batch_converter = self.init_spired_stab()
    
    def init_spired_stab(self):
        '''
        Initialize and configure the SPIRED-Stab model along with pre-trained ESM2 models for protein stability prediction.
        
        Returns:
            tuple: (model, esm2_650M, esm2_3B, esm2_batch_converter) where:
                - model: Initialized SPIRED-Stab stability prediction model
                - esm2_650M: ESM2 model with 650M parameters (for feature extraction)
                - esm2_3B: ESM2 model with 3B parameters (for feature extraction)
                - esm2_batch_converter: Tokenizer/batch processor for ESM2 inputs
        '''
        model = SPIRED_Stab(device_list=["cpu", "cpu", "cpu", "cpu"])
        model.load_state_dict(torch.load(f"./Spired/model/SPIRED-Stab.pth"))
        model.eval()

        esm2_650M, _ = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
        esm2_650M.eval()

        esm2_3B, esm2_alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t36_3B_UR50D")
        esm2_3B.eval()
        esm2_batch_converter = esm2_alphabet.get_batch_converter()
        
        return model, esm2_650M, esm2_3B, esm2_batch_converter
    
    def get_mutation_effect(self,
                            wt_seq: str,
                            mutations: dict):
        '''
        Calculate stability metrics (ddG and dTm) for a mutation by comparing wild-type and mutant sequences.
        
        Args:
            wt_seq (str): Wild-type protein sequence
            mut_seq (str): Mutated protein sequence
            
        Returns:
            tuple: (ddG, dTm) where ddG is the change in Gibbs free energy and dTm is the change in melting temperature
        '''
        wt_seq_list = list(wt_seq)
        mut_seq_list = deepcopy(wt_seq_list)
        for s, m in mutations.items():
            mut_seq_list[s-1] = m
        
        mut_pos_torch_list = torch.tensor((np.array(wt_seq_list) != np.array(mut_seq_list)).astype(int).tolist())
        
        with torch.no_grad():
            
            f1d_esm2_3B, f1d_esm2_650M, target_tokens = getStabDataTest(wt_seq, self.esm2_3B, self.esm2_650M, self.esm2_batch_converter)
            wt_data = {"target_tokens": target_tokens, "esm2-3B": f1d_esm2_3B, "embedding": f1d_esm2_650M}
            
            f1d_esm2_3B, f1d_esm2_650M, target_tokens = getStabDataTest("".join(mut_seq_list), self.esm2_3B, self.esm2_650M, self.esm2_batch_converter)
            mut_data = {"target_tokens": target_tokens, "esm2-3B": f1d_esm2_3B, "embedding": f1d_esm2_650M}
            
            ddG, dTm, _, _ = self.model(wt_data, mut_data, mut_pos_torch_list)
        
        return round(ddG.item(), 3), round(dTm.item(), 3)