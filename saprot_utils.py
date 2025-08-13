from SaProt.utils.foldseek_util import get_struc_seq
from SaProt.model.saprot.saprot_foldseek_mutation_model import SaprotFoldseekMutationModel
from Bio import SeqIO
import torch

import warnings
warnings.filterwarnings("ignore")

class SaProt_funcs:

    def __init__(self) -> None:
        self.model = self.load_model()

    def load_model(self):

        config = {
            "foldseek_path": None,
            "config_path": "./SaProt/weights/PLMs",
            "load_pretrained": True,
        }
        model = SaprotFoldseekMutationModel(**config)

        device = "cuda" if torch.cuda.is_available() else 'cpu'
        model.eval()
        model.to(device)

        return model

    def align_and_insert_gaps(self, seq1: str, seq2: str) -> str:
        i, j = 0, 0
        aligned_seq2 = []
        
        while i < len(seq1):
            if j < len(seq2) and seq1[i] == seq2[j]:
                aligned_seq2.append(seq2[j])
                j += 1
            else:
                aligned_seq2.append('-')
            i += 1
        
        return ''.join(aligned_seq2)

    def parse_foldseek_ss(self,
                          fasta_file: str,
                          structure_file: str,
                          chain_id: str = 'A',
                          plddt_mask: bool = False,
                          plddt_threshold: float = None) -> str:
        
        seq_fasta = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
        parsed_seqs = get_struc_seq(foldseek="./SaProt/foldseek/foldseek", path=structure_file, chains=chain_id, plddt_mask=plddt_mask, plddt_threshold=plddt_threshold)[chain_id]
        seq_strcut, seq_foldseek, _ = parsed_seqs
        seq_strcut_aligned = self.align_and_insert_gaps(seq_fasta, seq_strcut)
        seq_foldseek_iter = iter(seq_foldseek)
        seq_foldseek = ['#' if aa == '-' else next(seq_foldseek_iter).lower() for aa in seq_strcut_aligned]
        ss = [f"{a}{b}" for a, b in zip(list(seq_fasta), seq_foldseek)]
        
        return ''.join(ss)
        
    def mutation_score(self,
                       sequence_file: str,
                       structure_file: str,
                       chain_id: str,
                       mutations: dict) -> float:
        '''
        mutate_position: Position to be mutated, starting from 1.
        mutation: AA after mutation
        
        '''
        seq = list(self.parse_foldseek_ss(fasta_file=sequence_file, structure_file=structure_file, chain_id=chain_id))

        # No need to modify '#' due to marginalization
        seq = ''.join(seq)
        mut_info = ":".join([f"{seq[2 * (s - 1)]}{s}{m}" for s, m in mutations.items()])

        return round(self.model.predict_mut(seq, mut_info), 3)