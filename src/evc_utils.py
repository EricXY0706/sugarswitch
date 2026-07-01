from Bio import SeqIO
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os

from evcouplings.align.alignment import Alignment
from evcouplings.couplings.mapping import Segment
from evcouplings.align.protocol import describe_frequencies
from evcouplings.couplings.protocol import run

# from src.util import InteractionAnalyzer, SS_TAG, AA_INTERACTIONS

class EVC_funcs:

    def __init__(
        self, 
        alignment_file: str,
        structure_file: str,
        chain_id: str,
        out_dir: str
    ) -> None:
        
        self.aa_list = "-ACDEFGHIKLMNPQRSTVWYX"
        self.chain = chain_id
        self.struct_file = structure_file
        self.a3m_file = alignment_file
        self.a2m_file = str(Path(self.a3m_file).with_suffix(".a2m"))
        self.a3m_aln = list(SeqIO.parse(alignment_file, 'fasta'))
        self.query_seq = str(self.a3m_aln[0].seq)
        os.makedirs(out_dir, exist_ok=True)
        self.out_dir = out_dir
        self.segment = Segment("aa", "101", 1, len(self.query_seq), range(1, len(self.query_seq) + 1)).to_list()
        self.plmc = "./plmc/bin/plmc"

        self.a2m_aln = self.a3m_to_a2m()
        self.freq_file = self.describe_freq()
        self.ec_file = f"{self.out_dir}/_ECs.txt"

    def a3m_to_a2m(
        self, 
        gap_threshold: float = 0.25
    ) -> Alignment:
        '''
        Convert a3m file to a2m file with equal length
        '''
        with open(self.a2m_file, "a") as f:
            for record in self.a3m_aln:
                des = record.description
                seq = "".join([aa for aa in str(record.seq) if aa in self.aa_list])
                if (seq.count("-") / len(seq)) <= gap_threshold:
                    f.write(f">{des}\n{seq}\n")
        
        with open(self.a2m_file, "r") as f:
            aln = Alignment.from_file(f, format="fasta")
        
        return aln
    
    def describe_freq(self) -> Path:
        '''
        Generate conservation statistics file
        '''
        freq = describe_frequencies(alignment=self.a2m_aln, first_index=1, target_seq_index=0)
        freq.to_csv(f"{self.out_dir}/evc_aa_freq.csv", float_format='%.3f', index=False)

        return f"{self.out_dir}/evc_aa_freq.csv"

    def run_evc(
        self,
        protocol: str = "standard",
        scoring_model: str = "skewnormal",
        min_sequence_distance: int = 6,
        theta: float = 0.8,
        focus_mode: bool = True,
        focus_sequence: str = "101",
        segment: Segment = None,
        ignore_gaps: bool = True,
        iterations: int = 100,
        lambda_h: float = 0.01,
        lambda_J: float = 0.01,
        lambda_group: float = None,
        lambda_J_times_Lq: bool = True,
        scale_clusters: float = None,
        cpu: int = 10,
        reuse_ecs: bool = True
    ) -> Path:
        '''
        Run EVcouplings to generate pairwise coupling strength
        '''
        if not os.path.exists(self.ec_file):
            result = run(
                alignment_file=self.a2m_file,
                frequencies_file=self.freq_file,
                prefix=self.out_dir,
                protocol=protocol,                 
                scoring_model=scoring_model,
                min_sequence_distance=min_sequence_distance,
                theta=theta,
                focus_mode=focus_mode,
                focus_sequence=focus_sequence,
                alphabet=self.aa_list,
                segments=segment,
                ignore_gaps=ignore_gaps,
                iterations=iterations,
                lambda_h=lambda_h,
                lambda_J=lambda_J,
                lambda_group=lambda_group,
                lambda_J_times_Lq=lambda_J_times_Lq,
                scale_clusters=scale_clusters,
                cpu=cpu,
                plmc=self.plmc,
                reuse_ecs=reuse_ecs
                )
        else:
            pass
    
    def run_evc_filters(
        self,
        secondary_structure: dict,
        conservation_threshold: float = 0.8,
        evc_threshold: float = 0.5,
        pic_format: str = "pdf"
    ):
        '''
        Run filters
        '''
        lines = open(self.ec_file, 'r').readlines()
        coupling_strength = np.zeros((len(self.query_seq), len(self.query_seq)))
        for line in lines:
            line = line.strip('\n').split(' ')
            res1, res2 = int(line[0]) - 1, int(line[2]) - 1
            score = float(line[-1])
            coupling_strength[res1][res2] = score
            coupling_strength[res2][res1] = score
        plt.figure()
        plt.imshow(coupling_strength, cmap="RdBu_r", vmin=-1, vmax=1)
        plt.title("Evolutionary Coupling Matrix")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(f"{self.out_dir}/evc.{pic_format}", bbox_inches="tight")

        # conservation
        conserve_df = pd.read_csv(self.freq_file)
        conserved_sites = conserve_df.loc[conserve_df.to_numpy()[np.arange(len(conserve_df)), conserve_df.columns.get_indexer(conserve_df["A_i"])] > conservation_threshold, 'i'].tolist()

        # coupling
        coupling = (coupling_strength >= evc_threshold).astype(int)
        coupling_sites = (np.flatnonzero(np.sum(coupling, axis=1) != 0) + 1).tolist()

        return set(conserved_sites) | set(coupling_sites), conserve_df, np.mean(coupling_strength, axis=1)
