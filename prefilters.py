from utils import FastaLoader, StructureLoader, InteractionCheck, ClashCheck, GlycanMover, BordaCount, SS_TAG
from evc_utils import EVC_funcs
from rosetta_utils import Rosetta_funcs
from saprot_utils import SaProt_funcs
from spired_utils import Spired_funcs

from config import basic_configs, ranker_configs
from pathlib import Path
from tqdm import *
import pandas as pd
import numpy as np
import os
import warnings

def update_infer(
        input_fasta_file: str,
) -> None:
    """
    Update the infer result dir with the input fasta file.
    """
    if not os.path.exists(input_fasta_file):
        raise FileNotFoundError(f"Input fasta file `{input_fasta_file}` not found.")
    query_seq = FastaLoader.get_sequence(sequence_file=input_fasta_file)


def run_prefilters(
        input_fasta_file: str,
        input_structure_file: str,
        input_alignment_file: str,
        output_dir: str,
) -> None:
    """
    Run prefilters on the input fasta file.
    """
    warnings.filterwarnings("ignore")
    os.makedirs(output_dir, exist_ok=True)
    chain = basic_configs["protein_chain_id"]
    query_sequence = FastaLoader.get_sequence(sequence_file=input_fasta_file)
    filename = Path(input_structure_file).name.split('.')[0]
    ss = StructureLoader.get_secondary_structure(structure_file=input_structure_file)

    # Filtering out interacting sites with the given hotspots
    interaction_checker = InteractionCheck()
    hotspots_interacting_sites = interaction_checker.get_interaction_aa(
        structure_file=input_structure_file,
        positions=basic_configs["functional_hotspots"],
        dist_threshold=basic_configs["Cb_interaction_threshold"],
        is_self_included=True,
        num_neighbors=basic_configs["num_neighbors_to_shield"],
    )
    
    # Filtering out the strong-coupling and conserved sites
    evc = EVC_funcs(alignment_file=input_alignment_file, structure_file=input_structure_file, chain_id=chain, out_dir=f"{output_dir}/")
    evc.run_evc(
        min_sequence_distance=basic_configs["evc_min_sequence_distance"],
        theta=basic_configs["evc_theta"],
        iterations=basic_configs["evc_num_iterations"],
        lambda_h=basic_configs["evc_lambda_h"],
        lambda_J=basic_configs["evc_lambda_J"],
        cpu=basic_configs["evc_num_cpu"],
    )
    conserverd_coupling_sites, conservation_df, coupling_stength = evc.run_evc_filters(
        secondary_structure=ss,
        conservation_thresholds=basic_configs["conservation_threshold"],
        evc_threshold=basic_configs["evc_coupling_threshold"],
    )
    
    # Filtering out the low SASA sites
    rosetta = Rosetta_funcs()
    _, low_sasa_sites, sasa_index_dict = rosetta.get_SASA(
        structure_file=input_structure_file,
        cutoff=basic_configs["sasa_cutoff"],
    )

    non_editable_regions = hotspots_interacting_sites | conserverd_coupling_sites | low_sasa_sites
    editable_regions = set(list(range(1, len(query_sequence)+1))) - non_editable_regions

    # Modification pipeline
    results = []
    saprot = SaProt_funcs()
    spired = Spired_funcs()
    for s in tqdm(editable_regions, dynamic_ncols=True):

        conservation_score = conservation_df.loc[conservation_df["i"] == s, "conservation"].values[0]
        coupling_score = round(coupling_stength[s-1], 3)
        sasa_value = round(sasa_index_dict[s], 3)
        sasa_value_before1 = round(sasa_index_dict[s-1], 3) if s != 1 else 0.
        sasa_value_next1 = round(sasa_index_dict[s+1], 3) if s != len(query_sequence) else 0.
        sasa_value_next2 = round(sasa_index_dict[s+2], 3) if s != len(query_sequence) and s != (len(query_sequence) - 1) else 0.
        sasa_around_mean = round((sasa_value_before1 + sasa_value + sasa_value_next1) / 3, 3)
        sasa_next_mean = round((sasa_value + sasa_value_next1 + sasa_value_next2) / 3, 3)
        glycoprotein_structure_file = f"{output_dir}/{filename}_{list(query_sequence)[s-1]}{s}N.pdb"

        mut_score_s = saprot.mutation_score(
            sequence_file=input_fasta_file,
            structure_file=input_structure_file,
            chain_id=chain,
            mutations={s: "N"},
        )
        mut_score_s_next2_S = saprot.mutation_score(
            sequence_file=input_fasta_file,
            structure_file=input_structure_file,
            chain_id=chain,
            mutations={s: "N", s+2: "S"},
        ) if s != len(query_sequence) and s != (len(query_sequence) - 1) else None
        mut_score_s_next2_T = saprot.mutation_score(
            sequence_file=input_fasta_file,
            structure_file=input_structure_file,
            chain_id=chain,
            mutations={s: "N", s+2: "T"},
        ) if s != len(query_sequence) and s != (len(query_sequence) - 1) else None

        ddG_s, dTm_s = spired.get_mutation_effect(
            wt_seq=query_sequence,
            mutations={s: "N"},
        )
        ddG_next2_S, dTm_next2_S = spired.get_mutation_effect(
            wt_seq=query_sequence,
            mutations={s: "N", s+2: "S"},
        ) if s != len(query_sequence) and s != (len(query_sequence) - 1) else None
        ddG_next2_T, dTm_next2_T = spired.get_mutation_effect(
            wt_seq=query_sequence,
            mutations={s: "N", s+2: "S"},
        ) if s != len(query_sequence) and s != (len(query_sequence) - 1) else None
        
        rosetta.mutate(
            structure_file=input_structure_file,
            output_file=glycoprotein_structure_file,
            mutate_position=s,
            mutation="N",
        )
        ss_dict = StructureLoader.get_secondary_structure(structure_file=glycoprotein_structure_file)
        ss_site = SS_TAG[ss_dict[(chain, s)]][-1]
        
        glycan_mover = GlycanMover(
            bond_length=basic_configs["bond_length_C1_ND2"],
            angle_C1=basic_configs["angle_C1_ND2_CG"],
            angle_C2=basic_configs["angle_C2_C1_ND2"],
            angle_O5=basic_configs["angle_O5_C1_ND2"],
            dihedral_C1=basic_configs["dihedral_C1_ND2_CG_CB"].get(ss_site, -120.0),
            dihedral_C2=basic_configs["dihedral_C2_C1_ND2_CG"].get(ss_site, 120.0),
            dihedral_O5=basic_configs["dihedral_O5_C1_ND2_CG"].get(ss_site, -120.0),
        )
        glycan_mover.move(
            protein_structure_file=glycoprotein_structure_file,
            glycan_structure_file="./G51766DQ.pdb",
            output_pdb=glycoprotein_structure_file,
            asn_res_id=s,
            protein_chain=chain,
        )
        clash_checker = ClashCheck()
        clash_residues = clash_checker.has_clash(
            chain_id=chain,
            secondary_structure=ss_dict,
            structure_file=glycoprotein_structure_file,
        )

        results.append([s, SS_TAG[ss_dict[(chain, s)]][0], f"{list(query_sequence)[s-1]}{s}N", conservation_score, coupling_score, 
                        sasa_value, sasa_value_next1, sasa_value_next2, sasa_around_mean, sasa_next_mean, 
                        ddG_s, dTm_s, ddG_next2_S, dTm_next2_S, ddG_next2_T, dTm_next2_T, 
                        mut_score_s, mut_score_s_next2_S, mut_score_s_next2_T, clash_residues])

    df = pd.DataFrame(results)
    df.columns = ["Site", "SS", "Mutation", "Conservation_score", "Coupling_score", 
                  "SASA", "SASA_next1", "SASA_next2", "SASA_around_mean", "SASA_next_mean", 
                  "ddG", "dTm", "ddG_S", "dTm_S", "ddG_T", "dTm_T", 
                  "Mut_score", "Mut_score_S", "Mut_score_T", "Clash"]
    ranker = BordaCount(**ranker_configs)
    df = ranker.compute_score(df)
    df.to_csv(f"{output_dir}/uricase.xy0804.csv", index=False)