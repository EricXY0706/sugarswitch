from src.evc_utils import EVC_funcs
from src.rosetta_utils import Rosetta_funcs
from src.saprot_utils import SaProt_funcs
from src.spired_utils import Spired_funcs
from src.gvpbind_utils import gvpbind_predict

from src.util import *
from config import basic_configs, ranker_configs

from pathlib import Path
from tqdm import *
import ast
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
import warnings

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))

def run_prefilters(
    input_fasta_file: str,
    input_structure_file: str,
    output_dir: str,
    name: str = None,
    suffix: str = "",
    protein_chain_id: str = "A",
    functional_hotspots: list = [],
) -> None:
    """
    Run prefilters on the input fasta file.
    """
    warnings.filterwarnings("ignore")
    filename = name if name else Path(input_fasta_file).name.split(".")[0]
    structure_file = update_infer(input_fasta_file=input_fasta_file, input_structure_file=input_structure_file, output_dir=output_dir, filename=filename, suffix=suffix)
    query_sequence = FastaLoader.get_sequence(sequence_file=input_fasta_file, chain_id=protein_chain_id)
    chains_nums = [int(rec.description.split(" ")[-1].split(":")[-1]) for rec in SeqIO.parse(input_fasta_file, "fasta")]
    chains_nums.append(len(CHAIN_IDS) - sum(chains_nums))
    seq_chain_ids = []
    pos = 0
    aln_file_id = 1
    for i, l in enumerate(chains_nums):
        seq_chain_id = "".join(CHAIN_IDS[pos: pos + l])
        seq_chain_ids.append(seq_chain_id)
        pos += l
        if protein_chain_id in seq_chain_id:
            aln_file_id = i + 1

    ss = StructureLoader.get_secondary_structure(structure_file=structure_file, chain_id=protein_chain_id)
    
    # Scanning all the inter- and intra- chain interactions (denote them as functional sites)
    analyzer = InteractionAnalyzer(structure_file=structure_file, chain_id=protein_chain_id)
    interaction_dict = analyzer.analyze()
    ppi_sites = analyzer.extract_ppi_sites(chain_id=protein_chain_id, result=interaction_dict)
    hotspots_sites = analyzer.extract_sites_interacting_with_hotspots(chain_id=protein_chain_id, hotspots=ast.literal_eval(functional_hotspots), result=interaction_dict)
    
    # Filtering out the strong-coupling and conserved sites
    skip_evc = False
    input_alignment_file = f"{output_dir}/msa/uniref_{aln_file_id}.a3m"
    if len(list(SeqIO.parse(input_alignment_file, "fasta"))) >= 10:
        evc = EVC_funcs(alignment_file=input_alignment_file, structure_file=structure_file, chain_id=protein_chain_id, out_dir=f"{output_dir}/evc/")
        evc.run_evc(
            focus_sequence=filename,
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
    else:
        skip_evc = True
        conserverd_coupling_sites = set()

    # Filtering out the low SASA sites
    rosetta = Rosetta_funcs()
    sasa_cutoff, low_sasa_sites, sasa_index_dict = rosetta.get_SASA(
        structure_file=structure_file,
        cutoff=basic_configs["sasa_cutoff"],
        chain=protein_chain_id,
    )
    
    # Filtering out the highly-potencial PPI interface sites
    gvpbind_score_by_res = gvpbind_predict(
        structure_file=structure_file,
        chain_id=protein_chain_id,
        out_dir=output_dir,
        filename=filename,
    )
    gvpbind_sites = set([site for site, score in gvpbind_score_by_res.items() if score > basic_configs["gvpbind_cutoff"] and site <= len(query_sequence)])

    non_editable_regions = ppi_sites | hotspots_sites | conserverd_coupling_sites | low_sasa_sites | gvpbind_sites
    editable_regions = set(list(range(1, len(query_sequence)+1))) - non_editable_regions
    editable_regions = set([s for s in editable_regions if not (list(query_sequence)[s-1] == "P" and SS_TAG[ss[(protein_chain_id, s)]][0] == "turn")])
    
    # Modification pipeline
    results = []
    saprot = SaProt_funcs()
    spired = Spired_funcs()
    
    for s in tqdm(editable_regions, dynamic_ncols=True):
        
        if not skip_evc:
            conservation_score = conservation_df.loc[conservation_df["i"] == s, "conservation"].values[0]
            coupling_score = round(coupling_stength[s-1], 3)
        else:
            conservation_score, coupling_score = 0., 0.
        sasa_value = round(sasa_index_dict[s], 3)
        sasa_value_before1 = round(sasa_index_dict[s-1], 3) if s != 1 else sasa_value
        sasa_value_next1 = round(sasa_index_dict[s+1], 3) if s != len(query_sequence) else sasa_value
        sasa_value_next2 = round(sasa_index_dict[s+2], 3) if s != len(query_sequence) and s != (len(query_sequence) - 1) else sasa_value
        sasa_around_mean = round((sasa_value_before1 + sasa_value + sasa_value_next1) / 3, 3)
        sasa_next_mean = round((sasa_value + sasa_value_next1 + sasa_value_next2) / 3, 3)
        gvp_unbind_score = round(1 - gvpbind_score_by_res[s], 3)
        os.makedirs(f"{output_dir}/glycans/", exist_ok=True)
        glycoprotein_structure_file = f"{output_dir}/glycans/{filename}_{list(query_sequence)[s-1]}{s}N.pdb"

        mut_score_s = saprot.mutation_score(
            query_seq=query_sequence,
            structure_file=structure_file,
            chain_id=protein_chain_id,
            mutations={s: "N"},
        )
        if s != len(query_sequence) and s != (len(query_sequence) - 1):
            mut_score_s_next2_S = saprot.mutation_score(
                query_seq=query_sequence,
                structure_file=structure_file,
                chain_id=protein_chain_id,
                mutations={s: "N", s+2: "S"},
            )
            mut_score_s_next2_T = saprot.mutation_score(
                query_seq=query_sequence,
                structure_file=structure_file,
                chain_id=protein_chain_id,
                mutations={s: "N", s+2: "T"},
            )
        else:
            mut_score_s_next2_S, mut_score_s_next2_T = mut_score_s, mut_score_s

        ddG_s, dTm_s = spired.get_mutation_effect(
            wt_seq=query_sequence,
            mutations={s: "N"},
        )
        if s != len(query_sequence) and s != (len(query_sequence) - 1):
            ddG_next2_S, dTm_next2_S = spired.get_mutation_effect(
                wt_seq=query_sequence,
                mutations={s: "N", s+2: "S"},
            )
            ddG_next2_T, dTm_next2_T = spired.get_mutation_effect(
                wt_seq=query_sequence,
                mutations={s: "N", s+2: "T"},
            )
        else:
            ddG_next2_S, dTm_next2_S, ddG_next2_T, dTm_next2_T = ddG_s, dTm_s, ddG_s, dTm_s
        
        rosetta.mutate(
            structure_file=structure_file,
            output_file=glycoprotein_structure_file,
            mutations={f"{protein_chain_id}{s}": "N"},
        )
        ss_dict = StructureLoader.get_secondary_structure(structure_file=glycoprotein_structure_file, chain_id=protein_chain_id)
        ss_site = SS_TAG[ss_dict[(protein_chain_id, s)]][-1]
        
        glycanmover = GlycanMover(
            bond_length=basic_configs["bond_length_C1_ND2"],
            angle_C1=basic_configs["angle_C1_ND2_CG"],
            angle_C2=basic_configs["angle_C2_C1_ND2"],
            angle_O5=basic_configs["angle_O5_C1_ND2"],
            dihedral_C1=basic_configs["dihedral_C1_ND2_CG_CB"].get(ss_site, 178.5),
            dihedral_C2=basic_configs["dihedral_C2_C1_ND2_CG"].get(ss_site, 90.0),
            dihedral_O5=basic_configs["dihedral_O5_C1_ND2_CG"].get(ss_site, -95.0),
        )
        glycanmover.move(
            protein_structure_file=glycoprotein_structure_file,
            glycan_structure_file="./src/G51766DQ.pdb",
            output_pdb=glycoprotein_structure_file,
            glycan_positions={protein_chain_id: [s]},
        )
        clash_checker = ClashCheck()
        clash_residues = clash_checker.has_clash(
            chain_id=protein_chain_id,
            secondary_structure=ss_dict,
            structure_file=glycoprotein_structure_file,
        )

        results.append([s, SS_TAG[ss_dict[(protein_chain_id, s)]][0], f"{list(query_sequence)[s-1]}{s}N", conservation_score, coupling_score, 
                        sasa_value, sasa_value_next1, sasa_value_next2, sasa_around_mean, sasa_next_mean, gvp_unbind_score,
                        ddG_s, dTm_s, ddG_next2_S, dTm_next2_S, ddG_next2_T, dTm_next2_T, 
                        mut_score_s, mut_score_s_next2_S, mut_score_s_next2_T, clash_residues])

    df = pd.DataFrame(results)
    df.columns = ["Site", "SS", "Mutation", "ConservationScore", "CouplingScore",
                  "SASA_i", "SASA_i+1", "SASA_i+2", "SASA_(i-1:i+1)", "SASA_(i:i+2)", "GVP_unbind_score",
                  "ddG", "dTm", "ddG_NXS", "dTm_NXS", "ddG_NXT", "dTm_NXT", 
                  "MutScore", "MutScore_NXS", "MutScore_NXT", "Clash"]
    ranker = BordaCount(**ranker_configs)
    df = ranker.compute_score(df)
    df_file = f"{output_dir}/{filename}_single_points.csv"
    df.to_csv(df_file, index=False)
    pose = rosetta.get_pose(structure_file)
    StructureFileEditor.write_score_as_bfactor(
        pose=pose,
        structure_file=structure_file,
        chain_id=protein_chain_id,
        df_file=df_file,
    )
    reporter = prefilter_report(
        input_fasta_file=input_fasta_file,
        query_sequence=query_sequence,
        editable_regions=editable_regions,
        df_file=df_file,
        output_html=f"{output_dir}/{filename}_prefilter_report.html",
    )
    # shutil.rmtree(f"{output_dir}/msa{suffix}")
    # reporter.plot_heatmap(
    #     out_file=f"{output_dir}/{filename}_single_points_heatmap.pdf",
    # )
    reporter.generate_prefilter_report()