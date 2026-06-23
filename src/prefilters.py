from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stderr
import multiprocessing as mp
from pathlib import Path
from tqdm import *
import ast
import gc
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
import warnings
warnings.filterwarnings("ignore")

from src.evc_utils import EVC_funcs
from src.rosetta_utils import Rosetta_funcs
from src.util import *
from config import basic_configs, ranker_configs

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))

def _submit_on_gpu(executor, gpu_id, fn, *args):
    
    worker_environment = {
        "CUDA_VISIBLE_DEVICES": str(gpu_id),
        "OMP_NUM_THREADS": "4",
        "MKL_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "4",
        "NUMEXPR_NUM_THREADS": "4",
        "TOKENIZERS_PARALLELISM": "false",
        "TRANSFORMERS_VERBOSITY": "error",
        "HF_HUB_DISABLE_PROGRESS_BARS": "1",
    }
    previous_environment = {
        key: os.environ.get(key) for key in worker_environment
    }
    os.environ.update(worker_environment)
    try:
        return executor.submit(fn, gpu_id, *args)
    finally:
        for key, previous_value in previous_environment.items():
            if previous_value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = previous_value

def _configure_worker_runtime(gpu_id):
    
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    import logging
    import torch
    from transformers.utils import logging as transformers_logging

    transformers_logging.set_verbosity_error()
    logging.getLogger("transformers").setLevel(logging.ERROR)
    logging.getLogger("torch.hub").setLevel(logging.ERROR)
    torch.set_num_threads(4)
    try:
        torch.set_num_interop_threads(4)
    except RuntimeError:
        pass

def _mutation_effects_on_gpu(
    gpu_id,
    site_inputs,
    query_sequence,
    parsed_foldseek_seq,
    output_dir,
    filename,
    structure_file,
    protein_chain_id,
    enable_glycan_grafting,
):
    _configure_worker_runtime(gpu_id)
    with open(os.devnull, "w") as devnull, redirect_stderr(devnull):
        return _calculate_mutation_effects(
            site_inputs,
            query_sequence,
            parsed_foldseek_seq,
            output_dir,
            filename,
            structure_file,
            protein_chain_id,
            enable_glycan_grafting,
        )

def _calculate_mutation_effects(
    site_inputs,
    query_sequence,
    parsed_foldseek_seq,
    output_dir,
    filename,
    structure_file,
    protein_chain_id,
    enable_glycan_grafting,
):
    from src.saprot_utils import SaProt_funcs
    from src.spired_utils import Spired_funcs

    saprot = SaProt_funcs()
    spired = Spired_funcs()
    glycan_rosetta = Rosetta_funcs() if enable_glycan_grafting else None
    sequence_length = len(query_sequence)
    shard_results = []

    for s, ss_tag, mutation, sasa_value, gvp_unbind_score in site_inputs:
        mut_score_s = saprot.mutation_score(
            parsed_foldseek_seq=parsed_foldseek_seq,
            mutations={s: "N"},
        )
        if s + 2 <= sequence_length:
            mut_score_s_next2_S = saprot.mutation_score(
                parsed_foldseek_seq=parsed_foldseek_seq,
                mutations={s: "N", s + 2: "S"},
            )
            mut_score_s_next2_T = saprot.mutation_score(
                parsed_foldseek_seq=parsed_foldseek_seq,
                mutations={s: "N", s + 2: "T"},
            )
        else:
            mut_score_s_next2_S = mut_score_s
            mut_score_s_next2_T = mut_score_s

        ddG_s, dTm_s = spired.get_mutation_effect(
            wt_seq=query_sequence,
            mutations={s: "N"},
        )
        if s + 2 <= sequence_length:
            ddG_next2_S, dTm_next2_S = spired.get_mutation_effect(
                wt_seq=query_sequence,
                mutations={s: "N", s + 2: "S"},
            )
            ddG_next2_T, dTm_next2_T = spired.get_mutation_effect(
                wt_seq=query_sequence,
                mutations={s: "N", s + 2: "T"},
            )
        else:
            ddG_next2_S, dTm_next2_S = ddG_s, dTm_s
            ddG_next2_T, dTm_next2_T = ddG_s, dTm_s

        if enable_glycan_grafting:
            glycoprotein_structure_file = (
                f"{output_dir}/glycans/{filename}_"
                f"{list(query_sequence)[s - 1]}{s}N.pdb"
            )
            glycan_rosetta.mutate(
                structure_file=structure_file,
                output_file=glycoprotein_structure_file,
                mutations={f"{protein_chain_id}{s}": "N"},
            )
            ss_dict = StructureLoader.get_secondary_structure(
                structure_file=glycoprotein_structure_file,
                chain_id=protein_chain_id,
            )
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
        else:
            clash_residues = None
        shard_results.append([
            s, ss_tag, mutation, sasa_value, gvp_unbind_score,
            ddG_s, dTm_s, ddG_next2_S, dTm_next2_S,
            ddG_next2_T, dTm_next2_T, mut_score_s,
            mut_score_s_next2_S, mut_score_s_next2_T, clash_residues,
        ])

    return shard_results

def run_prefilters(
    input_fasta_file: str,
    input_structure_file: str,
    output_dir: str,
    name: str = None,
    suffix: str = "",
    protein_chain_id: str = "A",
    functional_hotspots: list = [],
    gpu_ids: tuple = (0, 1, 2, 3, 4, 5, 6, 7),
    enable_glycan_grafting: bool = True,
) -> None:
    """
    Run prefilters on the input fasta file.
    """
    
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
    
    # ---------------------------------------------------- Filter ----------------------------------------------------
    # 1. Interaction (intra- and inter- chain)
    analyzer = InteractionAnalyzer(structure_file=structure_file, chain_id=protein_chain_id)
    interaction_dict = analyzer.analyze()
    ppi_sites = analyzer.extract_ppi_sites(chain_id=protein_chain_id, result=interaction_dict)
    intrachain_sites = analyzer.extract_intra_interaction_sites(chain_id=protein_chain_id, result=interaction_dict)
    intrachain_sites = set([s for s in intrachain_sites if SS_TAG[ss[(protein_chain_id, s)]][1] != "loop"])
    hotspots_sites = analyzer.extract_sites_interacting_with_hotspots(chain_id=protein_chain_id, hotspots=ast.literal_eval(functional_hotspots), result=interaction_dict)
    
    # 2. Conservation and evolutionary coupling
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
        conserverd_coupling_sites = set()

    # 3. Solvent accessible surface area
    rosetta = Rosetta_funcs()
    sasa_cutoff, low_sasa_sites, sasa_index_dict = rosetta.get_SASA(
        structure_file=structure_file,
        cutoff=basic_configs["sasa_cutoff"],
        chain=protein_chain_id,
    )
    
    # 4. Exclude: Sites with inter- and intra- chain interactions, and highly conserved and/or evo-coupled
    non_editable_regions = ppi_sites | intrachain_sites | hotspots_sites | conserverd_coupling_sites | low_sasa_sites
    editable_regions = set(list(range(1, len(query_sequence)+1))) - non_editable_regions
    
    # ---------------------------------------------------- Ranker ----------------------------------------------------
    gpu_ids = tuple(int(gpu_id) for gpu_id in gpu_ids)
    if not gpu_ids:
        raise ValueError("gpu_ids must contain at least one GPU ID")
    if len(set(gpu_ids)) != len(gpu_ids):
        raise ValueError("gpu_ids must not contain duplicate GPU IDs")

    from src.saprot_utils import SaProt_funcs
    with open(os.devnull, "w") as devnull, redirect_stderr(devnull):
        saprot = SaProt_funcs()
        parsed_foldseek_seq = saprot.parse_foldseek_ss(
            query_seq=query_sequence,
            structure_file=structure_file,
            chain_id=protein_chain_id,
        )
    del saprot
    gc.collect()
    
    import torch
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    from src.gvpbind_utils import gvpbind_predict
    gvpbind_score_by_res = gvpbind_predict(
        structure_file=structure_file,
        chain_id=protein_chain_id,
        out_dir=output_dir,
        filename=filename,
    )

    spawn_context = mp.get_context("spawn")    
    os.makedirs(f"{output_dir}/glycans/", exist_ok=True)
    site_inputs = [
        (
            s,
            SS_TAG[ss[(protein_chain_id, s)]][0],
            f"{query_sequence[s - 1]}{s}N",
            round(sasa_index_dict[s], 3),
            round(1 - gvpbind_score_by_res[s], 3),
        )
        for s in sorted(editable_regions)
    ]

    results = []
    if site_inputs:
        active_gpu_ids = gpu_ids[:min(len(gpu_ids), len(site_inputs))]
        shards = [
            site_inputs[index::len(active_gpu_ids)]
            for index in range(len(active_gpu_ids))
        ]
        executors = []
        futures = {}
        try:
            for gpu_id, shard in zip(active_gpu_ids, shards):
                executor = ProcessPoolExecutor(max_workers=4, mp_context=spawn_context)
                executors.append(executor)
                future = _submit_on_gpu(
                    executor,
                    gpu_id,
                    _mutation_effects_on_gpu,
                    shard,
                    query_sequence,
                    parsed_foldseek_seq,
                    output_dir,
                    filename,
                    structure_file,
                    protein_chain_id,
                    enable_glycan_grafting,
                )
                futures[future] = len(shard)

            with tqdm(total=len(site_inputs), dynamic_ncols=True) as progress:
                for future in as_completed(futures):
                    results.extend(future.result())
                    progress.update(futures[future])
        finally:
            for executor in executors:
                executor.shutdown(wait=True, cancel_futures=True)

    results.sort(key=lambda row: row[0])
    df = pd.DataFrame(results)
    df.columns = ["Site", "SS", "Mutation", "SASA_i", "GVP_unbind_score",
                  "ddG", "dTm", "ddG_NXS", "dTm_NXS", "ddG_NXT", "dTm_NXT", 
                  "MutScore", "MutScore_NXS", "MutScore_NXT", "Clash"
                  ]
    ranker = BordaCount(**ranker_configs)
    df = ranker.compute_score(df)
    df_file = f"{output_dir}/{filename}_prefilter_result.csv"
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
    reporter.generate_prefilter_report()