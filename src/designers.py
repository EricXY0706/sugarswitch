import time
import re
import warnings
import os
import subprocess
from pathlib import Path
import shutil
import yaml
import pandas as pd
import numpy as np

from Bio import SeqIO
import torch
from transformers import EsmTokenizer, EsmForMaskedLM
from peft import PeftModel
import gc

from src.design_utils import set_seed, sample_sites, add_single_mutation
from src.esm_model import EsmModelClassification, ESM_TOKENS
from src.util import *

eps = -1e9
BASE_MODEL_NAME = "facebook/esm2_t30_150M_UR50D"

LORA_MODEL_NAME = "./ESM-LoRA-Gly/checkpoints/N-linked/ESM-150M/checkpoint"
HUMAN_LORA_MODEL_NAME = "./ESM-LoRA-Gly/checkpoints/N-linked/ESM-150M/checkpoint-human"

def _sample_gumbel(shape, device, eps=1e-9):
    u = torch.rand(shape, device=device)
    return -torch.log(-torch.log(u + eps) + eps)

def prepare_seq(
    input_fasta_file: str,
    wt_structure_file: str,
    output_dir: str,
    conservation_df: pd.DataFrame,
    coupling_stength: np.array,
    interaction_dict: dict,
    rsasa_index_dict: dict,
    name: str = None,
    chain_id: str = "A",
    num_gly_sites: int = 3,
    trial_times: list = [0],
    combination_id: int = 0,
):
    filename = name if name else Path(input_fasta_file).name.split(".")[0]
    modify_seq_id = (ord(chain_id) - ord("A") + 1)
    
    seqs = {}
    wt_seq = ""
    seq_to_design = ""
    sampled_sites = []
    wt_sequons = {}
    mut_sequons = {}
    states = {}
    
    count_l = count_r = 0
    for rec in SeqIO.parse(input_fasta_file, "fasta"):
        seq_num = int(rec.description.split("copies:")[1])
        count_r += seq_num
        
        if count_l < modify_seq_id and count_r >= modify_seq_id:
            seq = str(rec.seq)
            wt_seq = seq
            seq_to_design = seq
            sampled_sites = sample_sites(
                structure_file=wt_structure_file,
                scoring_df=f"{output_dir}/{filename}_prefilter_result.csv",
                conservation_df=conservation_df,
                coupling_stength=coupling_stength,
                interaction_dict=interaction_dict,
                rsasa_index_dict=rsasa_index_dict,
                chain_id=chain_id,
                wt_seq=wt_seq,
                num_sites_per_comb=num_gly_sites,
                combination_id=combination_id,
            )
            sampled_sites = sorted(sampled_sites, reverse=True)
            
            for i, s in enumerate(sampled_sites):
                gly_site_seqid = len(sampled_sites) - 1 - i
                seq_to_design, wt_sequon, mut_sequon, state = add_single_mutation(
                    chain_id=chain_id,
                    site=s, 
                    wt_seq=seq_to_design, 
                    conservation_df=conservation_df,
                    coupling_stength=coupling_stength,
                    interaction_dict=interaction_dict,
                    rsasa_index_dict=rsasa_index_dict,
                    N_minus_1_aa=wt_seq[s-2] if s > 1 else "X",
                    N_plus_1_aa=wt_seq[s] if s <= (len(wt_seq) - 1) else "X",
                    N_plus_2_aa=wt_seq[s+1] if s <= (len(wt_seq) - 2) else "X",
                    trial_times=trial_times,
                    gly_site_seqid=gly_site_seqid,
                )
                print(s, wt_sequon, mut_sequon, state, flush=True)
                wt_sequons[s] = wt_sequon
                mut_sequons[s] = mut_sequon
                states[s] = state
            seqs[rec.description] = seq_to_design
        else:
            seq = str(rec.seq)
            seqs[rec.description] = seq        
        count_l += seq_num
        
    asn_sites = {chain_id: [m.start() + 1 for m in re.finditer(r"N[^P][TS]", seq_to_design)]}
    wt_sequons = dict(sorted(wt_sequons.items(), key=lambda x: x[0]))
    mut_sequons = dict(sorted(mut_sequons.items(), key=lambda x: x[0]))
    states = dict(sorted(states.items(), key=lambda x: x[0]))
    
    return seq_to_design, asn_sites, seqs, wt_seq, sampled_sites, wt_sequons, mut_sequons, states
    
def _load_model(
    base_model_name: str,
    lora_model_name: str,
    device: torch.device,
    optional_lora_model_name: Optional[str] = None,
    load_masked_esm: bool = False,
):
    base_model = EsmModelClassification.from_pretrained(
        base_model_name,
        num_labels=2,
        torch_dtype=torch.float16,
    )

    model = PeftModel.from_pretrained(base_model, lora_model_name)
    model = model.merge_and_unload()

    if optional_lora_model_name is not None:
        model = PeftModel.from_pretrained(model, optional_lora_model_name)
        model = model.merge_and_unload()

    model.to(device)
    model.eval()
    
    if load_masked_esm:
        masked_lm = EsmForMaskedLM.from_pretrained(base_model_name, torch_dtype=torch.float16).to(device)
        masked_lm.eval()
        for p in masked_lm.parameters():
            p.requires_grad = False
        
        return model, masked_lm
    
    else:
        return model

def predict_seq(
    model: Any,
    tokenizer: EsmTokenizer,
    sequence: str,
    batch_size: int = 8,
):
    set_seed(seed=int(time.time()))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    candidate_positions_for_model = [m.start() + 1 for m in re.finditer(r"N[^P][ST]", sequence)]
    inputs = tokenizer(sequence, return_tensors="pt")
    
    all_predictions = []
    with torch.no_grad():
        for i in range(0, len(candidate_positions_for_model), batch_size):
            batch_positions = candidate_positions_for_model[i : i + batch_size]
            num_in_batch = len(batch_positions)
            
            batch_input_ids = inputs['input_ids'].repeat(num_in_batch, 1).to(device)
            batch_attention_mask = inputs['attention_mask'].repeat(num_in_batch, 1).to(device)
            batch_pos_tensor = torch.tensor(batch_positions, dtype=torch.long).to(device)
            
            outputs = model(input_ids=batch_input_ids, attention_mask=batch_attention_mask, pos=batch_pos_tensor)
            probs = torch.softmax(outputs.logits, dim=-1)
            predictions = torch.argmax(outputs.logits, dim=-1).cpu().numpy()
            all_predictions.extend(predictions)

    pos_probs = {}
    for i, original_pos in enumerate(candidate_positions_for_model):
        pos_probs[original_pos] = round(probs[i][1].item(), 4)
        
    return pos_probs

def hallucinate(
    model: Any,
    masked_lm: Any,
    tokenizer: EsmTokenizer,
    sequence: str,
    num_steps: int = 200,
    lr: float = 1e-2,
    temperature: float = 1.0,
    add_pll_loss: bool = True,
    pll_weight: float = 1.0,
    device: torch.device = None,
) -> str:
    """
    Hallucinate amino acids only for positions equal to `X`.

    Non-X residues are fixed, including the middle residue in N[^P][ST]
    motifs if it is not X.
    """
    set_seed(seed=int(time.time()))
    sequence = sequence.strip()
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    L = len(sequence)
    A = len(ESM_TOKENS)

    x_positions = {i for i, aa in enumerate(sequence) if aa == "X"}
    opt_positions = sorted(x_positions)

    if len(opt_positions) == 0:
        return sequence

    candidate_positions = [m.start() for m in re.finditer(r"N[^P][ST]", sequence)]
    if len(candidate_positions) == 0:
        return sequence

    pos_ids = [p + 1 for p in candidate_positions]

    natural_aas = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                   "L", "K", "M", "F", "S", "T", "W", "Y", "V"]
    token_keys = list(ESM_TOKENS.keys())
    allowed_idx = [ESM_TOKENS[a] for a in natural_aas if a in ESM_TOKENS]
    K = len(allowed_idx)

    fixed_logits = torch.full((L, A), -5.0, device=device, dtype=torch.float32)
    for i, aa in enumerate(sequence):
        if i not in opt_positions and aa in ESM_TOKENS:
            fixed_logits[i, ESM_TOKENS[aa]] = 5.0

    init_scale = 1.0
    gumbel_init = _sample_gumbel((len(opt_positions), K), device=device) * init_scale
    opt_params = torch.nn.Parameter(gumbel_init.to(dtype=torch.float32))
    optimizer = torch.optim.Adam([opt_params], lr=lr)

    with torch.no_grad():
        embedding_weight = model.esm.embeddings.word_embeddings.weight
        aa_token_ids = torch.tensor(
            [tokenizer._convert_token_to_id(aa) for aa in token_keys],
            device=device,
        )
        E = embedding_weight[aa_token_ids]
        mask_token_id = tokenizer.mask_token_id
        mask_embed = embedding_weight[mask_token_id]

    for step in range(num_steps):
        optimizer.zero_grad()

        seq_logits = fixed_logits.clone()
        for j, seq_pos in enumerate(opt_positions):
            p = opt_params[j]
            full = torch.full((A,), eps, device=device, dtype=p.dtype)
            full[allowed_idx] = p
            seq_logits[seq_pos] = full

        seq_probs = torch.softmax(seq_logits / temperature, dim=-1).to(E.dtype)
        seq_embeds = seq_probs @ E

        cls_id = tokenizer.cls_token_id
        eos_id = tokenizer.eos_token_id
        cls_embed = embedding_weight[cls_id].unsqueeze(0)
        eos_embed = embedding_weight[eos_id].unsqueeze(0)

        inputs_embeds_single = torch.cat(
            [cls_embed, seq_embeds, eos_embed],
            dim=0,
        ).unsqueeze(0)

        batch_size = len(pos_ids)
        inputs_embeds = inputs_embeds_single.repeat(batch_size, 1, 1)
        attention_mask = torch.ones(batch_size, L + 2, device=device)
        pos_tensor = torch.tensor(pos_ids, dtype=torch.long, device=device)

        try:
            model.esm.embeddings.mask_token_id = None
            model.esm.embeddings.token_dropout = False
            masked_lm.esm.embeddings.token_dropout = False
        except Exception:
            pass

        outputs = model(
            inputs_embeds=inputs_embeds,
            attention_mask=attention_mask,
            pos=pos_tensor,
        )
        logits = outputs.logits
        probs = torch.softmax(logits, dim=-1)

        gly_loss = -torch.log(probs[:, 1] + 1e-8).mean()
        loss = gly_loss

        if add_pll_loss and pll_weight > 0.0:
            P = len(opt_positions)
            inputs_embeds_batch = inputs_embeds.repeat(P, 1, 1)

            masked_indices = torch.tensor(
                [pos + 1 for pos in opt_positions],
                device=device,
                dtype=torch.long,
            )

            for b in range(P):
                inputs_embeds_batch[b, masked_indices[b], :] = mask_embed.to(
                    inputs_embeds_batch.dtype
                )

            attention_mask_batch = attention_mask.repeat(P, 1)

            lm_outputs = masked_lm(
                inputs_embeds=inputs_embeds_batch,
                attention_mask=attention_mask_batch,
            )
            lm_logits = lm_outputs.logits

            batch_idx = torch.arange(P, device=device)
            logits_at_mask = lm_logits[batch_idx, masked_indices, :]

            logits_at_mask_reordered = logits_at_mask[:, aa_token_ids]
            lm_log_probs = torch.log_softmax(logits_at_mask_reordered, dim=-1)

            seq_probs_opt = seq_probs[opt_positions, :]
            pll_per_pos = -(seq_probs_opt * lm_log_probs).sum(dim=-1)
            pll_loss = pll_per_pos.mean()

            loss = gly_loss + pll_weight * pll_loss

        loss.backward()
        optimizer.step()

    with torch.no_grad():
        final_seq_logits = fixed_logits.clone()
        for j, seq_pos in enumerate(opt_positions):
            p = opt_params[j]
            full = torch.full((A,), eps, device=device, dtype=p.dtype)
            full[allowed_idx] = p
            final_seq_logits[seq_pos] = full

        final_probs = torch.softmax(final_seq_logits, dim=-1)
        final_idx = torch.argmax(final_probs, dim=-1).cpu().tolist()

    designed_seq = "".join(token_keys[i] for i in final_idx)

    return designed_seq

def _predict_designed_structures(
    design_inputs,
    output_dir: str,
    filename: str,
):
    """
    Run one Boltz batch for all accepted designs.
    """
    if not design_inputs:
        return []

    os.makedirs(output_dir, exist_ok=True)
    infer_globals = getattr(update_infer, "__globals__", {})
    msa_cls = infer_globals.get("MsaFileGenerator", globals().get("MsaFileGenerator"))
    flow_list_cls = infer_globals.get("FlowList", globals().get("FlowList", list))
    chain_ids = infer_globals.get("CHAIN_IDS", globals().get("CHAIN_IDS"))
    if msa_cls is None:
        raise RuntimeError("MsaFileGenerator is required for batch structure prediction.")
    if chain_ids is None:
        chain_ids = [chr(i) for i in range(ord("A"), ord("Z") + 1)]

    batch_name = f"{filename}_designed_batch"
    batch_fasta_dir = Path(output_dir) / f"{batch_name}_fastas"
    batch_input_dir = Path(output_dir) / f"{batch_name}_inputs"
    batch_fasta_dir.mkdir(parents=True, exist_ok=True)
    batch_input_dir.mkdir(parents=True, exist_ok=True)

    for design in design_inputs:
        suffix = design["suffix"]
        target_name = f"{filename}{suffix}"
        fasta_file = batch_fasta_dir / f"{target_name}.fasta"
        msa_dir = Path(output_dir) / f"msa{suffix}"

        with open(fasta_file, "w") as f:
            for seq_des, seq in design["seqs"].items():
                if "X" in seq:
                    f.write(f">{seq_des}\n{design['designed_seq']}\n")
                else:
                    f.write(f">{seq_des}\n{seq}\n")

        msa = msa_cls(input_fasta_file=str(fasta_file))
        with open(fasta_file, "r") as f:
            query = f.read()
        msa.run_mmseqs2(x=query, prefix=str(msa_dir))

        seqs = {
            rec.id: (str(rec.seq), int(rec.description.split(" ")[-1].split(":")[-1]))
            for rec in SeqIO.parse(str(fasta_file), "fasta")
        }
        chain_ptr = 0
        out = []
        for _, (seq, num) in seqs.items():
            ids = flow_list_cls(chain_ids[chain_ptr : chain_ptr + num])
            out.append({
                "protein": {
                    "id": ids,
                    "sequence": seq,
                    "msa": str(msa_dir / f"uniref_{len(out) + 1}.a3m"),
                }
            })
            chain_ptr += num

        with open(batch_input_dir / f"{target_name}.yaml", "w") as f:
            yaml.safe_dump({"version": 1, "sequences": out}, f, sort_keys=False, indent=2)

    cmd = [
        "boltz", "predict", str(batch_input_dir),
        "--cache", f"./boltz_ckpt",
        "--output_format", "pdb",
        "--out_dir", output_dir,
    ]
    subprocess.run(cmd, check=True)

    result_roots = sorted(Path(output_dir).glob(f"boltz_results_{batch_input_dir.name}*"))
    pdb_files = []
    for design in design_inputs:
        target_name = f"{filename}{design['suffix']}"
        candidates = [
            result_root / "predictions" / target_name / f"{target_name}_model_0.pdb"
            for result_root in result_roots
        ]
        pdb_file = next((path for path in candidates if path.exists()), None)
        if pdb_file is None:
            pdb_file = next(
                Path(output_dir).glob(f"boltz_results_*/predictions/{target_name}/{target_name}_model_0.pdb"),
                None,
            )
        if pdb_file is None:
            raise FileNotFoundError(f"Boltz prediction for `{target_name}` was not found.")

        output_pdb = Path(output_dir) / f"{target_name}.pdb"
        shutil.move(str(pdb_file), str(output_pdb))
        pdb_files.append(str(output_pdb))

    for result_root in result_roots:
        shutil.rmtree(result_root, ignore_errors=True)
    for design in design_inputs:
        shutil.rmtree(Path(output_dir) / f"msa{design['suffix']}", ignore_errors=True)
    shutil.rmtree(batch_fasta_dir, ignore_errors=True)
    shutil.rmtree(batch_input_dir, ignore_errors=True)

    return pdb_files

def halludesign_esm(
    input_fasta_file: str,
    output_dir: str,
    conservation_df: pd.DataFrame,
    coupling_stength: np.array,
    interaction_dict: dict,
    rsasa_index_dict: dict,
    name: str,
    chain_id: str = "A",
    num_patterns: int = 1,
    num_candidates_per_pattern: int = 3,
    num_designs_per_pattern: int = 1,
    num_gly_sites: int = 3,
    n_steps: int = 100,
    learning_rate: float = 1e-2,
    temperature: float = 1.0,
    add_pll_loss: bool = True,
    pll_weight: float = 1.0,
    predict_structure: bool = True,
):
    warnings.filterwarnings("ignore")
    filename = name if name else Path(input_fasta_file).name.split(".")[0]
    wt_structure_file = f"{output_dir}/{filename}.pdb"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = EsmTokenizer.from_pretrained(BASE_MODEL_NAME)
    model, masked_llm = _load_model(
        base_model_name=BASE_MODEL_NAME,
        lora_model_name=LORA_MODEL_NAME,
        optional_lora_model_name=HUMAN_LORA_MODEL_NAME,
        device=device,
        load_masked_esm=True,
    )

    designed_seqs_list = []
    pos_probs_list = []
    pos_probs_plot_list = []
    structure_design_inputs = []
    sampled_sites_list = []
    wt_sequons_list = []
    mut_sequons_list = []
    states_list = []
    wt_seq = ""

    pattern_digits = max(1, len(str(num_patterns)))
    design_digits = max(1, len(str(num_designs_per_pattern)))

    for i in range(num_patterns):
        trial_times = [0] * num_gly_sites
        selected_candidates = []
        top_pos_probs = [{0:0.}]
        while not all(v > 0.5 for pos_prob in top_pos_probs for v in pos_prob.values()):
            seq_to_design, asn_sites, seqs, wt_seq, sampled_sites, wt_sequons, mut_sequons, states = prepare_seq(
                input_fasta_file=input_fasta_file,
                wt_structure_file=wt_structure_file,
                output_dir=output_dir,
                conservation_df=conservation_df,
                coupling_stength=coupling_stength,
                interaction_dict=interaction_dict,
                rsasa_index_dict=rsasa_index_dict,
                name=filename,
                chain_id=chain_id,
                num_gly_sites=num_gly_sites,
                trial_times=trial_times,
                combination_id=i,
            )

            natural_asn_sites = [
                m.start() + 1
                for m in re.finditer(r"N[^P][TS]", wt_seq)
            ]

            pattern_candidates = []

            for candidate_id in range(num_candidates_per_pattern):
                designed_seq = hallucinate(
                    model=model,
                    masked_lm=masked_llm,
                    tokenizer=tokenizer,
                    sequence=seq_to_design,
                    num_steps=n_steps,
                    lr=learning_rate,
                    temperature=temperature,
                    add_pll_loss=add_pll_loss,
                    pll_weight=pll_weight,
                    device=device,
                )

                designed_asn_sites = [
                    s
                    for s in [m.start() + 1 for m in re.finditer(r"N[^P][TS]", designed_seq)]
                    if s not in natural_asn_sites
                ]

                gc.collect()
                torch.cuda.empty_cache()
                torch.cuda.reset_max_memory_allocated()

                pos_probs = predict_seq(
                    model=model,
                    tokenizer=tokenizer,
                    sequence=designed_seq,
                    batch_size=8,
                )

                pos_prob_sum = float(np.sum(list(pos_probs.values()))) if pos_probs else 0.0

                pattern_candidates.append({
                    "designed_seq": designed_seq,
                    "designed_asn_sites": designed_asn_sites,
                    "pos_probs": pos_probs,
                    "pos_prob_sum": pos_prob_sum,
                    "candidate_id": candidate_id,
                })

            pattern_candidates = sorted(
                pattern_candidates,
                key=lambda x: x["pos_prob_sum"],
                reverse=True,
            )

            selected_candidates = pattern_candidates[:num_designs_per_pattern]
            top_pos_probs = [candidate["pos_probs"] for candidate in selected_candidates]
            if any(v <= 0.5 for pos_prob in top_pos_probs for v in pos_prob.values()):
                trial_times = [t + 1 if v <= 0.5 else t for t, v in zip(trial_times, top_pos_probs[0].values())]
            else:
                for selected_rank, candidate in enumerate(selected_candidates):
                    designed_seq = candidate["designed_seq"]
                    designed_asn_sites = candidate["designed_asn_sites"]
                    pos_probs = candidate["pos_probs"]

                    designed_seqs_list.append(designed_seq)
                    pos_probs_list.append({
                        s: p
                        for s, p in zip(states.keys(), pos_probs.values())
                    })
                    pos_probs_plot_list.append(pos_probs)
                    sampled_sites_list.append(sampled_sites)
                    wt_sequons_list.append({
                        s: seq
                        for s, seq in zip(states.keys(), wt_sequons.values())
                    })
                    mut_sequons_list.append({
                        s: designed_seq[pos - 2:min(pos - 2 + len(mut_sequon), len(designed_seq))]
                        for s, pos, mut_sequon in zip(states.keys(), designed_asn_sites, mut_sequons.values())
                    })
                    states_list.append(states)

                    if predict_structure:
                        suffix = (
                            f"_pattern_{i + 1:0{pattern_digits}d}"
                            f"_design_{selected_rank + 1:0{design_digits}d}"
                        )
                        structure_design_inputs.append({
                            "suffix": suffix,
                            "seqs": seqs.copy(),
                            "designed_seq": designed_seq,
                        })

    reporter = designer_report(
        input_fasta_file=input_fasta_file,
        wt_seq=wt_seq,
        pos_list=sampled_sites_list,
        designed_seqs=designed_seqs_list,
        pos_probs_list=pos_probs_list,
        pos_probs_plot_list=pos_probs_plot_list,
        wt_sequons_list=wt_sequons_list,
        mut_sequons_list=mut_sequons_list,
        states_list=states_list,
        output_html=f"{output_dir}/{filename}_designer_report.html",
    )
    reporter.generate_designer_report()
    
    if predict_structure:
        _predict_designed_structures(
            design_inputs=structure_design_inputs,
            output_dir=output_dir,
            filename=filename,
        )

    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.reset_max_memory_allocated()
