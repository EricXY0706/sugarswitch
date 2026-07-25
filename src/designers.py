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
MASKED_LLM_MODEL_NAME = "facebook/esm2_t30_150M_UR50D"

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
    designed_sites = {}
    max_trials = [0] * num_gly_sites
    num_combinations = 0
    
    count_l = count_r = 0
    for rec in SeqIO.parse(input_fasta_file, "fasta"):
        seq_num = int(rec.description.split("copies:")[1])
        count_r += seq_num
        
        if count_l < modify_seq_id and count_r >= modify_seq_id:
            seq = str(rec.seq)
            wt_seq = seq
            seq_to_design = seq
            sampled_sites, num_combinations = sample_sites(
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
                return_num_combinations=True,
            )
            sampled_sites = sorted(sampled_sites, reverse=True)
            
            for i, s in enumerate(sampled_sites):
                gly_site_seqid = len(sampled_sites) - 1 - i
                prev_len = len(seq_to_design)
                window_start = max(s - 2, 0)
                seq_to_design, wt_sequon, mut_sequon, state, n_offset, num_motifs = add_single_mutation(
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
                # Splicing this motif shifts every already-recorded (downstream,
                # higher) designed site by the length change.
                delta = len(seq_to_design) - prev_len
                for k, (npos_k, start_k, len_k) in list(designed_sites.items()):
                    designed_sites[k] = (npos_k + delta, start_k + delta, len_k)
                designed_sites[s] = (window_start + n_offset + 1, window_start, len(mut_sequon))
                max_trials[gly_site_seqid] = num_motifs - 1
                # print(s, wt_sequon, mut_sequon, state, flush=True)
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
    designed_sites = dict(sorted(designed_sites.items(), key=lambda x: x[0]))

    return seq_to_design, asn_sites, seqs, wt_seq, sampled_sites, wt_sequons, mut_sequons, states, designed_sites, max_trials, num_combinations
    
def _load_model(
    base_model_name: str,
    masked_llm_model_name: str,
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
        masked_lm = EsmForMaskedLM.from_pretrained(masked_llm_model_name, torch_dtype=torch.float16).to(device)
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
    tokenizer_gly: EsmTokenizer,
    tokenizer_llm: EsmTokenizer,
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

    tokenizer_gly corresponds to model.
    tokenizer_llm corresponds to masked_lm.
    """
    set_seed(seed=int(time.time()))
    sequence = sequence.strip()

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    L = len(sequence)
    A = len(ESM_TOKENS)
    eps = -20.0

    x_positions = {i for i, aa in enumerate(sequence) if aa == "X"}
    opt_positions = sorted(x_positions)

    if len(opt_positions) == 0:
        return sequence

    candidate_positions = [m.start() for m in re.finditer(r"N[^P][ST]", sequence)]
    if len(candidate_positions) == 0:
        return sequence

    pos_ids = [p + 1 for p in candidate_positions]

    natural_aas = [
        "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
        "L", "K", "M", "F", "S", "T", "W", "Y", "V",
    ]
    token_keys = list(ESM_TOKENS.keys())
    allowed_idx = [ESM_TOKENS[a] for a in natural_aas if a in ESM_TOKENS]
    K = len(allowed_idx)

    fixed_logits = torch.full((L, A), -5.0, device=device, dtype=torch.float32)
    for i, aa in enumerate(sequence):
        if i not in opt_positions and aa in ESM_TOKENS:
            fixed_logits[i, ESM_TOKENS[aa]] = 5.0

    gumbel_init = _sample_gumbel((len(opt_positions), K), device=device)
    opt_params = torch.nn.Parameter(gumbel_init.to(dtype=torch.float32))
    optimizer = torch.optim.Adam([opt_params], lr=lr)

    with torch.no_grad():
        gly_embedding_weight = model.esm.embeddings.word_embeddings.weight
        llm_embedding_weight = masked_lm.esm.embeddings.word_embeddings.weight

        aa_token_ids_gly = torch.tensor(
            [tokenizer_gly._convert_token_to_id(aa) for aa in token_keys],
            device=device,
        )
        aa_token_ids_llm = torch.tensor(
            [tokenizer_llm._convert_token_to_id(aa) for aa in token_keys],
            device=device,
        )

        E_gly = gly_embedding_weight[aa_token_ids_gly]
        E_llm = llm_embedding_weight[aa_token_ids_llm]

        gly_cls_embed = gly_embedding_weight[tokenizer_gly.cls_token_id].unsqueeze(0)
        gly_eos_embed = gly_embedding_weight[tokenizer_gly.eos_token_id].unsqueeze(0)

        llm_cls_embed = llm_embedding_weight[tokenizer_llm.cls_token_id].unsqueeze(0)
        llm_eos_embed = llm_embedding_weight[tokenizer_llm.eos_token_id].unsqueeze(0)
        llm_mask_embed = llm_embedding_weight[tokenizer_llm.mask_token_id]

    try:
        model.esm.embeddings.mask_token_id = None
        model.esm.embeddings.token_dropout = False
        masked_lm.esm.embeddings.token_dropout = False
    except Exception:
        pass

    for step in range(num_steps):
        optimizer.zero_grad()

        seq_logits = fixed_logits.clone()
        for j, seq_pos in enumerate(opt_positions):
            p = opt_params[j]
            full = torch.full((A,), eps, device=device, dtype=p.dtype)
            full[allowed_idx] = p
            seq_logits[seq_pos] = full

        seq_probs = torch.softmax(seq_logits / temperature, dim=-1)

        # gly model branch
        seq_embeds_gly = seq_probs.to(E_gly.dtype) @ E_gly
        inputs_embeds_gly_single = torch.cat(
            [gly_cls_embed, seq_embeds_gly, gly_eos_embed],
            dim=0,
        ).unsqueeze(0)

        gly_batch_size = len(pos_ids)
        inputs_embeds_gly = inputs_embeds_gly_single.repeat(gly_batch_size, 1, 1)
        attention_mask_gly = torch.ones(gly_batch_size, L + 2, device=device)
        pos_tensor = torch.tensor(pos_ids, dtype=torch.long, device=device)

        outputs = model(
            inputs_embeds=inputs_embeds_gly,
            attention_mask=attention_mask_gly,
            pos=pos_tensor,
        )

        probs = torch.softmax(outputs.logits, dim=-1)
        gly_loss = -torch.log(probs[:, 1] + 1e-8).mean()
        loss = gly_loss

        # masked LM PLL branch
        if add_pll_loss and pll_weight > 0.0:
            P = len(opt_positions)

            seq_embeds_llm = seq_probs.to(E_llm.dtype) @ E_llm
            inputs_embeds_llm_single = torch.cat(
                [llm_cls_embed, seq_embeds_llm, llm_eos_embed],
                dim=0,
            ).unsqueeze(0)

            inputs_embeds_llm = inputs_embeds_llm_single.repeat(P, 1, 1)
            attention_mask_llm = torch.ones(P, L + 2, device=device)

            masked_indices = torch.tensor(
                [pos + 1 for pos in opt_positions],
                device=device,
                dtype=torch.long,
            )

            for b in range(P):
                inputs_embeds_llm[b, masked_indices[b], :] = llm_mask_embed.to(
                    inputs_embeds_llm.dtype
                )

            lm_outputs = masked_lm(
                inputs_embeds=inputs_embeds_llm,
                attention_mask=attention_mask_llm,
            )

            lm_logits = lm_outputs.logits
            batch_idx = torch.arange(P, device=device)
            logits_at_mask = lm_logits[batch_idx, masked_indices, :]

            logits_at_mask_reordered = logits_at_mask[:, aa_token_ids_llm]
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
    tokenizer_gly = EsmTokenizer.from_pretrained(BASE_MODEL_NAME)
    tokenizer_llm = EsmTokenizer.from_pretrained(MASKED_LLM_MODEL_NAME)
    model, masked_llm = _load_model(
        base_model_name=BASE_MODEL_NAME,
        masked_llm_model_name=MASKED_LLM_MODEL_NAME,
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

    def _store_design(candidate, pattern_idx, design_rank):
        context = candidate["context"]
        states = context["states"]
        wt_sequons = context["wt_sequons"]
        seqs = context["seqs"]
        sampled_sites = context["sampled_sites"]
        designed_sites = context["designed_sites"]

        designed_seq = candidate["designed_seq"]
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
            s: designed_seq[start:start + mlen]
            for s, (npos, start, mlen) in designed_sites.items()
        })
        states_list.append(states)

        if predict_structure:
            suffix = (
                f"_pattern_{pattern_idx + 1:0{pattern_digits}d}"
                f"_design_{design_rank + 1:0{design_digits}d}"
            )
            structure_design_inputs.append({
                "suffix": suffix,
                "seqs": seqs.copy(),
                "designed_seq": designed_seq,
            })

    for i in range(num_patterns):
        combination_id = i
        trial_times = [0] * num_gly_sites
        accepted_candidates = []
        candidate_pool = []
        while len(accepted_candidates) < num_designs_per_pattern:
            seq_to_design, asn_sites, seqs, wt_seq, sampled_sites, wt_sequons, mut_sequons, states, designed_sites, max_trials, num_combinations = prepare_seq(
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
                combination_id=combination_id,
            )

            designed_positions = [npos for npos, _, _ in designed_sites.values()]

            context = {
                "seqs": seqs,
                "sampled_sites": sampled_sites,
                "wt_sequons": wt_sequons,
                "states": states,
                "designed_sites": designed_sites,
            }

            trial_candidates = []

            for candidate_id in range(num_candidates_per_pattern):
                designed_seq = hallucinate(
                    model=model,
                    masked_lm=masked_llm,
                    tokenizer_gly=tokenizer_gly,
                    tokenizer_llm=tokenizer_llm,
                    sequence=seq_to_design,
                    num_steps=n_steps,
                    lr=learning_rate,
                    temperature=temperature,
                    add_pll_loss=add_pll_loss,
                    pll_weight=pll_weight,
                    device=device,
                )

                gc.collect()
                torch.cuda.empty_cache()
                torch.cuda.reset_max_memory_allocated()

                pos_probs = predict_seq(
                    model=model,
                    tokenizer=tokenizer_gly,
                    sequence=designed_seq,
                    batch_size=8,
                )
                pos_probs = {pos: pos_probs.get(pos, 0.0) for pos in designed_positions}

                pos_prob_sum = float(np.sum(list(pos_probs.values()))) if pos_probs else 0.0

                trial_candidates.append({
                    "designed_seq": designed_seq,
                    "pos_probs": pos_probs,
                    "pos_prob_sum": pos_prob_sum,
                    "candidate_id": candidate_id,
                    "context": context,
                })

            trial_candidates = sorted(
                trial_candidates,
                key=lambda x: x["pos_prob_sum"],
                reverse=True,
            )
            candidate_pool.extend(trial_candidates)

            need = num_designs_per_pattern - len(accepted_candidates)
            selected_candidates = trial_candidates[:need]

            failing_candidate = None
            for candidate in selected_candidates:
                if all(v > 0.5 for v in candidate["pos_probs"].values()):
                    accepted_candidates.append(candidate)
                elif failing_candidate is None:
                    failing_candidate = candidate

            if len(accepted_candidates) >= num_designs_per_pattern:
                break

            if failing_candidate is None:
                continue

            failing = [v <= 0.5 for v in failing_candidate["pos_probs"].values()]
            exhausted = all(
                (not fail) or (t >= mt)
                for t, mt, fail in zip(trial_times, max_trials, failing)
            )
            if not exhausted:
                trial_times = [
                    t + 1 if fail else t
                    for t, fail in zip(trial_times, failing)
                ]
            else:
                accepted_ids = {id(c) for c in accepted_candidates}
                for candidate in sorted(
                    candidate_pool,
                    key=lambda x: x["pos_prob_sum"],
                    reverse=True,
                ):
                    if len(accepted_candidates) >= num_designs_per_pattern:
                        break
                    if id(candidate) in accepted_ids:
                        continue
                    accepted_candidates.append(candidate)
                    accepted_ids.add(id(candidate))
                break

        for design_rank, candidate in enumerate(accepted_candidates):
            _store_design(candidate, i, design_rank)

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
