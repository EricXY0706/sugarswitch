import time
import re
import warnings
import tempfile
from pathlib import Path
import shutil

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
    structure_unfav_sites: set,
    sequence_unfav_sites: set,
    low_rsasa_sites: set,
    fav_sites: set,
    name: str = None,
    chain_id: str = "A",
    num_gly_sites: int = 3,
    trial_times: int = 0,
):
    filename = name if name else Path(input_fasta_file).name.split(".")[0]
    sampled_sites = sample_sites(
        structure_file=wt_structure_file,
        scoring_df=f"{output_dir}/{filename}_prefilter_result.csv",
        chain_id=chain_id,
        num_sites_per_comb=num_gly_sites,
    )
    modify_seq_id = (ord(chain_id) - ord("A") + 1)
    count_l = count_r = 0
    
    seqs = {}
    wt_seq = ""
    seq_to_design = ""
    for rec in SeqIO.parse(input_fasta_file, "fasta"):
        seq_num = int(rec.description.split("copies:")[1])
        count_r += seq_num
        
        if count_l < modify_seq_id and count_r >= modify_seq_id:
            sampled_sites = sorted(sampled_sites, reverse=True)
            seq = str(rec.seq)
            wt_seq += seq
            seq_to_design += seq
            for s in sampled_sites:
                seq_to_design = add_single_mutation(
                    site=s, 
                    wt_seq=seq_to_design, 
                    structure_unfav_sites=structure_unfav_sites, 
                    sequence_unfav_sites=sequence_unfav_sites, 
                    low_rsasa_sites=low_rsasa_sites, 
                    fav_sites=fav_sites, 
                    N_minus_1_aa=wt_seq[s-2] if s > 1 else "X",
                    N_plus_1_aa=wt_seq[s] if s <= (len(wt_seq) - 1) else "X",
                    N_plus_2_aa=wt_seq[s+1] if s <= (len(wt_seq) - 2) else "X",
                    trial_times=trial_times,
                )
            seqs[rec.description] = seq_to_design
        else:
            seq = str(rec.seq)
            seqs[rec.description] = seq        
        count_l += seq_num
        
    asn_sites = {chain_id: [m.start() + 1 for m in re.finditer(r"N[^P][TS]", seq_to_design)]}
    
    return seq_to_design, asn_sites, seqs, wt_seq
    
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

    # print("\n--- Predicted Glycosylation Sites ---")
    pos_probs = {}
    for i, original_pos in enumerate(candidate_positions_for_model):
        pos_probs[original_pos] = round(probs[i][1].item(), 4)
        # if all_predictions[i] == 1:
        #     print(f"Positive prediction -> Position: {original_pos:<4}, motif: {sequence[original_pos-1:original_pos+2]}, probability: {probs[i][1].item():.4f}")
        # else:
        #     print(f"Negative prediction -> Position: {original_pos:<4}, motif: {sequence[original_pos-1:original_pos+2]}, probability: {probs[i][1].item():.4f}")
    
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
    pll_weight: float = 0.1,
    device: torch.device = None,
) -> str:
    """Hallucinate amino-acids for `X` in `sequence` so that classifier predicts glycosylation for candidate positions.

    - Only tokens equal to 'X' are optimized. Other tokens are fixed.
    - Candidate positions are located by the motif `N.[ST]` (same as original design).
    - The model is evaluated with repeated inputs (one batch item per candidate position), providing `pos` indices.

    Returns the designed sequence (string of same length as input).
    """
    set_seed(seed=int(time.time()))
    sequence = sequence.strip()
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    L = len(sequence)
    A = len(ESM_TOKENS)

    x_positions = {i for i, aa in enumerate(sequence) if aa == "X"}
    motif_matches = [m.start() for m in re.finditer(r"N[^P][ST]", sequence)]
    p_positions = {m + 1 for m in motif_matches if m + 1 < L}
    opt_positions = sorted(list(x_positions.union(p_positions)))

    if len(opt_positions) == 0:
        print("No positions (X or NX[ST] middle) to optimize; returning input sequence.")
        return sequence

    natural_aas = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","S","T","W","Y","V"]
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
        aa_token_ids = torch.tensor([tokenizer._convert_token_to_id(aa) for aa in token_keys], device=device)
        E = embedding_weight[aa_token_ids]  # [A, D]
        mask_token_id = tokenizer.mask_token_id
        mask_embed = embedding_weight[mask_token_id]

    candidate_positions = [m.start() for m in re.finditer(r"N[^P][ST]", sequence)]
    if len(candidate_positions) == 0:
        print("No NXS/T motifs found; returning input sequence.")
        return sequence

    pos_ids = [p + 1 for p in candidate_positions]

    for step in range(num_steps):
        optimizer.zero_grad()

        seq_logits = fixed_logits.clone()
        for j, seq_pos in enumerate(opt_positions):
            p = opt_params[j]
            full = torch.full((A,), eps, device=device, dtype=p.dtype)
            full[allowed_idx] = p
            seq_logits[seq_pos] = full

        seq_probs = torch.softmax(seq_logits / temperature, dim=-1).to(E.dtype)  # [L, A]
        seq_embeds = seq_probs @ E  # [L, D]

        cls_id = tokenizer.cls_token_id
        eos_id = tokenizer.eos_token_id
        cls_embed = embedding_weight[cls_id].unsqueeze(0)
        eos_embed = embedding_weight[eos_id].unsqueeze(0)

        inputs_embeds_single = torch.cat([cls_embed, seq_embeds, eos_embed], dim=0).unsqueeze(0)  # [1, L+2, D]

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

        outputs = model(inputs_embeds=inputs_embeds, attention_mask=attention_mask, pos=pos_tensor)
        logits = outputs.logits  # [batch_size, num_labels]
        probs = torch.softmax(logits, dim=-1)

        gly_loss = -torch.log(probs[:, 1] + 1e-8).mean()
        
        if add_pll_loss:
            if len(opt_positions) > 0 and pll_weight > 0.0:
                P = len(opt_positions)
                inputs_embeds_batch = inputs_embeds.repeat(P, 1, 1)

                masked_indices = torch.tensor([pos + 1 for pos in opt_positions], device=device, dtype=torch.long)

                for b in range(P):
                    inputs_embeds_batch[b, masked_indices[b], :] = mask_embed.to(inputs_embeds_batch.dtype)

                attention_mask_batch = attention_mask.repeat(P, 1)

                lm_outputs = masked_lm(inputs_embeds=inputs_embeds_batch, attention_mask=attention_mask_batch)
                lm_logits = lm_outputs.logits  # [P, L+2, vocab_size]

                batch_idx = torch.arange(P, device=device)
                logits_at_mask = lm_logits[batch_idx, masked_indices, :]

                logits_at_mask_reordered = logits_at_mask[:, aa_token_ids]
                lm_log_probs = torch.log_softmax(logits_at_mask_reordered, dim=-1)

                seq_probs_opt = seq_probs[opt_positions, :]

                pll_per_pos = - (seq_probs_opt * lm_log_probs).sum(dim=-1)
                pll_loss = pll_per_pos.mean()

                loss = gly_loss + pll_weight * pll_loss
        else:
            loss = gly_loss

        loss.backward()
        optimizer.step()

        # if (step + 1) % 10 == 0 or step == 0:
        #     print(f"step {step + 1} | GLY loss {gly_loss.item():.4f} | PLL loss {pll_loss.item():.4f} | total loss {loss.item():.4f}", flush=True)
            
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
    # print(f"Original sequence:\n{sequence}\nDesigned sequence:\n{designed_seq}", flush=True)

    return designed_seq

def halludesign_esm(
    input_fasta_file: str,
    output_dir: str,
    structure_unfav_sites: set,
    sequence_unfav_sites: set,
    low_rsasa_sites: set,
    fav_sites: set,
    name: str,
    chain_id: str = "A",
    num_designs: int = 1,
    num_gly_sites: int = 3,
    n_steps: int = 100,
    learning_rate: float = 1e-2,
    temperature: float = 1.0,
    add_pll_loss: bool = True,
    pll_weight: float = 0.1,
    predict_structure: bool = True,
):
    warnings.filterwarnings("ignore")
    filename = name if name else Path(input_fasta_file).name.split(".")[0]
    wt_structure_file = f"{output_dir}/{filename}.pdb"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = EsmTokenizer.from_pretrained(BASE_MODEL_NAME)
    model, masked_llm = _load_model(base_model_name=BASE_MODEL_NAME, lora_model_name=LORA_MODEL_NAME, optional_lora_model_name=HUMAN_LORA_MODEL_NAME, device=device, load_masked_esm=True)
    
    num_success = 0
    trial_times = -1
    designed_seqs_list_total = []
    pos_probs_list_total = []
    while num_success < num_designs:
        trial_times += 1
        seq_to_design, asn_sites, seqs, wt_seq = prepare_seq(
            input_fasta_file=input_fasta_file,
            wt_structure_file=wt_structure_file,
            output_dir=output_dir,
            structure_unfav_sites=structure_unfav_sites,
            sequence_unfav_sites=sequence_unfav_sites,
            low_rsasa_sites=low_rsasa_sites,
            fav_sites=fav_sites,
            name=filename,
            chain_id=chain_id,
            num_gly_sites=num_gly_sites,
            trial_times=trial_times,
        )
    
        designed_seqs_list = []
        pos_probs_list = []
    
        for i in range(num_designs):
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
            gc.collect()
            torch.cuda.empty_cache()
            torch.cuda.reset_max_memory_allocated()
            
            pos_probs = predict_seq(
                model=model,
                tokenizer=tokenizer,
                sequence=designed_seq,
                batch_size=8,
            )
            
            if all(v > 0.5 for v in pos_probs.values()):
                designed_seqs_list.append(designed_seq)
                pos_probs_list.append(pos_probs)
                designed_seqs_list_total.extend(designed_seqs_list)
                pos_probs_list_total.extend(pos_probs_list)
                
                if predict_structure:
                    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
                        for seq_des, seq in seqs.items():
                            if "X" in seq:
                                f.write(f">{seq_des}\n{designed_seq}\n")
                            else:
                                f.write(f">{seq_des}\n{seq}\n")    
                        f.flush()
                        suffix = f"_designed_{i+1:0{int(len(str(num_designs)))}d}"
                        glycoprotein_structure_file = update_infer(
                            input_fasta_file=f.name,
                            input_structure_file=None,
                            output_dir=output_dir,
                            filename=filename,
                            suffix=suffix,
                        )
                        shutil.rmtree(f"{output_dir}/msa{suffix}")
                    glycanmover = GlycanMover()
                    glycanmover.move(
                        protein_structure_file=glycoprotein_structure_file,
                        glycan_structure_file="./src/G67828VR.pdb",
                        output_pdb=glycoprotein_structure_file,
                        glycan_positions=asn_sites,
                    )
        num_success = len(designed_seqs_list_total)
        
    reporter = designer_report(
        input_fasta_file=input_fasta_file,
        wt_seq=seq_to_design,
        asn_sites=asn_sites,
        designed_seqs=designed_seqs_list_total,
        pos_probs_list=pos_probs_list_total,
        output_html=f"{output_dir}/{filename}_designer_report.html",
    )
    reporter.generate_designer_report()
    
    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.reset_max_memory_allocated()