import argparse
import random
import re
import warnings

from Bio import SeqIO
import numpy as np
import torch
from transformers import EsmTokenizer, EsmForMaskedLM
from peft import PeftModel
import gc

from config import basic_configs
from design_utils import set_seed, sample_sites, GLY_MOTIFS
from esm_model import EsmModelClassification, ESM_TOKENS

eps = -1e9

def prepare_seq(
    input_fasta_file: str,
    wt_structure_file: str,
    output_dir: str,
    num_gly_sites: int = 5,
):
    modify_chain_id = basic_configs["protein_chain_id"]
    sampled_sites = sample_sites(
        structure_file=wt_structure_file,
        scoring_df=f"{output_dir}/{basic_configs["name"]}_single_points.csv",
        chain_id=modify_chain_id,
        num_sites_per_comb=num_gly_sites,
    )
    modify_seq_id = (ord(modify_chain_id) - ord("A") + 1)
    count_l = count_r = 0

    for rec in SeqIO.parse(input_fasta_file, "fasta"):
        seq_num = int(rec.description.split("copies:")[1])
        count_r += seq_num
        
        if count_l < modify_seq_id and count_r >= modify_seq_id:
            sampled_sites = sorted(sampled_sites, reverse=True)
            seq = str(rec.seq)
            for s in sampled_sites:
                seq = [seq[:s-1-2], seq[s-1-2:s-1+3], seq[s-1+3:]]
                seq[1] = random.sample(GLY_MOTIFS, 1)[0]
                seq = "".join(seq)
        else:
            seq = str(rec.seq)
    seq = seq.replace("NX", "NP")
    
    return seq
    
def _load_model(
    base_model_name: str,
    lora_model_name: str,
    device: torch.device
):
    base_model = EsmModelClassification.from_pretrained(base_model_name, num_labels=2, torch_dtype=torch.float16)
    model = PeftModel.from_pretrained(base_model, lora_model_name)
    model = model.merge_and_unload()
    model.to(device)
    model.eval()
    return model

def predict(
    sequence: str,
    base_model_name: str,
    lora_model_name: str,
    batch_size: int = 8,
):
    set_seed(42)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = _load_model(base_model_name, lora_model_name, device)
    tokenizer = EsmTokenizer.from_pretrained(base_model_name)
    candidate_positions_for_model = [m.start() + 1 for m in re.finditer(r"N.[ST]", sequence)]
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

    print("\n--- Predicted Glycosylation Sites ---")
    for i, original_pos in enumerate(candidate_positions_for_model):
        if all_predictions[i] == 1:
            print(f"Positive prediction -> Position: {original_pos:<4}, motif: {sequence[original_pos-1:original_pos+2]}, probability: {probs[i][1].item():.4f}")

def hallucinate(
    sequence: str,
    base_model_name: str,
    lora_model_name: str,
    num_steps: int = 300,
    lr: float = 1e-2,
    temperature: float = 1.0,
    pll_weight: float = 0.1,
    device: torch.device = None,
) -> str:
    """Hallucinate amino-acids for `X` in `sequence` so that classifier predicts glycosylation for candidate positions.

    - Only tokens equal to 'X' are optimized. Other tokens are fixed.
    - Candidate positions are located by the motif `N.[ST]` (same as original design).
    - The model is evaluated with repeated inputs (one batch item per candidate position), providing `pos` indices.

    Returns the designed sequence (string of same length as input).
    """
    set_seed(42)
    sequence = sequence.strip()
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    tokenizer = EsmTokenizer.from_pretrained(base_model_name)
    model = _load_model(base_model_name, lora_model_name, device)

    # load masked LM model for PLL regularization (kept in eval, no grad)
    masked_lm = EsmForMaskedLM.from_pretrained(base_model_name, torch_dtype=torch.float16).to(device)
    masked_lm.eval()
    for p in masked_lm.parameters():
        p.requires_grad = False

    L = len(sequence)
    A = len(ESM_TOKENS)

    # identify positions to optimize: all 'X' residues plus the middle residue of NP[ST] motifs
    x_positions = {i for i, aa in enumerate(sequence) if aa == "X"}
    # find NP[ST] motifs (we look for N followed by any then [ST]?) original used NP[ST]
    motif_matches = [m.start() for m in re.finditer(r"NP[ST]", sequence)]
    p_positions = {m + 1 for m in motif_matches if m + 1 < L}

    opt_positions = sorted(list(x_positions.union(p_positions)))

    if len(opt_positions) == 0:
        print("No positions (X or NP[ST] middle) to optimize; returning input sequence.")
        return sequence

    # allowed natural amino acids (single-letter), exclude 'P'
    natural_aas = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","S","T","W","Y","V"]
    token_keys = list(ESM_TOKENS.keys())
    allowed_idx = [ESM_TOKENS[a] for a in natural_aas if a in ESM_TOKENS]
    K = len(allowed_idx)

    # prepare fixed logits for positions we do NOT optimize
    fixed_logits = torch.full((L, A), -5.0, device=device, dtype=torch.float32)
    for i, aa in enumerate(sequence):
        if i not in opt_positions and aa in ESM_TOKENS:
            fixed_logits[i, ESM_TOKENS[aa]] = 5.0

    # create optimizable parameters only over allowed tokens for each opt position
    opt_params = torch.nn.Parameter(torch.zeros((len(opt_positions), K), device=device, dtype=torch.float32))
    optimizer = torch.optim.Adam([opt_params], lr=lr)

    # prepare embedding lookup for all ESM tokens we enumerated
    with torch.no_grad():
        embedding_weight = model.esm.embeddings.word_embeddings.weight
        aa_token_ids = torch.tensor([tokenizer._convert_token_to_id(aa) for aa in token_keys], device=device)
        E = embedding_weight[aa_token_ids]  # [A, D]
        mask_token_id = tokenizer.mask_token_id
        mask_embed = embedding_weight[mask_token_id]

    # candidate positions for glycosylation (motif N.[ST]) -> positions of N (0-based)
    candidate_positions = [m.start() for m in re.finditer(r"NP[ST]", sequence)]
    if len(candidate_positions) == 0:
        print("No NXS/T motifs found; returning input sequence.")
        return sequence

    pos_ids = [p + 1 for p in candidate_positions]

    for step in range(num_steps):
        optimizer.zero_grad()

        # assemble full logits: copy fixed + place optimized logits (mapped into full A) at opt positions
        seq_logits = fixed_logits.clone()
        # map opt_params (K) into full A logits, disallowed indices get very negative value
        for j, seq_pos in enumerate(opt_positions):
            p = opt_params[j]
            full = torch.full((A,), eps, device=device, dtype=p.dtype)
            full[allowed_idx] = p
            seq_logits[seq_pos] = full

        # compute token probabilities and sequence embeddings
        seq_probs = torch.softmax(seq_logits / temperature, dim=-1).to(E.dtype)  # [L, A]
        seq_embeds = seq_probs @ E  # [L, D]

        # add CLS and EOS embeddings
        cls_id = tokenizer.cls_token_id
        eos_id = tokenizer.eos_token_id
        cls_embed = embedding_weight[cls_id].unsqueeze(0)
        eos_embed = embedding_weight[eos_id].unsqueeze(0)

        inputs_embeds_single = torch.cat([cls_embed, seq_embeds, eos_embed], dim=0).unsqueeze(0)  # [1, L+2, D]

        # batch by candidate positions
        batch_size = len(pos_ids)
        inputs_embeds = inputs_embeds_single.repeat(batch_size, 1, 1)
        attention_mask = torch.ones(batch_size, L + 2, device=device)
        pos_tensor = torch.tensor(pos_ids, dtype=torch.long, device=device)

        # disable masking behavior in embeddings when providing inputs_embeds
        try:
            model.esm.embeddings.mask_token_id = None
            model.esm.embeddings.token_dropout = False
            masked_lm.esm.embeddings.token_dropout = False
        except Exception:
            pass

        outputs = model(inputs_embeds=inputs_embeds, attention_mask=attention_mask, pos=pos_tensor)
        logits = outputs.logits  # [batch_size, num_labels]
        probs = torch.softmax(logits, dim=-1)

        # objective: maximize positive class probability across all candidate positions
        gly_loss = -torch.log(probs[:, 1] + 1e-8).mean()

        # ESM PLL loss computed only over optimized positions to save memory
        if len(opt_positions) > 0 and pll_weight > 0.0:
            P = len(opt_positions)
            # inputs_embeds_single: [1, L+2, D] -> expand to [P, L+2, D]
            inputs_embeds_batch = inputs_embeds.repeat(P, 1, 1)

            # masked positions indices in inputs_embeds (offset by 1 for CLS)
            masked_indices = torch.tensor([pos + 1 for pos in opt_positions], device=device, dtype=torch.long)

            # replace the embedding at each masked index in the batch with mask_embed
            for b in range(P):
                inputs_embeds_batch[b, masked_indices[b], :] = mask_embed.to(inputs_embeds_batch.dtype)

            attention_mask_batch = attention_mask.repeat(P, 1)

            # run masked LM on the smaller batch
            lm_outputs = masked_lm(inputs_embeds=inputs_embeds_batch, attention_mask=attention_mask_batch)
            lm_logits = lm_outputs.logits  # [P, L+2, vocab_size]

            # gather logits at masked positions for each batch element
            batch_idx = torch.arange(P, device=device)
            logits_at_mask = lm_logits[batch_idx, masked_indices, :]

            # reorder logits to match token_keys / aa_token_ids ordering
            logits_at_mask_reordered = logits_at_mask[:, aa_token_ids]
            lm_log_probs = torch.log_softmax(logits_at_mask_reordered, dim=-1)

            # take seq_probs only at optimized positions (in token_keys order)
            seq_probs_opt = seq_probs[opt_positions, :]

            pll_per_pos = - (seq_probs_opt * lm_log_probs).sum(dim=-1)
            pll_loss = pll_per_pos.mean()

            loss = gly_loss + pll_weight * pll_loss

        loss.backward()
        optimizer.step()

        if (step + 1) % 10 == 0 or step == num_steps - 1:
            print(f"step {step + 1} | GLY loss {gly_loss.item():.4f} | PLL loss {pll_loss.item():.4f} | total loss {loss.item():.4f}", flush=True)
            
    # construct final sequence: pick argmax tokens at each position
    with torch.no_grad():
        final_seq_logits = fixed_logits.clone()
        for j, seq_pos in enumerate(opt_positions):
            p = opt_params[j]
            full = torch.full((A,), eps, device=device, dtype=p.dtype)
            full[allowed_idx] = p
            final_seq_logits[seq_pos] = full
        final_probs = torch.softmax(final_seq_logits, dim=-1)
        final_idx = torch.argmax(final_probs, dim=-1).cpu().tolist()

    # map indices back to token characters using token_keys ordering
    designed_seq = "".join(token_keys[i] for i in final_idx)
    print(f"Original sequence:\n{sequence}\nDesigned sequence:\n{designed_seq}", flush=True)

    return designed_seq

def halludesign_esm(
    input_fasta_file: str,
    wt_structure_file: str,
    output_dir: str,
    n_steps: int = 200,
    learning_rate: float = 1e-2,
    temperature: float = 1.0,
):
    warnings.filterwarnings("ignore")
    wt_seq = prepare_seq(
        input_fasta_file=input_fasta_file,
        wt_structure_file=wt_structure_file,
        output_dir=output_dir,
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    designed_seq = hallucinate(
        sequence=wt_seq,
        base_model_name="facebook/esm2_t30_150M_UR50D",
        lora_model_name="../ESM-LoRA-Gly/checkpoints/N-linked/ESM-150M/checkpoint",
        num_steps=n_steps,
        lr=learning_rate,
        temperature=temperature,
        device=device,
    )
    
    predict(
        sequence=designed_seq,
        base_model_name="facebook/esm2_t30_150M_UR50D",
        lora_model_name="../ESM-LoRA-Gly/checkpoints/N-linked/ESM-150M/checkpoint",
        batch_size=8,
    )
    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.reset_max_memory_allocated()