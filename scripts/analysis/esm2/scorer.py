import logging
import sys
import time
import pandas as pd
import torch
import torch.nn.functional as F
from transformers import AutoModelForMaskedLM, AutoTokenizer
STANDARD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

def setup_logging(log_file=None, verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=level, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%H:%M:%S', handlers=handlers)
    return logging.getLogger(__name__)

def load_model(model_path, device, logger):
    logger.info(f'Loading ESM2 from: {model_path}')
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(model_path, torch_dtype=torch.bfloat16, local_files_only=True)
    model = model.to(device)
    model.eval()
    logger.info(f'  Model loaded in {time.time() - t0:.1f}s on {device}')
    if device == 'cuda':
        mem = torch.cuda.memory_allocated() / 1000000000.0
        logger.info(f'  GPU memory used: {mem:.1f} GB')
    return (model, tokenizer)

def compute_position_logprobs(model, tokenizer, sequence, position_0based, device):
    idx = position_0based
    masked_seq = sequence[:idx] + tokenizer.mask_token + sequence[idx + 1:]
    inputs = tokenizer(masked_seq, return_tensors='pt').to(device)
    with torch.no_grad():
        logits = model(**inputs).logits
    mask_tok_pos = idx + 1
    log_probs = F.log_softmax(logits[0, mask_tok_pos, :], dim=-1)
    result = {}
    for aa in STANDARD_AAS:
        tok_id = tokenizer.convert_tokens_to_ids(aa)
        if tok_id != tokenizer.unk_token_id:
            result[aa] = log_probs[tok_id].item()
    return result

def compute_llr(log_probs, wt_aa, mut_aa):
    if wt_aa not in log_probs or mut_aa not in log_probs:
        return None
    return log_probs[mut_aa] - log_probs[wt_aa]

def precompute_all_positions(model, tokenizer, protein_name, sequence, device, logger, positions=None):
    seq_len = len(sequence)
    if positions is None:
        positions = list(range(seq_len))
    cache = {}
    n = len(positions)
    logger.info(f'  {protein_name}: computing logprobs for {n} positions...')
    t0 = time.time()
    for i, pos in enumerate(positions):
        if pos >= seq_len:
            continue
        cache[pos] = compute_position_logprobs(model, tokenizer, sequence, pos, device)
        if (i + 1) % 50 == 0 or i == n - 1:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (n - i - 1) / rate if rate > 0 else 0
            logger.info(f'    {i + 1}/{n} positions  ({rate:.1f} pos/s, ETA {eta:.0f}s)')
    logger.info(f'  {protein_name}: done in {time.time() - t0:.1f}s')
    return cache

def run_full_dms_mode(model, tokenizer, sequences, device, logger, proteins_to_scan=None):
    logger.info('\n=== FULL DMS MODE ===')
    rows = []
    if proteins_to_scan is None:
        proteins_to_scan = list(sequences.keys())
    for protein in proteins_to_scan:
        if protein not in sequences:
            logger.warning(f'  {protein}: not in FASTA, skipping')
            continue
        seq = sequences[protein]
        logger.info(f'\n  {protein}: full DMS scan ({len(seq)} positions x 19 AAs = {len(seq) * 19:,} mutations)')
        pos_cache = precompute_all_positions(model, tokenizer, protein, seq, device, logger)
        n_computed = 0
        for pos_0, log_probs in pos_cache.items():
            wt_aa = seq[pos_0]
            pos_1 = pos_0 + 1
            if wt_aa not in STANDARD_AAS:
                continue
            for mut_aa in STANDARD_AAS:
                if mut_aa == wt_aa:
                    continue
                llr = compute_llr(log_probs, wt_aa, mut_aa)
                if llr is None:
                    continue
                rows.append({'mutation_id': f'{protein}_{wt_aa}{pos_1}{mut_aa}', 'protein': protein, 'position': pos_1, 'wt_aa': wt_aa, 'mut_aa': mut_aa, 'llr': round(llr, 6), 'scan_mode': 'full_dms'})
                n_computed += 1
        logger.info(f'  {protein}: {n_computed} mutations computed')
    return pd.DataFrame(rows)
