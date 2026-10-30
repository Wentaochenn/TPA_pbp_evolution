#!/usr/bin/env python3
"""Extract the study PBP codon alignments and estimate site-specific FUBAR support.

The four gene alignments retain all surveyed codons. The output includes the
complete site table and a focal/supported-site table, using posterior P(beta >
alpha) > 0.90 as the study's reporting threshold.
"""

import argparse
import json
import math
import os
import subprocess
import sys


# SS14 reference coordinates (1-based, from GFF)
GENES = {
    "mrcA":      (772313, 774967, "-", 884),   # (start, end, strand, expected_aa)
    "pbp1":      (535746, 537620, "+", 624),
    "pbp2":      (823707, 825665, "+", 652),
    "tp47":      (623741, 625045, "-", 434),
}

# Focal sites to highlight in results (1-based aa positions)
FOCAL_SITES = {
    "mrcA":      [506, 708, 625, 44, 487],
    "pbp1":      [564, 16],
    "pbp2":      [366, 415, 351, 549],
    "tp47":      [394, 107],
}

ALIGNMENT = "results/01_masked_alignment/tpa_ss14_masked.fasta"
TREE      = "results/03_iqtree/tpa_ss14_rooted.treefile"
OUTDIR    = "results/10_selection"


def extract_codon_alignment(aln_file, gene_name, start, end, strand,
                             expected_aa, outdir, max_seqs=None):
    """Extract a gene region from whole-genome alignment, codon-align."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    print(f"\n{'='*60}")
    print(f"Extracting {gene_name}: {start}-{end} ({strand})")
    print(f"{'='*60}")

    gene_len = end - start + 1
    print(f"  Expected CDS length: {gene_len} bp ({gene_len//3} codons, incl. stop)")

    STOP_CODONS = {"TAA", "TAG", "TGA"}

    seqs = {}
    n_total = 0
    n_skipped_ambig = 0
    for rec in SeqIO.parse(aln_file, "fasta"):
        n_total += 1
        # Extract gene region (0-based slicing)
        region = str(rec.seq)[start-1:end]

        if len(region) != gene_len:
            continue

        # Negative strand: reverse complement
        if strand == "-":
            region = str(Seq(region).reverse_complement())

        # Remove terminal stop codon (last 3 bp)
        last_codon = region[-3:].upper()
        if last_codon in STOP_CODONS or "N" in last_codon.upper():
            region = region[:-3]

        # Replace any codon containing N/n with gaps (---)
        # HyPhy cannot handle ambiguous bases in codons
        cleaned = []
        for i in range(0, len(region) - 2, 3):
            codon = region[i:i+3]
            if any(c not in "ACGTacgt" for c in codon):
                cleaned.append("---")
            else:
                cleaned.append(codon.upper())
        region = "".join(cleaned)

        # Check: skip if >50% gaps
        n_gaps = region.count("-")
        if n_gaps > len(region) * 0.5:
            n_skipped_ambig += 1
            continue

        # Verify no internal stop codons remain
        has_internal_stop = False
        for i in range(0, len(region) - 2, 3):
            codon = region[i:i+3]
            if codon in STOP_CODONS:
                has_internal_stop = True
                break
        if has_internal_stop:
            n_skipped_ambig += 1
            continue

        seqs[rec.id] = region

        if max_seqs and len(seqs) >= max_seqs:
            break

    cds_len = expected_aa * 3  # without stop codon
    print(f"  Sequences extracted: {len(seqs)} / {n_total}")
    print(f"  Skipped (>50% gaps or internal stops): {n_skipped_ambig}")
    if seqs:
        first_len = len(list(seqs.values())[0])
        print(f"  Final CDS length: {first_len} bp ({first_len//3} codons, no stop)")
        print(f"  Divisible by 3: {'YES' if first_len % 3 == 0 else 'NO'}")
        first_codon = list(seqs.values())[0][:3]
        print(f"  First codon: {first_codon} ({'start' if first_codon == 'ATG' else 'not ATG'})")

    if len(seqs) < 10:
        print(f"  [ERROR] Too few sequences (<10), skipping")
        return None

    # Write codon alignment
    out_fasta = os.path.join(outdir, f"{gene_name}_codon.fasta")
    with open(out_fasta, "w") as f:
        for sid, seq in seqs.items():
            f.write(f">{sid}\n{seq}\n")

    print(f"  [Output] {out_fasta}")
    return out_fasta


def prune_tree(tree_file, seq_ids, outdir, gene_name):
    """Prune tree to match alignment sample IDs using R/ape."""
    out_tree = os.path.join(outdir, f"{gene_name}_tree.nwk")
    ids_file = os.path.join(outdir, f"{gene_name}_ids.txt")

    with open(ids_file, "w") as f:
        for sid in seq_ids:
            f.write(sid + "\n")

    r_cmd = f"""
    suppressMessages(library(ape))
    tree <- read.tree("{tree_file}")
    keep <- scan("{ids_file}", what="", quiet=TRUE)
    drop <- setdiff(tree$tip.label, keep)
    if (length(drop) > 0) tree <- drop.tip(tree, drop)
    write.tree(tree, "{out_tree}")
    cat("Pruned to", length(tree$tip.label), "tips\\n")
    """

    result = subprocess.run(["Rscript", "-e", r_cmd],
                            capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  Pruned tree: {result.stdout.strip()} -> {out_tree}")
        return out_tree
    else:
        print(f"  [ERROR] Tree pruning failed: {result.stderr[:200]}")
        return None


def run_hyphy(method, aln_file, tree_file, outdir, gene_name):
    """Run the FUBAR site model on a prepared codon alignment."""
    method_lower = method.lower()
    out_json = os.path.join(outdir, f"{gene_name}_{method_lower}.json")

    print(f"\n  Running HyPhy {method} on {gene_name}...")
    cmd = [
        "hyphy", method_lower,
        "--alignment", aln_file,
        "--tree", tree_file,
        "--output", out_json,
        "--code", "Universal",
    ]

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)

    if result.returncode != 0:
        print(f"  [ERROR] {method} failed:")
        print(result.stderr[:500] if result.stderr else "no stderr")
        return None

    if os.path.exists(out_json):
        print(f"  [Output] {out_json}")
        return out_json
    else:
        print(f"  [ERROR] output file not created")
        return None


def parse_fubar_results(json_file, gene_name, focal_sites):
    """Read all covered codons using named probability and Bayes-factor columns."""
    with open(json_file) as handle:
        obj = json.load(handle)
    headers = [h[0] if isinstance(h, list) else h for h in obj['MLE']['headers']]
    required = ['alpha', 'beta', 'Prob[alpha<beta]', 'Prob[alpha>beta]', 'BayesFactor[alpha<beta]']
    if any(headers.count(key) != 1 for key in required):
        raise ValueError('Missing or ambiguous FUBAR column names')
    results = {}
    for partition, rows in obj['MLE']['content'].items():
        coverage = obj['data partitions'][str(partition)]['coverage']
        while len(coverage) == 1 and isinstance(coverage[0], list):
            coverage = coverage[0]
        if len(coverage) != len(rows):
            raise ValueError('FUBAR coverage differs from its result-row count')
        for zero_based, row in zip(coverage, rows):
            site = int(zero_based) + 1
            if site in results or len(row) < len(headers):
                raise ValueError('Duplicate codon or incomplete FUBAR row')
            values = dict(zip(headers, row))
            positive = float(values['Prob[alpha<beta]'])
            negative = float(values['Prob[alpha>beta]'])
            if any(not math.isfinite(v) or not 0 <= v <= 1 for v in (positive, negative)):
                raise ValueError('Invalid FUBAR posterior probability')
            results[site] = dict(gene=gene_name, site=site,
                alpha=float(values['alpha']), beta=float(values['beta']),
                prob_positive=positive, prob_purifying=negative,
                bayes_factor_positive=float(values['BayesFactor[alpha<beta]']),
                supported_positive=positive > .9,
                is_focal=site in focal_sites.get(gene_name, []))
    expected = int(obj['input']['number of sites'])
    if set(results) != set(range(1, expected + 1)):
        raise ValueError('Incomplete FUBAR site coverage')
    return [results[i] for i in range(1, expected + 1)]


def main():
    import pandas as pd
    from Bio import SeqIO
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genes', nargs='+', choices=list(GENES), default=list(GENES))
    parser.add_argument('--max_seqs', type=int)
    parser.add_argument('--alignment', default=ALIGNMENT)
    parser.add_argument('--tree', default=TREE)
    parser.add_argument('--outdir', default=OUTDIR)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    all_sites = []
    for gene_name in args.genes:
        start, end, strand, expected_aa = GENES[gene_name]
        alignment = extract_codon_alignment(args.alignment, gene_name, start, end,
            strand, expected_aa, args.outdir, max_seqs=args.max_seqs)
        if alignment is None:
            raise RuntimeError('Codon alignment was not generated: ' + gene_name)
        ids = [rec.id for rec in SeqIO.parse(alignment, 'fasta')]
        tree = prune_tree(args.tree, ids, args.outdir, gene_name)
        if tree is None:
            raise RuntimeError('Gene tree could not be pruned: ' + gene_name)
        result = run_hyphy('fubar', alignment, tree, args.outdir, gene_name)
        if result is None:
            raise RuntimeError('FUBAR did not complete: ' + gene_name)
        records = parse_fubar_results(result, gene_name, FOCAL_SITES)
        if len(records) != expected_aa:
            raise ValueError('FUBAR codon count differs from the reference: ' + gene_name)
        all_sites.extend(records)
    frame = pd.DataFrame(all_sites)
    frame.to_csv(os.path.join(args.outdir, 'fubar_all_sites.tsv'), sep='\t', index=False)
    frame[frame['is_focal'] | frame['supported_positive']].to_csv(
        os.path.join(args.outdir, 'fubar_selected_sites.tsv'), sep='\t', index=False)
    print('FUBAR codons:', len(frame))
    print('Posterior probability > 0.90:', int(frame['supported_positive'].sum()))


if __name__ == '__main__':
    main()
