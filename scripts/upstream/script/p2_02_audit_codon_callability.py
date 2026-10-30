#!/usr/bin/env python3
"""Audit p2_01 mutation absences against Nichols-coordinate, indexed BAMs.

Keep p2_01 outputs unchanged. Run from the project root after p2_01 has finished.
Requires pandas and pysam in the Python environment used to run this script.
"""

import argparse
import collections
import csv
import importlib.util
import json
import re
import sys
from pathlib import Path

import pandas as pd


def load_p2(path):
    spec = importlib.util.spec_from_file_location("project_p2_01", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manifest(paths):
    result = {}
    categories = collections.Counter()
    for raw in Path(paths).read_text().splitlines():
        path = Path(raw.strip())
        if not raw.strip():
            continue
        value = str(path)
        if "/bam/ss14_align2_nichols/postfilter/" in value and path.name.endswith(".ss14_a2n.final.bam"):
            sample = path.name.removesuffix(".ss14_a2n.final.bam")
            kind = "reads"
        elif "/bam/nichols/postfilter/" in value and path.name.endswith(".nichols.final.bam"):
            sample = path.name.removesuffix(".nichols.final.bam")
            kind = "reads"
        elif "/results/04_complete_genome_variant_calls/bam/" in value and path.name.endswith("_vs_nichols.bam"):
            sample = path.name.removesuffix("_vs_nichols.bam")
            kind = "assembly"
        else:
            continue
        if sample in result:
            raise ValueError("More than one selected Nichols BAM for %s: %s and %s" % (sample, result[sample][0], path))
        result[sample] = (path, kind)
        categories[kind] += 1
    return result, categories


def columns_to_codons(columns, p2, ref_seq):
    by_codon = collections.defaultdict(list)
    for column in columns:
        match = re.fullmatch(r"([^:]+):p\.([A-Za-z*]+)([0-9]+)([A-Za-z*]+)", column)
        if not match:
            raise ValueError("Unsupported locus column: " + column)
        gene, ref_three, number, alt_three = match.groups()
        if gene not in p2.RESISTANCE_GENES:
            raise ValueError("Gene not present in Nichols coordinates: " + gene)
        aa_pos = int(number)
        _, start, end, strand, _ = p2.RESISTANCE_GENES[gene]
        if strand == "+":
            lo = start + 3 * (aa_pos - 1)
        else:
            lo = end - 3 * aa_pos + 1
        hi = lo + 2
        if lo < start or hi > end:
            raise ValueError("Codon outside annotated gene: " + column)
        ref_codon = ref_seq[lo - 1:hi]
        if strand == "-":
            ref_codon = p2.revcomp(ref_codon)
        ref_aa = p2.aa3(p2.CODON_TABLE.get(ref_codon, "X"))
        if ref_aa != ref_three:
            raise ValueError("Reference/protein mismatch for %s: reference says %s" % (column, ref_aa))
        key = (gene, aa_pos, lo, hi, strand, ref_codon)
        by_codon[key].append((column, alt_three))
    return dict(by_codon)


def codon_from_read(read, lo, hi, min_bq):
    bases = {}
    insertion = False
    sequence = read.query_sequence
    if not sequence:
        return None
    qualities = read.query_qualities
    ref_cursor = read.reference_start
    query_cursor = 0
    for operation, length in read.cigartuples:
        if operation in (0, 7, 8):
            first = max(ref_cursor, lo - 1)
            last = min(ref_cursor + length, hi)
            for ref_idx in range(first, last):
                query_idx = query_cursor + ref_idx - ref_cursor
                if qualities is not None and qualities[query_idx] < min_bq:
                    return None
                base = sequence[query_idx].upper()
                if base not in "ACGT":
                    return None
                bases[ref_idx] = base
            ref_cursor += length
            query_cursor += length
        elif operation == 1:
            if lo - 1 < ref_cursor < hi:
                insertion = True
            query_cursor += length
        elif operation in (2, 3):
            ref_cursor += length
        elif operation == 4:
            query_cursor += length
        elif operation not in (5, 6):
            return None
        if ref_cursor >= hi:
            break
    if insertion or len(bases) != 3:
        return None
    return "".join(bases[pos] for pos in range(lo - 1, hi))


def call_codon(bam, chrom, key, kind, min_dp, min_mq, min_bq, min_fraction, p2):
    gene, aa_pos, lo, hi, strand, ref_codon = key
    molecule = collections.defaultdict(set)
    for read in bam.fetch(chrom, lo - 1, hi):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
            continue
        if read.mapping_quality < min_mq:
            continue
        codon = codon_from_read(read, lo, hi, min_bq if kind == "reads" else 0)
        if codon is not None:
            molecule[read.query_name].add(codon)
    counts = collections.Counter(next(iter(codons)) for codons in molecule.values() if len(codons) == 1)
    depth = sum(counts.values())
    threshold = min_dp if kind == "reads" else 1
    if depth < threshold:
        return None, depth, 0.0, "low_full_codon_depth"
    top_codon, top_count = counts.most_common(1)[0]
    fraction = top_count / depth
    if fraction < min_fraction or (kind == "assembly" and (depth != 1 or len(counts) != 1)):
        return None, depth, fraction, "mixed_or_ambiguous"
    if strand == "-":
        top_codon = p2.revcomp(top_codon)
    aa = p2.aa3(p2.CODON_TABLE.get(top_codon, "X"))
    if aa == "Xaa":
        return None, depth, fraction, "invalid_translation"
    return aa, depth, fraction, "callable"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--p2-script", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--bam-paths", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--exclude-samples", help="One sample ID per line; e.g. remove_samples.txt")
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--min-mq", type=int, default=20)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-fraction", type=float, default=0.90)
    ap.add_argument("--check-only", action="store_true")
    args = ap.parse_args()
    if args.min_dp < 1 or not 0.5 < args.min_fraction <= 1:
        ap.error("min-dp must be positive; min-fraction must be >0.5 and <=1")
    p2 = load_p2(args.p2_script)
    chrom, seq = p2.read_fasta_first_record(args.ref_fasta)
    if chrom != "NC_021490.2":
        raise ValueError("Expected Nichols NC_021490.2 reference; got " + chrom)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0, dtype="string")
    if matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise ValueError("Duplicate samples or loci in p2 matrix")
    if not matrix.isin(["0", "1"]).all().all():
        raise ValueError("p2 input matrix must contain only 0/1; missing entries require manual review")
    if args.exclude_samples:
        excluded = {line.strip() for line in Path(args.exclude_samples).read_text().splitlines() if line.strip()}
        matrix = matrix.drop(index=list(excluded & set(matrix.index)))
    paths, categories = manifest(args.bam_paths)
    missing = sorted(set(matrix.index) - paths.keys())
    if missing:
        raise ValueError("No uniquely selected Nichols BAM for %d samples, first: %s" % (len(missing), ", ".join(missing[:10])))
    codons = columns_to_codons(matrix.columns, p2, seq)
    print("samples=%d codons=%d loci=%d selected_Nichols_BAMs=%s" % (len(matrix), len(codons), len(matrix.columns), dict(categories)), flush=True)
    if args.check_only:
        print("PASS: manifest, Nichols FASTA, p2 matrix and codon coordinates; BAM headers and indexes not checked", flush=True)
        return
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required; check the cluster Python environment") from exc
    outdir = Path(args.outdir)
    if outdir.resolve() == Path(args.matrix).parent.resolve():
        raise ValueError("Output must be separate from the original p2 directory")
    outdir.mkdir(parents=True, exist_ok=True)
    out_matrix = matrix.copy()
    audit_path = outdir / "codon_callability_audit.tsv"
    totals = collections.Counter()
    with audit_path.open("w", newline="") as handle:
        fields = ["sample", "gene", "aa_pos", "nichols_start", "nichols_end", "source", "bam", "reference_codon", "consensus_aa", "full_codon_depth", "major_fraction", "status", "original_states", "revised_states"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for idx, sample in enumerate(matrix.index, 1):
            path, kind = paths[sample]
            if not path.is_file():
                raise FileNotFoundError("Missing BAM: " + str(path))
            with pysam.AlignmentFile(str(path), "rb") as bam:
                if not bam.has_index():
                    raise ValueError("Missing BAM index for " + str(path))
                if chrom not in bam.references or bam.get_reference_length(chrom) != len(seq):
                    raise ValueError("BAM reference ID/length does not match Nichols FASTA: " + str(path))
                for key, loci in codons.items():
                    aa, depth, fraction, status = call_codon(bam, chrom, key, kind, args.min_dp, args.min_mq, args.min_bq, args.min_fraction, p2)
                    ref_aa = p2.aa3(p2.CODON_TABLE[key[5]])
                    if aa is not None and aa != ref_aa and aa not in {alt for _, alt in loci}:
                        status = "unrepresented_aa_state"
                    before = []
                    after = []
                    for column, alt_aa in loci:
                        old = matrix.at[sample, column]
                        if aa is None or status == "unrepresented_aa_state":
                            revised = pd.NA
                        elif old == "1" and aa != alt_aa:
                            revised = pd.NA
                            status = "p2_bam_discordance"
                        elif old == "0" and aa == alt_aa:
                            revised = pd.NA
                            status = "p2_bam_discordance"
                        else:
                            revised = old
                        out_matrix.at[sample, column] = revised
                        before.append("%s=%s" % (column, old))
                        after.append("%s=%s" % (column, "NA" if pd.isna(revised) else revised))
                        totals["NA" if pd.isna(revised) else str(revised)] += 1
                    writer.writerow(dict(sample=sample, gene=key[0], aa_pos=key[1], nichols_start=key[2], nichols_end=key[3], source=kind, bam=str(path), reference_codon=key[5], consensus_aa=aa or "NA", full_codon_depth=depth, major_fraction="%.4f" % fraction, status=status, original_states=";".join(before), revised_states=";".join(after)))
            if idx % 25 == 0 or idx == len(matrix):
                print("audited %d/%d samples" % (idx, len(matrix)), flush=True)
    out_matrix.to_csv(outdir / "sample_locus_ternary_matrix.tsv", sep="\t", na_rep="NA")
    summary = {"input_matrix": args.matrix, "reference": args.ref_fasta, "min_dp_reads": args.min_dp, "min_mq": args.min_mq, "min_bq_reads": args.min_bq, "min_fraction": args.min_fraction, "sample_count": len(matrix), "codon_count": len(codons), "locus_count": len(matrix.columns), "state_counts": dict(totals), "interpretation": "0=callable and specified alternate amino acid absent; 1=callable and p2 alternate supported; NA=uncallable or discordant. Assembly calls require one unambiguous aligned codon."}
    (outdir / "callability_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("COMPLETE: " + str(outdir), flush=True)


if __name__ == "__main__":
    main()
