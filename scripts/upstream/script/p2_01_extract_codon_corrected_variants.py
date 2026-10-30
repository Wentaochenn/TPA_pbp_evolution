#!/usr/bin/env python3
"""
p2_01_extract_codon_corrected_variants.py

Extract target-gene variants from annotated TSV files and generate codon-corrected amino-acid matrices.

Inputs:
  VCF/ss14_align2_nichols/*.ann.tsv
  VCF/nichols/*.ann.tsv
  results/04_complete_genome_variant_calls/ann_tsv/*.ann.tsv
  refs/TPA_Nichols.fasta

Outputs:
  results/05_extracted_variants/all_target_gene_variants.tsv
  results/05_extracted_variants/sample_locus_binary_matrix.tsv
  results/05_extracted_variants/codon_corrected_amino_acid_changes.tsv

Usage:
  python3 script/p2_01_extract_codon_corrected_variants.py \
    --sample_ann_dir VCF/ss14_align2_nichols \
    --sample_ann_dir VCF/nichols \
    --complete_genome_ann_dir results/04_complete_genome_variant_calls/ann_tsv \
    --outdir results/05_extracted_variants \
    --ref_fasta refs/TPA_Nichols.fasta
"""


import os
import re
import sys
import argparse
import logging
import datetime as dt
from io import StringIO
from collections import defaultdict

import pandas as pd


# Target gene coordinates (Nichols reference NC_021490.2), inclusive 1-based coordinates
RESISTANCE_GENES = {
    "gyrA":      ("TPANIC_RS00025",    4391,    6832, "+", "DNA gyrase subunit A"),
    "murA":      ("TPANIC_RS00145",   36917,   38194, "+", "UDP-N-acetylglucosamine transferase (murA)"),
    "rpoB":      ("TPANIC_RS01220",  246741,  250277, "+", "RNA polymerase subunit beta (rpoB)"),
    "rpoD":      ("TPANIC_RS02405",  529477,  531312, "+", "RNA polymerase sigma factor RpoD"),
    "pbp1":      ("TPANIC_RS02435",  535881,  537755, "+", "penicillin-binding protein 1 (pbp1)"),
    "lytB":      ("TPANIC_RS02685",  591996,  593126, "-", "4-hydroxy-3-methylbut-2-enyl diphosphate reductase (lytB/ispH)"),
    "tp47":      ("TPANIC_RS02830",  623858,  625162, "-", "penicillin-binding Tp47 domain A"),
    "dnaE":      ("TPANIC_RS03315",  732207,  735719, "-", "DNA polymerase III subunit alpha (dnaE)"),
    "mrcA":      ("TPANIC_RS03500",  772433,  775087, "-", "penicillin-binding protein (mrcA)"),
    "pbp2":      ("TPANIC_RS03760",  823827,  825785, "+", "penicillin-binding protein (pbp2)"),
    "gyrB":      ("TPANIC_RS04965", 1096629, 1098542, "-", "DNA gyrase subunit B (gyrB)"),
}

DNA_BASES = set("ACGT")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa", "?": "Xaa", "-": "Del",
}

AA3_TO_1 = {v: k for k, v in AA1_TO_3.items()}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--sample_ann_dir", "--vcf_dir1",
        dest="sample_ann_dirs",
        action="append",
        default=None,
        help="WGS ann.tsv input directory. Can be supplied multiple times."
    )
    p.add_argument(
        "--complete_genome_ann_dir", "--vcf_dir2",
        dest="complete_genome_ann_dirs",
        action="append",
        default=None,
        help="Optional complete-genome ann.tsv input directory. Can be supplied multiple times."
    )
    p.add_argument("--outdir", default="results/05_extracted_variants")
    p.add_argument("--min_dp", type=int, default=5)
    p.add_argument("--ref_fasta", required=True,
                   help="Nichols reference FASTA matching RESISTANCE_GENES coordinates, e.g. NC_021490.2")
    p.add_argument("--allow_missing_dp", action="store_true", default=True,
                   help="Keep rows with missing DP. Default: true, matching the old script behavior.")
    p.add_argument("--no_allow_missing_dp", dest="allow_missing_dp", action="store_false")
    p.add_argument("--a506_check", action="store_true", default=True,
                   help="Print explicit mrcA A506T/A506V/A506M QC. Default: true.")
    p.add_argument("--fail_on_multiple_506_states", action="store_true", default=False,
                   help="Exit with nonzero status if corrected matrix still has multiple 506 states per sample.")
    args = p.parse_args()
    if args.sample_ann_dirs is None:
        args.sample_ann_dirs = ["VCF/ss14_align2_nichols", "VCF/nichols"]
    if args.complete_genome_ann_dirs is None:
        args.complete_genome_ann_dirs = ["results/04_complete_genome_variant_calls/ann_tsv"]
    return args


def setup_logging(outdir):
    os.makedirs(outdir, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(outdir, f"p2_01_extract_codon_corrected_variants.{timestamp}.log")

    logger = logging.getLogger("p2_codon_corrected")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    fh = logging.FileHandler(log_file)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger, log_file


def read_fasta_first_record(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Reference FASTA not found: {path}")
    name = None
    seqs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].split()[0]
                else:
                    # Use first record; ignore additional records.
                    break
            else:
                if name is not None:
                    seqs.append(line)
    if name is None or not seqs:
        raise ValueError(f"No FASTA sequence found in {path}")
    seq = "".join(seqs).upper()
    return name, seq


def revcomp(seq):
    return seq.translate(COMP)[::-1].upper()


def complement_base(base):
    return base.translate(COMP).upper()


def translate_cds(seq):
    aas = []
    seq = seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aas.append(CODON_TABLE.get(codon, "X"))
    return "".join(aas)


def aa3(aa1):
    return AA1_TO_3.get(aa1, "Xaa")


def hgvs_p_from_aa(ref_aa1, pos1, alt_aa1):
    return f"p.{aa3(ref_aa1)}{pos1}{aa3(alt_aa1)}"


def parse_hgvs_p(hgvs_p):
    m = re.match(r"p\.([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3}|\*)", str(hgvs_p))
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    return None, None, None


def get_tsv_files(directory):
    if not os.path.isdir(directory):
        return []
    return sorted([
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".ann.tsv")
    ])


def read_ann_tsv(tsv_file):
    # Skip JVM warning lines from snpEff, e.g. [0.026s][warning]...
    with open(tsv_file, "r") as fh:
        lines = [l for l in fh if not l.startswith("[")]
    return pd.read_csv(StringIO("".join(lines)), sep="\t", low_memory=False)


def extract_from_tsv(tsv_file, genes, min_dp, source_tag, allow_missing_dp, logger):
    sample = os.path.basename(tsv_file).replace(".ann.tsv", "")
    records = []
    try:
        df = read_ann_tsv(tsv_file)
    except Exception as e:
        logger.warning("%s: cannot read %s: %s", sample, tsv_file, e)
        return records

    if df.empty or "POS" not in df.columns:
        logger.warning("%s: empty file or missing POS column, skipped", sample)
        return records

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    dp_col = "GEN[0].DP" if "GEN[0].DP" in df.columns else None
    gt_col = "GEN[0].GT" if "GEN[0].GT" in df.columns else None

    for gene_name, (locus, start, end, strand, product) in genes.items():
        sub = df[(df["POS"] >= start) & (df["POS"] <= end)].copy()
        if sub.empty:
            continue

        if dp_col:
            dp_vals = pd.to_numeric(sub[dp_col], errors="coerce")
            # DP==1 indicates reference-genome FASTA alignment (single copy); skip DP filter.
            if not (dp_vals == 1).all():
                if allow_missing_dp:
                    sub = sub[(dp_vals >= min_dp) | dp_vals.isna()]
                else:
                    sub = sub[(dp_vals >= min_dp)]
        if sub.empty:
            continue

        for _, row in sub.iterrows():
            hgvs_p = str(row.get("ANN[*].HGVS_P", "")).split(",")[0].strip()
            hgvs_c = str(row.get("ANN[*].HGVS_C", "")).split(",")[0].strip()
            effect = str(row.get("ANN[*].EFFECT", "")).split(",")[0].strip()
            impact = str(row.get("ANN[*].IMPACT", "")).split(",")[0].strip()

            try:
                pos_int = int(row["POS"])
            except Exception:
                continue

            records.append({
                "sample":    sample,
                "source":    source_tag,
                "gene":      gene_name,
                "locus_tag": locus,
                "product":   product,
                "strand":    strand,
                "CHROM":     str(row.get("CHROM", "")),
                "POS":       pos_int,
                "REF":       str(row.get("REF", "")).upper(),
                "ALT":       str(row.get("ALT", "")).upper(),
                "EFFECT":    effect,
                "IMPACT":    impact,
                "HGVS_C":    hgvs_c,
                "HGVS_P":    hgvs_p,
                "DP":        str(row.get(dp_col, "")) if dp_col else "",
                "GT":        str(row.get(gt_col, "")) if gt_col else "",
            })
    return records


def extract_ref_cds(ref_seq, gene_name):
    locus, start, end, strand, product = RESISTANCE_GENES[gene_name]
    # start/end are inclusive 1-based genomic coordinates.
    gseq = ref_seq[start-1:end].upper()
    cds = gseq if strand == "+" else revcomp(gseq)
    return cds


def genomic_to_coding_index(pos, gene_name):
    locus, start, end, strand, product = RESISTANCE_GENES[gene_name]
    if strand == "+":
        return pos - start
    return end - pos


def alt_to_coding_base(alt, gene_name):
    strand = RESISTANCE_GENES[gene_name][3]
    alt = alt.upper()
    return alt if strand == "+" else complement_base(alt)


def ref_to_coding_base(ref, gene_name):
    strand = RESISTANCE_GENES[gene_name][3]
    ref = ref.upper()
    return ref if strand == "+" else complement_base(ref)


def make_raw_snpeff_outputs(df_all, all_samples, outdir, logger):
    """Preserve raw snpEff per-record outputs for audit, but do not use as final AA states."""
    ms_all = df_all[df_all["EFFECT"].str.contains("missense", case=False, na=False)].copy()
    valid_ms = ms_all[
        ms_all["HGVS_P"].notna() &
        ~ms_all["HGVS_P"].isin(["", "nan", ".", "None"])
    ].copy()
    valid_ms["locus_aa"] = valid_ms["gene"] + ":" + valid_ms["HGVS_P"]

    all_loci = sorted(valid_ms["locus_aa"].unique())
    mat_loci = pd.DataFrame(0, index=all_samples, columns=all_loci)
    mat_loci.index.name = "sample"
    for _, row in valid_ms.iterrows():
        s, locus = row["sample"], row["locus_aa"]
        if s in mat_loci.index and locus in mat_loci.columns:
            mat_loci.loc[s, locus] = 1
    out = os.path.join(outdir, "raw_snpeff_sample_locus_binary_matrix.tsv")
    mat_loci.to_csv(out, sep="\t")
    logger.info("[OUT] raw snpEff per-record locus matrix: %s (%d loci)", out, len(all_loci))

    hgvs = df_all[
        df_all["HGVS_P"].notna() &
        ~df_all["HGVS_P"].isin(["", "nan", ".", "None"])
    ][["gene", "POS", "REF", "ALT", "HGVS_C", "HGVS_P", "EFFECT", "IMPACT"]].drop_duplicates().sort_values(["gene", "POS"])
    out = os.path.join(outdir, "raw_snpeff_unique_amino_acid_changes.tsv")
    hgvs.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] raw snpEff unique amino-acid changes: %s (%d)", out, len(hgvs))

    return mat_loci


def make_nucleotide_wide(df_all, all_samples, outdir, logger):
    mask = (df_all["REF"].str.len() == 1) & (df_all["ALT"].str.len() == 1)
    snp_df = df_all[mask].copy()

    if snp_df.empty:
        logger.warning("no SNP variants, skipping nucleotide wide matrix")
        return

    snp_df["col_name"] = (
        snp_df["gene"] + "_" +
        snp_df["POS"].astype(str) +
        "(ref=" + snp_df["REF"] + ")"
    )

    locus_ref = {}
    for _, row in snp_df.drop_duplicates("col_name").iterrows():
        locus_ref[row["col_name"]] = row["REF"]

    def sort_key_nuc(col):
        m = re.search(r'_(\d+)\(ref=', col)
        pos = int(m.group(1)) if m else 0
        gene = col[:m.start()] if m else col
        return (gene, pos)

    all_cols = sorted(snp_df["col_name"].unique(), key=sort_key_nuc)

    mat = pd.DataFrame(index=all_samples, columns=all_cols)
    mat.index.name = "sample"
    for col, ref in locus_ref.items():
        mat[col] = ref

    for _, row in snp_df.iterrows():
        s = row["sample"]
        col = row["col_name"]
        if s in mat.index and col in mat.columns:
            mat.loc[s, col] = row["ALT"]

    out = os.path.join(outdir, "nucleotide_wide_matrix.tsv")
    mat.to_csv(out, sep="\t")
    logger.info("[OUT] nucleotide wide matrix: %s (%d samples x %d SNP sites)", out, len(all_samples), len(all_cols))


def detect_snpeff_conflicts(df_all, outdir, logger):
    ms = df_all[
        df_all["EFFECT"].str.contains("missense", case=False, na=False) &
        df_all["HGVS_P"].notna() &
        ~df_all["HGVS_P"].isin(["", "nan", ".", "None"])
    ].copy()

    if ms.empty:
        out = os.path.join(outdir, "snpeff_conflicting_aa_calls.tsv")
        pd.DataFrame().to_csv(out, sep="\t", index=False)
        logger.info("[OUT] snpEff conflicting AA calls: %s (0 rows)", out)
        return pd.DataFrame()

    parsed = ms["HGVS_P"].apply(lambda x: pd.Series(parse_hgvs_p(x), index=["ref_aa_snpeff", "aa_pos_snpeff", "alt_aa_snpeff"]))
    ms = pd.concat([ms, parsed], axis=1)
    ms = ms.dropna(subset=["aa_pos_snpeff", "alt_aa_snpeff"])
    if ms.empty:
        out = os.path.join(outdir, "snpeff_conflicting_aa_calls.tsv")
        pd.DataFrame().to_csv(out, sep="\t", index=False)
        logger.info("[OUT] snpEff conflicting AA calls: %s (0 rows)", out)
        return ms

    ms["aa_pos_snpeff"] = ms["aa_pos_snpeff"].astype(int)
    alt_counts = (
        ms.groupby(["sample", "gene", "aa_pos_snpeff"])["alt_aa_snpeff"]
          .nunique()
          .reset_index(name="n_alt_aa_snpeff")
    )
    conflict_keys = alt_counts[alt_counts["n_alt_aa_snpeff"] > 1]

    if conflict_keys.empty:
        conflict_detail = ms.iloc[0:0].copy()
    else:
        conflict_detail = ms.merge(
            conflict_keys[["sample", "gene", "aa_pos_snpeff", "n_alt_aa_snpeff"]],
            on=["sample", "gene", "aa_pos_snpeff"],
            how="inner"
        ).sort_values(["sample", "gene", "aa_pos_snpeff", "POS"])

    out = os.path.join(outdir, "snpeff_conflicting_aa_calls.tsv")
    conflict_detail.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] snpEff conflicting AA calls: %s (%d rows; %d sample-gene-site conflicts)",
                out, len(conflict_detail), len(conflict_keys))

    return conflict_detail


def build_codon_corrected_calls(df_all, all_samples, ref_seq, outdir, logger):
    """Build final codon-level amino-acid changes from raw SNP records."""
    records = []
    compound_records = []
    skipped_records = []
    ref_mismatch_records = []
    same_index_conflict_records = []

    ref_cds_by_gene = {}
    ref_prot_by_gene = {}
    for gene in RESISTANCE_GENES:
        cds = extract_ref_cds(ref_seq, gene)
        ref_cds_by_gene[gene] = cds
        ref_prot_by_gene[gene] = translate_cds(cds)
        if len(cds) % 3 != 0:
            logger.warning("%s CDS length is not divisible by 3: %d", gene, len(cds))
        logger.info("Reference gene %-10s CDS length=%d aa length=%d strand=%s",
                    gene, len(cds), len(ref_prot_by_gene[gene]), RESISTANCE_GENES[gene][3])

    # Only biallelic SNPs can be safely applied to the CDS here.
    df_all = df_all.copy()
    df_all["POS"] = pd.to_numeric(df_all["POS"], errors="coerce")

    grouped = df_all.groupby(["sample", "gene"], sort=False)

    for sample in all_samples:
        for gene in RESISTANCE_GENES:
            ref_cds = ref_cds_by_gene[gene]
            sample_cds = list(ref_cds)
            variants = grouped.get_group((sample, gene)) if (sample, gene) in grouped.groups else pd.DataFrame()

            applied_by_codon = defaultdict(list)
            applied_by_index = {}

            if not variants.empty:
                for _, row in variants.iterrows():
                    pos = int(row["POS"])
                    ref = str(row["REF"]).upper()
                    alt = str(row["ALT"]).upper()

                    if len(ref) != 1 or len(alt) != 1 or ref not in DNA_BASES or alt not in DNA_BASES:
                        skipped_records.append({
                            "sample": sample, "gene": gene, "POS": pos,
                            "REF": ref, "ALT": alt,
                            "reason": "non_biallelic_SNP_or_non_ACGT",
                            "HGVS_P": row.get("HGVS_P", ""),
                            "HGVS_C": row.get("HGVS_C", ""),
                            "EFFECT": row.get("EFFECT", ""),
                        })
                        continue

                    idx = genomic_to_coding_index(pos, gene)
                    if idx < 0 or idx >= len(ref_cds):
                        skipped_records.append({
                            "sample": sample, "gene": gene, "POS": pos,
                            "REF": ref, "ALT": alt,
                            "reason": "coding_index_out_of_range",
                            "HGVS_P": row.get("HGVS_P", ""),
                            "HGVS_C": row.get("HGVS_C", ""),
                            "EFFECT": row.get("EFFECT", ""),
                        })
                        continue

                    coding_ref_expected = ref_cds[idx]
                    coding_ref_from_vcf = ref_to_coding_base(ref, gene)
                    coding_alt = alt_to_coding_base(alt, gene)

                    if coding_ref_expected != coding_ref_from_vcf:
                        ref_mismatch_records.append({
                            "sample": sample, "gene": gene, "POS": pos,
                            "strand": RESISTANCE_GENES[gene][3],
                            "REF_vcf_plus": ref,
                            "ALT_vcf_plus": alt,
                            "coding_index_0based": idx,
                            "coding_pos_1based": idx + 1,
                            "coding_ref_expected": coding_ref_expected,
                            "coding_ref_from_vcf": coding_ref_from_vcf,
                            "coding_alt": coding_alt,
                            "HGVS_P": row.get("HGVS_P", ""),
                            "HGVS_C": row.get("HGVS_C", ""),
                        })
                        # Still apply the ALT because VCF/snpEff record passed original pipeline;
                        # mismatch is logged for audit.

                    if idx in applied_by_index and applied_by_index[idx] != coding_alt:
                        same_index_conflict_records.append({
                            "sample": sample, "gene": gene, "POS": pos,
                            "coding_index_0based": idx,
                            "previous_alt": applied_by_index[idx],
                            "new_alt": coding_alt,
                            "reason": "different_ALT_at_same_coding_index",
                        })
                        sample_cds[idx] = "N"
                    else:
                        sample_cds[idx] = coding_alt
                        applied_by_index[idx] = coding_alt

                    codon_index0 = idx // 3
                    applied_by_codon[codon_index0].append({
                        "sample": sample,
                        "gene": gene,
                        "POS": pos,
                        "strand": RESISTANCE_GENES[gene][3],
                        "REF_vcf_plus": ref,
                        "ALT_vcf_plus": alt,
                        "coding_index_0based": idx,
                        "coding_pos_1based": idx + 1,
                        "coding_ref": coding_ref_expected,
                        "coding_alt": coding_alt,
                        "HGVS_C_snpeff": row.get("HGVS_C", ""),
                        "HGVS_P_snpeff": row.get("HGVS_P", ""),
                        "DP": row.get("DP", ""),
                        "GT": row.get("GT", ""),
                    })

            sample_cds_str = "".join(sample_cds)
            sample_prot = translate_cds(sample_cds_str)
            ref_prot = ref_prot_by_gene[gene]

            # Report compound codon events (>=2 SNPs in same codon).
            for codon_index0, snps in applied_by_codon.items():
                if len(snps) >= 2:
                    codon_start = codon_index0 * 3
                    ref_codon = ref_cds[codon_start:codon_start+3]
                    sample_codon = sample_cds_str[codon_start:codon_start+3]
                    ref_aa = CODON_TABLE.get(ref_codon, "X")
                    alt_aa = CODON_TABLE.get(sample_codon, "X")
                    compound_records.append({
                        "sample": sample,
                        "gene": gene,
                        "aa_pos": codon_index0 + 1,
                        "ref_codon": ref_codon,
                        "sample_codon": sample_codon,
                        "ref_aa": aa3(ref_aa),
                        "sample_aa": aa3(alt_aa),
                        "HGVS_P_final": hgvs_p_from_aa(ref_aa, codon_index0 + 1, alt_aa),
                        "n_snp_in_codon": len(snps),
                        "genomic_POS_list": ";".join(str(x["POS"]) for x in snps),
                        "vcf_REF_ALT_list": ";".join(f'{x["REF_vcf_plus"]}>{x["ALT_vcf_plus"]}' for x in snps),
                        "coding_changes": ";".join(f'{x["coding_pos_1based"]}:{x["coding_ref"]}>{x["coding_alt"]}' for x in snps),
                        "snpeff_HGVS_P_list": ";".join(str(x["HGVS_P_snpeff"]) for x in snps),
                        "DP_list": ";".join(str(x["DP"]) for x in snps),
                        "GT_list": ";".join(str(x["GT"]) for x in snps),
                    })

            # Compare final protein states.
            min_len = min(len(ref_prot), len(sample_prot))
            for aa_index0 in range(min_len):
                ref_aa = ref_prot[aa_index0]
                alt_aa = sample_prot[aa_index0]
                if ref_aa == alt_aa:
                    continue
                # Avoid calling ambiguous codons as true missense states.
                if alt_aa == "X" or ref_aa == "X":
                    continue
                hgvs_p = hgvs_p_from_aa(ref_aa, aa_index0 + 1, alt_aa)
                codon_start = aa_index0 * 3
                records.append({
                    "sample": sample,
                    "gene": gene,
                    "locus_tag": RESISTANCE_GENES[gene][0],
                    "product": RESISTANCE_GENES[gene][4],
                    "strand": RESISTANCE_GENES[gene][3],
                    "aa_pos": aa_index0 + 1,
                    "ref_aa_1": ref_aa,
                    "alt_aa_1": alt_aa,
                    "ref_aa": aa3(ref_aa),
                    "alt_aa": aa3(alt_aa),
                    "HGVS_P": hgvs_p,
                    "locus_aa": f"{gene}:{hgvs_p}",
                    "ref_codon": ref_cds[codon_start:codon_start+3],
                    "sample_codon": sample_cds_str[codon_start:codon_start+3],
                    "n_snp_in_codon": len(applied_by_codon.get(aa_index0, [])),
                    "genomic_POS_list": ";".join(str(x["POS"]) for x in applied_by_codon.get(aa_index0, [])),
                    "snpeff_HGVS_P_list": ";".join(str(x["HGVS_P_snpeff"]) for x in applied_by_codon.get(aa_index0, [])),
                })

    corrected = pd.DataFrame(records)
    compound = pd.DataFrame(compound_records)
    skipped = pd.DataFrame(skipped_records)
    ref_mismatch = pd.DataFrame(ref_mismatch_records)
    same_index_conflict = pd.DataFrame(same_index_conflict_records)

    out = os.path.join(outdir, "codon_corrected_amino_acid_changes.tsv")
    corrected.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] codon-corrected amino-acid changes: %s (%d rows)", out, len(corrected))

    out = os.path.join(outdir, "compound_codon_events.tsv")
    compound.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] compound codon events: %s (%d rows)", out, len(compound))

    out = os.path.join(outdir, "skipped_non_snp_or_ambiguous_records.tsv")
    skipped.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] skipped non-SNP/ambiguous records: %s (%d rows)", out, len(skipped))

    out = os.path.join(outdir, "reference_mismatch_records.tsv")
    ref_mismatch.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] reference mismatch records: %s (%d rows)", out, len(ref_mismatch))

    out = os.path.join(outdir, "same_index_conflict_records.tsv")
    same_index_conflict.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] same-index conflict records: %s (%d rows)", out, len(same_index_conflict))

    return corrected, compound, skipped, ref_mismatch, same_index_conflict, ref_prot_by_gene


def make_codon_corrected_matrices(corrected, all_samples, outdir, logger):
    genes = list(RESISTANCE_GENES.keys())

    if corrected.empty:
        logger.warning("No codon-corrected amino-acid changes found; writing empty matrices.")
        mat_loci = pd.DataFrame(index=all_samples)
        mat_loci.index.name = "sample"
        mat_loci.to_csv(os.path.join(outdir, "sample_locus_binary_matrix.tsv"), sep="\t")
        return mat_loci, pd.DataFrame(index=all_samples), pd.DataFrame(index=all_samples)

    all_loci = sorted(corrected["locus_aa"].unique())
    mat_loci = pd.DataFrame(0, index=all_samples, columns=all_loci)
    mat_loci.index.name = "sample"
    for _, row in corrected.iterrows():
        s, locus = row["sample"], row["locus_aa"]
        if s in mat_loci.index and locus in mat_loci.columns:
            mat_loci.loc[s, locus] = 1
    out = os.path.join(outdir, "sample_locus_binary_matrix.tsv")
    mat_loci.to_csv(out, sep="\t")
    logger.info("[OUT] CODON-CORRECTED per-locus 0/1 matrix: %s (%d samples x %d loci)",
                out, len(all_samples), len(all_loci))

    # sample x gene binary matrix, any codon-corrected amino-acid change
    mat_gene_bin = pd.DataFrame(0, index=all_samples, columns=genes)
    mat_gene_bin.index.name = "sample"
    for _, row in corrected.iterrows():
        s, g = row["sample"], row["gene"]
        if s in mat_gene_bin.index and g in mat_gene_bin.columns:
            mat_gene_bin.loc[s, g] = 1
    out = os.path.join(outdir, "missense_binary_matrix.tsv")
    mat_gene_bin.to_csv(out, sep="\t")
    logger.info("[OUT] CODON-CORRECTED sample x gene 0/1 matrix: %s", out)

    # sample x gene amino-acid changes, pipe-delimited
    mat_gene_aa = pd.DataFrame("", index=all_samples, columns=genes)
    mat_gene_aa.index.name = "sample"
    for _, row in corrected.sort_values(["sample", "gene", "aa_pos"]).iterrows():
        s, g, aa = row["sample"], row["gene"], row["HGVS_P"]
        if s in mat_gene_aa.index and g in mat_gene_aa.columns:
            cur = mat_gene_aa.loc[s, g]
            if cur == "":
                mat_gene_aa.loc[s, g] = aa
            elif aa not in cur.split("|"):
                mat_gene_aa.loc[s, g] = cur + "|" + aa
    out = os.path.join(outdir, "sample_gene_aa_matrix.tsv")
    mat_gene_aa.to_csv(out, sep="\t")
    logger.info("[OUT] CODON-CORRECTED sample x gene AA-change matrix: %s", out)

    # unique corrected AA changes
    unique = corrected[["gene", "aa_pos", "ref_aa", "alt_aa", "HGVS_P", "ref_codon", "sample_codon", "locus_aa"]].drop_duplicates()
    unique = unique.sort_values(["gene", "aa_pos", "HGVS_P"])
    out = os.path.join(outdir, "unique_amino_acid_changes.tsv")
    unique.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] CODON-CORRECTED unique amino-acid changes: %s (%d)", out, len(unique))

    return mat_loci, mat_gene_bin, mat_gene_aa


def make_aminoacid_wide_corrected(corrected, all_samples, ref_prot_by_gene, outdir, logger):
    """
    Amino-acid wide matrix, corrected final states:
      rows  = samples
      cols  = gene_aaPos(ref=wildtype aa3)
      value = amino acid state at that site (aa3), reference if no corrected change
    """
    if corrected.empty:
        logger.warning("no corrected amino-acid changes, skipping aminoacid_wide_matrix.tsv")
        return

    corrected = corrected.copy()
    corrected["col_name"] = (
        corrected["gene"] + "_" +
        corrected["aa_pos"].astype(int).astype(str) +
        "(ref=" + corrected["ref_aa"] + ")"
    )

    def sort_key_aa(col):
        m = re.search(r'_(\d+)\(ref=', col)
        pos = int(m.group(1)) if m else 0
        gene = col[:m.start()] if m else col
        return (gene, pos)

    all_cols = sorted(corrected["col_name"].unique(), key=sort_key_aa)
    col_ref = corrected.drop_duplicates("col_name").set_index("col_name")["ref_aa"].to_dict()

    mat = pd.DataFrame(index=all_samples, columns=all_cols)
    mat.index.name = "sample"
    for col in all_cols:
        mat[col] = col_ref[col]

    for _, row in corrected.iterrows():
        s = row["sample"]
        col = row["col_name"]
        alt = row["alt_aa"]
        if s in mat.index and col in mat.columns:
            mat.loc[s, col] = alt

    out = os.path.join(outdir, "aminoacid_wide_matrix.tsv")
    mat.to_csv(out, sep="\t")
    logger.info("[OUT] CODON-CORRECTED amino-acid wide matrix: %s (%d samples x %d AA sites)",
                out, len(all_samples), len(all_cols))


def write_qc_summary(df_all, corrected, compound, snpeff_conflicts, skipped, ref_mismatch,
                      same_index_conflict, raw_matrix, corrected_matrix, all_samples,
                      outdir, logger):
    rows = []
    rows.append({"metric": "n_samples", "value": len(all_samples)})
    rows.append({"metric": "n_raw_snpeff_records", "value": len(df_all)})
    rows.append({"metric": "n_raw_snpeff_loci", "value": raw_matrix.shape[1] if raw_matrix is not None else 0})
    rows.append({"metric": "n_codon_corrected_aa_change_records", "value": len(corrected)})
    rows.append({"metric": "n_codon_corrected_loci", "value": corrected_matrix.shape[1] if corrected_matrix is not None else 0})
    rows.append({"metric": "n_compound_codon_events", "value": len(compound)})
    rows.append({"metric": "n_snpeff_conflicting_aa_rows", "value": len(snpeff_conflicts)})
    rows.append({"metric": "n_skipped_non_snp_or_ambiguous_records", "value": len(skipped)})
    rows.append({"metric": "n_reference_mismatch_records", "value": len(ref_mismatch)})
    rows.append({"metric": "n_same_index_conflict_records", "value": len(same_index_conflict)})

    # focal counts if columns exist
    focal_cols = [
        "mrcA:p.Ala506Thr",
        "mrcA:p.Ala506Val",
        "mrcA:p.Ala506Met",
        "mrcA:p.Met625Val",
        "mrcA:p.Gly708Ser",
        "pbp1:p.Pro564Leu",
        "tp47:p.Ser394Arg",
    ]
    for col in focal_cols:
        if corrected_matrix is not None and col in corrected_matrix.columns:
            rows.append({"metric": f"corrected_count__{col}", "value": int(corrected_matrix[col].sum())})
        else:
            rows.append({"metric": f"corrected_count__{col}", "value": "NA"})

    # A506 state exclusivity
    a506_cols = [c for c in ["mrcA:p.Ala506Thr", "mrcA:p.Ala506Val", "mrcA:p.Ala506Met"]
                 if corrected_matrix is not None and c in corrected_matrix.columns]
    if a506_cols:
        multi = (corrected_matrix[a506_cols].sum(axis=1) > 1).sum()
        any_state = (corrected_matrix[a506_cols].sum(axis=1) > 0).sum()
        rows.append({"metric": "corrected_A506_any_T_V_M", "value": int(any_state)})
        rows.append({"metric": "corrected_A506_multiple_states_T_V_M", "value": int(multi)})
        if "mrcA:p.Ala506Thr" in a506_cols and "mrcA:p.Ala506Val" in a506_cols:
            double_tv = ((corrected_matrix["mrcA:p.Ala506Thr"] == 1) &
                         (corrected_matrix["mrcA:p.Ala506Val"] == 1)).sum()
            rows.append({"metric": "corrected_A506T_A506V_double_positive", "value": int(double_tv)})

    summary = pd.DataFrame(rows)
    out = os.path.join(outdir, "codon_correction_qc_summary.tsv")
    summary.to_csv(out, sep="\t", index=False)
    logger.info("[OUT] codon-correction QC summary: %s", out)

    return summary


def print_a506_confirmation(corrected_matrix, compound, logger, fail_on_multiple):
    logger.info("\n================ A506 CODON-LEVEL CONFIRMATION ================")
    cols = ["mrcA:p.Ala506Thr", "mrcA:p.Ala506Val", "mrcA:p.Ala506Met"]
    for col in cols:
        if col in corrected_matrix.columns:
            logger.info("%-24s count = %d", col, int(corrected_matrix[col].sum()))
        else:
            logger.info("%-24s count = 0 (column absent)", col)

    present_cols = [c for c in cols if c in corrected_matrix.columns]
    if present_cols:
        multi_mask = corrected_matrix[present_cols].sum(axis=1) > 1
        n_multi = int(multi_mask.sum())
        logger.info("A506 multiple-state genomes among T/V/M = %d", n_multi)
        if "mrcA:p.Ala506Thr" in present_cols and "mrcA:p.Ala506Val" in present_cols:
            n_tv = int(((corrected_matrix["mrcA:p.Ala506Thr"] == 1) &
                        (corrected_matrix["mrcA:p.Ala506Val"] == 1)).sum())
            logger.info("A506T + A506V double-positive genomes = %d", n_tv)
        if n_multi > 0:
            logger.warning("Corrected matrix still has samples with multiple A506 states:")
            logger.warning("%s", corrected_matrix.loc[multi_mask, present_cols].to_string())
            if fail_on_multiple:
                raise SystemExit("ERROR: multiple A506 states remain after codon correction")

    if not compound.empty:
        a506m = compound[(compound["gene"] == "mrcA") & (compound["aa_pos"] == 506)].copy()
        if not a506m.empty:
            logger.info("\nCompound-codon events at mrcA position 506:")
            cols_show = ["sample", "gene", "aa_pos", "ref_codon", "sample_codon", "HGVS_P_final",
                         "genomic_POS_list", "vcf_REF_ALT_list", "snpeff_HGVS_P_list", "DP_list"]
            logger.info("\n%s", a506m[cols_show].to_string(index=False))
    logger.info("================================================================\n")


def main():
    args = parse_args()
    logger, log_file = setup_logging(args.outdir)

    logger.info("Starting p2_01_extract_codon_corrected_variants.py")
    logger.info("Command: %s", " ".join(sys.argv))
    logger.info("Output directory: %s", args.outdir)
    logger.info("Log file: %s", log_file)
    logger.info("min_dp: %s", args.min_dp)
    logger.info("allow_missing_dp: %s", args.allow_missing_dp)

    ref_name, ref_seq = read_fasta_first_record(args.ref_fasta)
    logger.info("Reference FASTA: %s", args.ref_fasta)
    logger.info("Reference record: %s length=%d", ref_name, len(ref_seq))

    # Basic reference bounds check.
    max_end = max(v[2] for v in RESISTANCE_GENES.values())
    if len(ref_seq) < max_end:
        raise ValueError(f"Reference sequence length {len(ref_seq)} shorter than max gene end {max_end}")

    sample_files_raw = []
    complete_files_raw = []
    total_sample_files = 0
    total_complete_files = 0

    for ann_dir in args.sample_ann_dirs:
        files = get_tsv_files(ann_dir)
        total_sample_files += len(files)
        sample_files_raw.extend((f, "sample") for f in files)
        logger.info("Sample ann.tsv files:          %d from %s", len(files), ann_dir)

    for ann_dir in args.complete_genome_ann_dirs:
        files = get_tsv_files(ann_dir)
        total_complete_files += len(files)
        complete_files_raw.extend((f, "complete_genome") for f in files)
        logger.info("Complete-genome ann.tsv files: %d from %s", len(files), ann_dir)

    # De-duplicate by sample ID. Complete-genome ann.tsv files are evaluated first
    # so complete genome FASTA-derived annotations are used when the same accession
    # is also present in a WGS ann.tsv directory. This keeps the final sample count
    # aligned to the tree while avoiding duplicate sample rows.
    all_files_raw = complete_files_raw + sample_files_raw
    seen_samples = {}
    all_files = []
    duplicate_rows = []
    for f, tag in all_files_raw:
        sample_id = os.path.basename(f).replace(".ann.tsv", "")
        if sample_id in seen_samples:
            duplicate_rows.append((sample_id, seen_samples[sample_id][0], f))
            continue
        seen_samples[sample_id] = (f, tag)
        all_files.append((f, tag))

    logger.info("Sample ann.tsv files total:          %d", total_sample_files)
    logger.info("Complete-genome ann.tsv files total: %d", total_complete_files)
    logger.info("Total ann.tsv files before de-duplication: %d", len(all_files_raw))
    logger.info("Total unique sample IDs after de-duplication: %d", len(all_files))
    logger.info("Target genes: %s", ", ".join(RESISTANCE_GENES.keys()))

    if duplicate_rows:
        logger.warning("Duplicate sample IDs found across input dirs and skipped: %d", len(duplicate_rows))
        logger.warning("First duplicates: %s", ", ".join(row[0] for row in duplicate_rows[:20]))

    if not all_files:
        raise SystemExit("No .ann.tsv files found. Check --sample_ann_dir inputs.")

    all_samples = sorted(seen_samples)

    all_records = []
    for i, (tsv_file, tag) in enumerate(all_files):
        if i % 400 == 0:
            logger.info("processing: %d/%d (%.1f%%)", i, len(all_files), i / len(all_files) * 100)
        all_records.extend(extract_from_tsv(tsv_file, RESISTANCE_GENES, args.min_dp, tag, args.allow_missing_dp, logger))

    logger.info("Extracted %d raw variant records from target genes", len(all_records))
    if not all_records:
        raise SystemExit("No variants found in target genes.")

    df_all = pd.DataFrame(all_records)

    # Output raw variants.
    out1 = os.path.join(args.outdir, "all_target_gene_variants.tsv")
    df_all.to_csv(out1, sep="\t", index=False)
    logger.info("[OUT] all raw target-gene variants: %s (%d records)", out1, len(df_all))

    # Per-gene raw files + raw stats.
    logger.info("\nRaw snpEff per-record gene summary:")
    logger.info("%-15s %13s %10s  %s", "gene", "n_mut_samples", "missense", "top amino-acid changes")
    logger.info("%s", "-" * 75)
    for gene in RESISTANCE_GENES:
        sub = df_all[df_all["gene"] == gene]
        if not sub.empty:
            sub.to_csv(os.path.join(args.outdir, f"{gene}_variants.tsv"), sep="\t", index=False)
        ms = sub[sub["EFFECT"].str.contains("missense", case=False, na=False)] if not sub.empty else sub
        n_s = sub["sample"].nunique() if not sub.empty else 0
        n_m = len(ms)
        top = [x for x in (ms["HGVS_P"].value_counts().head(3).index if not ms.empty else []) if x not in ["nan", ".", ""]]
        logger.info("%-15s %13d %10d  %s", gene, n_s, n_m, ", ".join(top))

    # Raw snpEff matrices and QC conflicts.
    raw_matrix = make_raw_snpeff_outputs(df_all, all_samples, args.outdir, logger)
    snpeff_conflicts = detect_snpeff_conflicts(df_all, args.outdir, logger)

    # Nucleotide wide matrix from raw SNPs for audit.
    make_nucleotide_wide(df_all, all_samples, args.outdir, logger)

    # Build codon-level corrected AA calls.
    corrected, compound, skipped, ref_mismatch, same_index_conflict, ref_prot_by_gene = build_codon_corrected_calls(
        df_all, all_samples, ref_seq, args.outdir, logger
    )

    # Corrected matrices.
    corrected_matrix, mat_gene_bin, mat_gene_aa = make_codon_corrected_matrices(corrected, all_samples, args.outdir, logger)
    make_aminoacid_wide_corrected(corrected, all_samples, ref_prot_by_gene, args.outdir, logger)

    # QC summary and confirmation.
    write_qc_summary(df_all, corrected, compound, snpeff_conflicts, skipped, ref_mismatch,
                     same_index_conflict, raw_matrix, corrected_matrix, all_samples,
                     args.outdir, logger)

    if args.a506_check:
        print_a506_confirmation(corrected_matrix, compound, logger, args.fail_on_multiple_506_states)

    logger.info("Done. Output directory: %s", args.outdir)
    logger.info("Log saved to: %s", log_file)


if __name__ == "__main__":
    main()