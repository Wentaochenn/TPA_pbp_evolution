#!/usr/bin/env bash
# p4_03_prepare_beast2_inputs.sh
#
# Prepare BEAUti2 input files for the selected temporal subset.
#
# Inputs:
#   results/12_treetime/selected_subset/treetime_input_samples.tsv
#   results/12_treetime/selected_subset/strict_clock/outliers.tsv
#   results/02_gubbins/tpa_ss14.filtered_polymorphic_sites.fasta
#   results/01_masked_alignment/tpa_ss14_masked.fasta
#
# Outputs:
#   results/13_beast2_selected_subset/input/final_samples.tsv
#   results/13_beast2_selected_subset/input/beast2_subset_snps_variable_only.fasta
#   results/13_beast2_selected_subset/input/fconst.txt
#   results/13_beast2_selected_subset/input/beast2_input_manifest.tsv
#   results/13_beast2_selected_subset/input/beauti2_configuration_guide.txt
#
# Manual step after this script:
#   Generate BEAST2 XML files in BEAUti2 and save them to results/13_beast2_selected_subset/input/.
#
# Usage:
#   bash SC/p4_03_prepare_beast2_inputs.sh

set -euo pipefail

# -----------------------------
# Config
# -----------------------------

SUB_LIST=${SUB_LIST:-"results/12_treetime/selected_subset/treetime_input_samples.tsv"}

TREETIME_DIR=${TREETIME_DIR:-"results/12_treetime/selected_subset/strict_clock"}
OUTLIER_LIST="${TREETIME_DIR}/outliers.tsv"
TREETIME_CLOCK="${TREETIME_DIR}/molecular_clock.txt"
TREETIME_TIMETREE="${TREETIME_DIR}/timetree.nexus"

GUBBINS_SNP="results/02_gubbins/tpa_ss14.filtered_polymorphic_sites.fasta"
MASKED_FA="results/01_masked_alignment/tpa_ss14_masked.fasta"
NICHOLS_REF="refs/TPA_Nichols.fasta"

OUTDIR=${OUTDIR:-"results/13_beast2_selected_subset/input"}

# Effective genome length after pre-mask + recombination removal.
# Keep this consistent with your upstream Gubbins/masking workflow.
GENOME_LENGTH=1083903

# If 1, stop when any final sample is missing from the SNP or masked alignment.
STRICT_SAMPLE_MATCH=1

mkdir -p "${OUTDIR}"

# Export for embedded Python.
export SUB_LIST
export TREETIME_DIR
export OUTLIER_LIST
export TREETIME_CLOCK
export TREETIME_TIMETREE
export GUBBINS_SNP
export MASKED_FA
export NICHOLS_REF
export OUTDIR
export GENOME_LENGTH
export STRICT_SAMPLE_MATCH

echo "=========================================="
echo "BEAST2 input preparation"
echo "=========================================="
echo "Subsample file: ${SUB_LIST}"
echo "TreeTime dir:   ${TREETIME_DIR}"
echo "Outlier file:   ${OUTLIER_LIST}"
echo "SNP alignment:  ${GUBBINS_SNP}"
echo "Masked FASTA:   ${MASKED_FA}"
echo "Output dir:     ${OUTDIR}"
echo ""

# -----------------------------
# Check required files
# -----------------------------

for f in "${SUB_LIST}" "${GUBBINS_SNP}" "${MASKED_FA}"; do
    if [[ ! -f "${f}" ]]; then
        echo "[ERROR] Missing required file: ${f}" >&2
        exit 1
    fi
done

if [[ ! -f "${OUTLIER_LIST}" ]]; then
    echo "[WARN] TreeTime outlier file not found: ${OUTLIER_LIST}"
    echo "       Continuing with zero outliers."
fi

if [[ ! -f "${TREETIME_CLOCK}" ]]; then
    echo "[WARN] TreeTime molecular clock file not found: ${TREETIME_CLOCK}"
else
    echo "TreeTime molecular clock summary:"
    grep -E "rate|r\^2|Root-Tip|--rate|--r" "${TREETIME_CLOCK}" || true
    echo ""
fi

# -----------------------------
# Main preparation
# -----------------------------

python3 << 'PYEOF'
import os
import sys
import math
from collections import OrderedDict, Counter
import pandas as pd

SUB_LIST = os.environ["SUB_LIST"]
TREETIME_DIR = os.environ["TREETIME_DIR"]
OUTLIER_LIST = os.environ["OUTLIER_LIST"]
TREETIME_CLOCK = os.environ["TREETIME_CLOCK"]
TREETIME_TIMETREE = os.environ["TREETIME_TIMETREE"]
GUBBINS_SNP = os.environ["GUBBINS_SNP"]
MASKED_FA = os.environ["MASKED_FA"]
NICHOLS_REF = os.environ["NICHOLS_REF"]
OUTDIR = os.environ["OUTDIR"]
GENOME_LENGTH = int(os.environ["GENOME_LENGTH"])
STRICT_SAMPLE_MATCH = os.environ.get("STRICT_SAMPLE_MATCH", "1") == "1"

os.makedirs(OUTDIR, exist_ok=True)

DNA = set("ACGT")


def log(msg):
    print(msg, flush=True)


def fail(msg):
    raise SystemExit(f"[ERROR] {msg}")


def read_fasta(path):
    records = OrderedDict()
    name = None
    chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()

                header = line[1:].strip()
                name = header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

    if name is not None:
        records[name] = "".join(chunks).upper()

    if not records:
        fail(f"No FASTA records read from {path}")

    return records


def write_fasta(records, path, width=80):
    with open(path, "w") as out:
        for name, seq in records.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def check_alignment_lengths(records, label):
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        length_counts = Counter(len(seq) for seq in records.values())
        fail(
            f"{label} has unequal sequence lengths: "
            + ", ".join(f"{k}:{v}" for k, v in sorted(length_counts.items()))
        )
    return next(iter(lengths))


def read_outliers(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return set()

    # Try normal TSV with header first.
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
        if df.shape[0] == 0:
            return set()

        lower_cols = {c.lower(): c for c in df.columns}
        for candidate in ["sample", "strain", "taxon", "name", "node"]:
            if candidate in lower_cols:
                col = lower_cols[candidate]
                vals = df[col].dropna().astype(str).tolist()
                return set(v.strip() for v in vals if v.strip())

        # Fallback: use first column.
        col = df.columns[0]
        vals = df[col].dropna().astype(str).tolist()
        return set(v.strip() for v in vals if v.strip())

    except Exception:
        pass

    # Fallback: no header.
    out = set()
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.add(line.split()[0])
    return out


def subset_records(records, sample_order, label):
    missing = [s for s in sample_order if s not in records]

    if missing:
        msg = (
            f"{len(missing)} final samples are missing from {label}. "
            f"First missing: {missing[:10]}"
        )
        if STRICT_SAMPLE_MATCH:
            fail(msg)
        else:
            log("[WARN] " + msg)

    kept = OrderedDict()
    for s in sample_order:
        if s in records:
            kept[s] = records[s]

    return kept, missing


def filter_variable_columns(records):
    """Return variable-only alignment and constant SNP-column counts.

    The input is the SNP alignment after subsetting final samples.
    Some columns can become constant after outlier removal/subsampling.
    These constant columns should not remain in a SNP-only BEAST alignment;
    they are added to fconst instead.
    """
    names = list(records.keys())
    seqs = [records[n] for n in names]
    L = check_alignment_lengths(records, "Subset SNP alignment")

    keep_idx = []
    constant_counts = Counter()
    non_acgt_constant = 0
    all_gap_or_unknown = 0

    for i in range(L):
        col = [seq[i] for seq in seqs]
        bases = [b for b in col if b in DNA]
        unique = sorted(set(bases))

        if len(unique) >= 2:
            keep_idx.append(i)
        elif len(unique) == 1:
            constant_counts[unique[0]] += 1
        else:
            all_gap_or_unknown += 1

    new_records = OrderedDict()
    for name, seq in records.items():
        new_records[name] = "".join(seq[i] for i in keep_idx)

    return new_records, constant_counts, {
        "original_snp_columns": L,
        "variable_snp_columns_kept": len(keep_idx),
        "constant_snp_columns_after_subset": sum(constant_counts.values()),
        "all_gap_or_unknown_snp_columns": all_gap_or_unknown
    }


def count_invariant_full_alignment(records):
    """Count invariant A/C/G/T columns in a full masked alignment subset.

    Columns with gaps/N/ambiguous symbols are ignored unless all valid bases
    are A/C/G/T and identical across all samples. This is conservative.
    """
    names = list(records.keys())
    seqs = [records[n] for n in names]
    L = check_alignment_lengths(records, "Subset masked full alignment")

    inv_counts = Counter()
    variable_cols = 0
    ambiguous_or_gap_cols = 0
    all_unknown_cols = 0

    for i in range(L):
        col = [seq[i] for seq in seqs]
        chars = set(col)

        if all(c in DNA for c in chars) and len(chars) == 1:
            inv_counts[next(iter(chars))] += 1
        else:
            bases = [c for c in col if c in DNA]
            unique_bases = set(bases)

            if len(unique_bases) >= 2:
                variable_cols += 1
            elif len(unique_bases) == 1:
                # Column has one DNA state plus N/gap/ambiguous symbols.
                # Do not count it as strict invariant.
                ambiguous_or_gap_cols += 1
            else:
                all_unknown_cols += 1

    return inv_counts, {
        "masked_alignment_length": L,
        "raw_strict_invariant_total": sum(inv_counts.values()),
        "raw_variable_columns_in_masked_subset": variable_cols,
        "ambiguous_or_gap_single_base_columns": ambiguous_or_gap_cols,
        "all_unknown_columns": all_unknown_cols
    }


def scale_counts_to_total(counts, target_total):
    """Scale A/C/G/T counts to sum exactly to target_total."""
    bases = ["A", "C", "G", "T"]
    raw_total = sum(counts.get(b, 0) for b in bases)

    if target_total < 0:
        fail(f"Target invariant total is negative: {target_total}")

    if target_total == 0:
        return Counter({b: 0 for b in bases})

    if raw_total == 0:
        # Fallback to equal distribution if no composition information exists.
        base_val = target_total // 4
        rem = target_total - base_val * 4
        scaled = Counter({b: base_val for b in bases})
        for b in bases[:rem]:
            scaled[b] += 1
        return scaled

    exact = {b: counts.get(b, 0) / raw_total * target_total for b in bases}
    floored = {b: int(math.floor(exact[b])) for b in bases}
    remainder = target_total - sum(floored.values())

    # Largest remainder method.
    remainders = sorted(
        bases,
        key=lambda b: exact[b] - floored[b],
        reverse=True
    )

    scaled = Counter(floored)
    for b in remainders[:remainder]:
        scaled[b] += 1

    return scaled


def parse_clock_file(path):
    if not os.path.exists(path):
        return {"clock_rate": None, "clock_r2": None}

    rate = None
    r2 = None

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if "--rate:" in line:
                try:
                    rate = float(line.split("--rate:")[1].strip().split()[0])
                except Exception:
                    pass
            if "--r^2:" in line:
                try:
                    r2 = float(line.split("--r^2:")[1].strip().split()[0])
                except Exception:
                    pass

    return {"clock_rate": rate, "clock_r2": r2}


# -----------------------------
# Step 1: final sample list
# -----------------------------

log("")
log("[Step 1] Reading selected subset and TreeTime outliers")

sub = pd.read_csv(SUB_LIST, sep="\t")

required_cols = {"sample", "year"}
missing_cols = required_cols - set(sub.columns)
if missing_cols:
    fail(f"Subsample file missing columns: {sorted(missing_cols)}")

sub["sample"] = sub["sample"].astype(str)
sub["year"] = pd.to_numeric(sub["year"], errors="coerce")
sub = sub.dropna(subset=["sample", "year"]).copy()

outliers = read_outliers(OUTLIER_LIST)

log(f"  Subsample rows: {len(sub)}")
log(f"  Outliers read:  {len(outliers)}")

if outliers:
    log(f"  First outliers: {sorted(outliers)[:10]}")

sub["treetime_outlier"] = sub["sample"].isin(outliers)
outlier_count = int(sub["treetime_outlier"].sum())
final = sub.loc[~sub["treetime_outlier"]].copy()

if final.empty:
    fail("All samples were removed as outliers.")

# Keep original order.
final = final.drop_duplicates(subset=["sample"], keep="first")

final_samples = final["sample"].tolist()

log(f"  Final samples after outlier removal: {len(final_samples)}")

final_samples_tsv = os.path.join(OUTDIR, "final_samples.tsv")
final_ids_txt = os.path.join(OUTDIR, "final_sample_ids.txt")

final.to_csv(final_samples_tsv, sep="\t", index=False)

with open(final_ids_txt, "w") as out:
    for s in final_samples:
        out.write(s + "\n")

log(f"  [Output] {final_samples_tsv}")
log(f"  [Output] {final_ids_txt}")

# -----------------------------
# Step 2: dates files
# -----------------------------

log("")
log("[Step 2] Writing BEAST2 date files")

date_simple = os.path.join(OUTDIR, "beast2_dates.tsv")
date_traits = os.path.join(OUTDIR, "beast2_date_traits.tsv")
date_traits_csv = os.path.join(OUTDIR, "beast2_date_traits.csv")

final[["sample", "year"]].to_csv(
    date_simple,
    sep="\t",
    index=False,
    header=True
)

trait_df = pd.DataFrame({
    "taxon": final["sample"],
    "date": final["year"]
})

trait_df.to_csv(
    date_traits,
    sep="\t",
    index=False,
    header=True
)

trait_df.to_csv(
    date_traits_csv,
    index=False,
    header=True
)

log(f"  [Output] {date_simple}")
log(f"  [Output] {date_traits}")
log(f"  [Output] {date_traits_csv}")

# -----------------------------
# Step 3: subset SNP alignment
# -----------------------------

log("")
log("[Step 3] Extracting SNP alignment for final samples")

snp_records_all = read_fasta(GUBBINS_SNP)
snp_subset_all, missing_snp = subset_records(
    snp_records_all,
    final_samples,
    "Gubbins SNP alignment"
)

snp_all_out = os.path.join(OUTDIR, "beast2_subset_snps_all_columns.fasta")
write_fasta(snp_subset_all, snp_all_out)

snp_all_len = check_alignment_lengths(snp_subset_all, "Subset SNP all-column alignment")

log(f"  Records written: {len(snp_subset_all)}")
log(f"  SNP columns before final variable filtering: {snp_all_len}")
log(f"  [Output] {snp_all_out}")

snp_variable_records, constant_snp_counts, snp_filter_stats = filter_variable_columns(snp_subset_all)

snp_var_out = os.path.join(OUTDIR, "beast2_subset_snps_variable_only.fasta")
write_fasta(snp_variable_records, snp_var_out)

snp_var_len = check_alignment_lengths(
    snp_variable_records,
    "Subset SNP variable-only alignment"
)

log(f"  Variable SNP columns kept: {snp_var_len}")
log(f"  Constant SNP columns after subset: {sum(constant_snp_counts.values())}")
log(f"  [Output] {snp_var_out}")

# -----------------------------
# Step 4: subset masked full alignment
# -----------------------------

log("")
log("[Step 4] Extracting masked full alignment for fconst composition")

masked_records_all = read_fasta(MASKED_FA)
masked_subset, missing_masked = subset_records(
    masked_records_all,
    final_samples,
    "masked full alignment"
)

masked_out = os.path.join(OUTDIR, "beast2_subset_masked_full.fasta")
write_fasta(masked_subset, masked_out)

masked_len = check_alignment_lengths(masked_subset, "Subset masked full alignment")

log(f"  Records written: {len(masked_subset)}")
log(f"  Masked alignment length: {masked_len}")
log(f"  [Output] {masked_out}")

# -----------------------------
# Step 5: fconst calculation
# -----------------------------

log("")
log("[Step 5] Computing fconst for BEAST2 ascertainment correction")

raw_inv_counts, masked_stats = count_invariant_full_alignment(masked_subset)

variable_sites_for_beast = snp_var_len
target_invariant_total = GENOME_LENGTH - variable_sites_for_beast

if target_invariant_total < 0:
    fail(
        f"GENOME_LENGTH ({GENOME_LENGTH}) is smaller than variable SNP sites "
        f"({variable_sites_for_beast}). Check GENOME_LENGTH."
    )

scaled_fconst = scale_counts_to_total(raw_inv_counts, target_invariant_total)

# Add diagnostic constant SNP columns after subsetting.
# The final scaled fconst already targets GENOME_LENGTH - variable_sites,
# so these constant SNP columns are represented in the target total.
# We still report their base composition separately for transparency.

fconst_tsv = os.path.join(OUTDIR, "fconst.tsv")
fconst_txt = os.path.join(OUTDIR, "fconst.txt")

bases = ["A", "C", "G", "T"]

fconst_rows = []
for b in bases:
    fconst_rows.append({
        "base": b,
        "raw_strict_invariant_count_from_masked_subset": raw_inv_counts.get(b, 0),
        "constant_snp_columns_after_subset": constant_snp_counts.get(b, 0),
        "scaled_fconst_for_beast2": scaled_fconst.get(b, 0)
    })

fconst_df = pd.DataFrame(fconst_rows)
fconst_df.to_csv(fconst_tsv, sep="\t", index=False)

with open(fconst_txt, "w") as out:
    out.write(
        " ".join(str(scaled_fconst.get(b, 0)) for b in bases) + "\n"
    )

log(f"  GENOME_LENGTH: {GENOME_LENGTH}")
log(f"  Variable SNP sites used for BEAST2: {variable_sites_for_beast}")
log(f"  Target invariant total: {target_invariant_total}")
log(
    "  Raw strict invariant counts from masked subset: "
    + " ".join(f"{b}={raw_inv_counts.get(b, 0)}" for b in bases)
)
log(
    "  Scaled fconst for BEAST2 A C G T: "
    + " ".join(str(scaled_fconst.get(b, 0)) for b in bases)
)
log(f"  [Output] {fconst_tsv}")
log(f"  [Output] {fconst_txt}")

# -----------------------------
# Step 6: BEAUti2 guide and summary
# -----------------------------

log("")
log("[Step 6] Writing BEAUti2 guide and summary")

clock = parse_clock_file(TREETIME_CLOCK)

guide_file = os.path.join(OUTDIR, "beauti2_configuration_guide.txt")
summary_file = os.path.join(OUTDIR, "beast2_input_summary.txt")
manifest_file = os.path.join(OUTDIR, "beast2_input_manifest.tsv")

fconst_string = " ".join(str(scaled_fconst.get(b, 0)) for b in bases)
final_n = len(final_samples)
variable_sites = snp_filter_stats["variable_snp_columns_kept"]
xml_strict_path = os.path.join(OUTDIR, f"tpa_ss14_selected_subset_n{final_n}_strict_{variable_sites}sites_ascertained.xml")
xml_relaxed_path = os.path.join(OUTDIR, f"tpa_ss14_selected_subset_n{final_n}_relaxed_{variable_sites}sites_ascertained.xml")

with open(guide_file, "w") as out:
    out.write("BEAST2 / BEAUti2 configuration guide for selected temporal subset\n")
    out.write("================================================================\n\n")

    out.write("Input alignment:\n")
    out.write(f"  {snp_var_out}\n\n")

    out.write("Date trait file:\n")
    out.write(f"  {date_traits}\n")
    out.write(f"  {date_traits_csv}\n\n")

    out.write("Final sample list:\n")
    out.write(f"  {final_samples_tsv}\n\n")

    out.write("Ascertainment correction / fconst:\n")
    out.write("  Use nucleotide order: A C G T\n")
    out.write(f"  fconst = {fconst_string}\n")
    out.write(f"  fconst file: {fconst_txt}\n\n")

    out.write("Recommended BEAUti2 setup:\n")
    out.write("  1. Import beast2_subset_snps_variable_only.fasta.\n")
    out.write("  2. Set tip dates using decimal years from beast2_date_traits.tsv/csv.\n")
    out.write("  3. Use a nucleotide substitution model appropriate for SNP data.\n")
    out.write("  4. Enable ascertainment correction / constant site weights if using SNP-only alignment.\n")
    out.write("  5. Add constant site weights in A C G T order:\n")
    out.write(f"       {fconst_string}\n")
    out.write("  6. Consider starting with a strict clock, because TreeTime found a reasonable clock rate.\n")
    out.write("  7. If strict-clock ESS and posterior fit are acceptable, relaxed clock can be a sensitivity analysis.\n")
    out.write("  8. Run enough MCMC length to achieve ESS > 200 for key parameters.\n")
    out.write("  9. Export the strict-clock XML to:\n")
    out.write(f"       {xml_strict_path}\n")
    out.write(" 10. Export the relaxed-clock XML, if used, to:\n")
    out.write(f"       {xml_relaxed_path}\n")
    out.write(" 11. Run p4_04_run_beast2_and_summarize.sh after XML export.\n\n")

    out.write("TreeTime selected subset:\n")
    out.write(f"  TreeTime directory: {os.path.dirname(TREETIME_CLOCK)}\n")
    if clock["clock_rate"] is not None:
        out.write(f"  TreeTime clock rate: {clock['clock_rate']:.6e} sub/site/year\n")
    if clock["clock_r2"] is not None:
        out.write(f"  TreeTime root-to-tip R2: {clock['clock_r2']:.4f}\n")
    out.write("\n")


summary_lines = []
summary_lines.append("BEAST2 selected subset input preparation summary")
summary_lines.append("======================================")
summary_lines.append(f"Subsample file: {SUB_LIST}")
summary_lines.append(f"TreeTime dir: {os.path.dirname(TREETIME_CLOCK)}")
summary_lines.append(f"Original subsample rows: {len(sub)}")
summary_lines.append(f"Outliers removed: {outlier_count}")
summary_lines.append(f"Final samples: {len(final_samples)}")
summary_lines.append(f"SNP records written: {len(snp_variable_records)}")
summary_lines.append(f"SNP original columns after sample subset: {snp_filter_stats['original_snp_columns']}")
summary_lines.append(f"SNP variable columns kept: {snp_filter_stats['variable_snp_columns_kept']}")
summary_lines.append(f"SNP columns constant after subset: {snp_filter_stats['constant_snp_columns_after_subset']}")
summary_lines.append(f"Masked alignment records written: {len(masked_subset)}")
summary_lines.append(f"Masked alignment length: {masked_stats['masked_alignment_length']}")
summary_lines.append(f"GENOME_LENGTH: {GENOME_LENGTH}")
summary_lines.append(f"Target invariant total: {target_invariant_total}")
summary_lines.append(f"fconst A C G T: {fconst_string}")
summary_lines.append(f"Expected strict XML: {xml_strict_path}")
summary_lines.append(f"Expected relaxed XML: {xml_relaxed_path}")

if clock["clock_rate"] is not None:
    summary_lines.append(f"TreeTime clock rate: {clock['clock_rate']:.6e}")
if clock["clock_r2"] is not None:
    summary_lines.append(f"TreeTime root-to-tip R2: {clock['clock_r2']:.4f}")

summary_lines.append("")
summary_lines.append("Output files:")
for p in [
    final_samples_tsv,
    final_ids_txt,
    date_simple,
    date_traits,
    date_traits_csv,
    snp_all_out,
    snp_var_out,
    masked_out,
    fconst_tsv,
    fconst_txt,
    guide_file
]:
    summary_lines.append(f"  {p}")

with open(summary_file, "w") as out:
    out.write("\n".join(summary_lines) + "\n")

manifest_rows = [
    ("source_temporal_subset", SUB_LIST),
    ("treetime_dir", TREETIME_DIR),
    ("original_subset_samples", len(sub)),
    ("treetime_outliers_removed", outlier_count),
    ("final_beast2_samples", final_n),
    ("variable_snp_sites", variable_sites),
    ("fconst_A_C_G_T", fconst_string),
    ("strict_xml", xml_strict_path),
    ("relaxed_xml", xml_relaxed_path),
]
with open(manifest_file, "w") as out:
    out.write("field\tvalue\n")
    for k, v in manifest_rows:
        out.write(f"{k}\t{v}\n")

log(f"  [Output] {guide_file}")
log(f"  [Output] {summary_file}")
log(f"  [Output] {manifest_file}")

log("")
log("==========================================")
log("BEAST2 input preparation complete")
log("==========================================")
log(f"Final samples: {len(final_samples)}")
log(f"Variable SNP sites for BEAST2: {snp_var_len}")
log(f"fconst A C G T: {fconst_string}")
log("")
log("Key files:")
log(f"  Alignment: {snp_var_out}")
log(f"  Dates:     {date_traits}")
log(f"  fconst:    {fconst_txt}")
log(f"  Guide:     {guide_file}")
PYEOF

echo ""
echo "=========================================="
echo "Done."
echo "=========================================="
echo "Check summary:"
echo "  ${OUTDIR}/beast2_input_summary.txt"
echo ""
echo "BEAUti2 guide:"
echo "  ${OUTDIR}/beauti2_configuration_guide.txt"