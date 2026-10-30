#!/usr/bin/env bash
# p6_01_nichols_temporal_signal_analysis.sh
#
# Data analysis pipeline testing whether the Nichols lineage's lack of a
# usable molecular-clock temporal signal can be explained by, or rescued
# by: (1) reference-genome choice, (2) subsampling strategy, or (3) mixing
# with the validated SS14 stratified subset. Five stages, run in order:
#
#   Stage 1 (reference substitution): rebuild the full 1,861-genome dataset
#     under Nichols coordinates (mirroring the existing SS14-coordinate
#     pipeline) and re-test root-to-tip / date-randomization signal for
#     the full dataset and the Nichols-lineage subset.
#   Stage 2 (Nichols subsampling): test whether stratified subsampling of
#     the Nichols lineage (analogous to the strategy that rescued signal
#     for SS14) can rescue a usable temporal signal. Reuses the Stage 1
#     Nichols-reference tree; no realignment or rebuilding. Per-sample
#     root-to-tip points are saved for every profile.
#   Stage 3 (mixed-lineage subsampling): test whether combining the
#     validated SS14 stratified subset (n=555) with various Nichols
#     subsampling profiles improves or worsens overall temporal signal.
#     Uses SS14 coordinates throughout, pruning the existing SS14-reference
#     tree; no realignment or rebuilding. Per-sample root-to-tip points are
#     saved for every profile.
#   Stage 4 (cross-reference validation): two additional checks under
#     Nichols coordinates, reusing the Stage 1 Nichols-reference tree:
#     (a) the validated SS14 stratified subset (n=555) tested alone, to
#     check whether its temporal signal is robust to reference-genome
#     choice; (b) the full Stage 3 mixed-lineage matrix repeated under
#     Nichols coordinates, to check whether the signal collapse seen
#     under SS14 coordinates is reference-independent.
#   Stage 5: consolidated console summary across all stages (no new
#     computation; reads back the summary tables written above).
#
# Inputs:
#   refs/TPA_Nichols.fasta
#   consensus/nichols/*.nichols.snp_only.consensus.fasta            (Nichols-lineage short reads vs Nichols ref)
#   consensus/ss14_align2_nichols/*.ss14_a2n.snp_only.consensus.fasta (SS14-lineage short reads vs Nichols ref)
#   results/04_complete_genome_variant_calls/vcf/*_vs_nichols.pass.norm.vcf.gz (complete genomes vs Nichols ref)
#   results/03_iqtree/tpa_ss14_rooted.treefile                      (existing SS14-coordinate tree, for Stage 3)
#   results/11_temporal_signal/recommended_temporal_subset.tsv      (validated SS14 subset, 555 samples, for Stage 3/4)
#   TPA_master_metadata.xlsx
#   sample_dates.tsv
#   nichols_premask_regions.bed
#
# Outputs (all under results/p6_01_nichols_temporal_signal/):
#   01_reference_substitution/01_msa/TPA_nicholsref_msa.fasta
#   01_reference_substitution/02_masked/tpa_nichols_masked.fasta
#   01_reference_substitution/03_gubbins/tpa_nichols.filtered_polymorphic_sites.fasta
#   01_reference_substitution/04_iqtree/tpa_nichols_rooted.treefile
#   01_reference_substitution/05_temporal/temporal_signal_summary.tsv
#   01_reference_substitution/05_temporal/rtd_points_*.tsv
#   02_nichols_subsampling/profiles/<profile>_sample_ids.txt
#   02_nichols_subsampling/profiles/<profile>_tree.nwk
#   02_nichols_subsampling/profiles/rtd_points_<profile>.tsv          (now written for every profile)
#   02_nichols_subsampling/temporal_signal_summary.tsv
#   03_mixed_lineage/profiles/<profile>_sample_ids.txt
#   03_mixed_lineage/profiles/<profile>_tree.nwk
#   03_mixed_lineage/profiles/rtd_points_<profile>.tsv                (now written for every profile)
#   03_mixed_lineage/temporal_signal_summary.tsv
#   04_crossref_validation/profiles/<profile>_sample_ids.txt
#   04_crossref_validation/profiles/<profile>_tree.nwk
#   04_crossref_validation/profiles/rtd_points_<profile>.tsv
#   04_crossref_validation/temporal_signal_summary.tsv
#
# Usage:
#   bash script/p6_01_nichols_temporal_signal_analysis.sh

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$(pwd)}
REFERENCE=${REFERENCE:-"refs/TPA_Nichols.fasta"}
COMPLETE_GENOME_VCF_DIR=${COMPLETE_GENOME_VCF_DIR:-"results/04_complete_genome_variant_calls/vcf"}
NICHOLS_CONSENSUS_DIR=${NICHOLS_CONSENSUS_DIR:-"consensus/nichols"}
SS14_A2N_CONSENSUS_DIR=${SS14_A2N_CONSENSUS_DIR:-"consensus/ss14_align2_nichols"}
META_FILE=${META_FILE:-"TPA_master_metadata.xlsx"}
SAMPLE_DATES_FILE=${SAMPLE_DATES_FILE:-"sample_dates.tsv"}
PREMASK_BED=${PREMASK_BED:-"nichols_premask_regions.bed"}
SS14_TREE=${SS14_TREE:-"results/03_iqtree/tpa_ss14_rooted.treefile"}
SS14_SUBSET_FILE=${SS14_SUBSET_FILE:-"results/11_temporal_signal/recommended_temporal_subset.tsv"}
OUTROOT=${OUTROOT:-"results/p6_01_nichols_temporal_signal"}
OUTBASE="${OUTROOT}/01_reference_substitution"
OUTBASE_SUBSAMPLE="${OUTROOT}/02_nichols_subsampling"
OUTBASE_MIXED="${OUTROOT}/03_mixed_lineage"
OUTBASE_CROSSREF="${OUTROOT}/04_crossref_validation"
THREADS=${THREADS:-24}
MIN_COVERAGE=${MIN_COVERAGE:-0.75}
GUBBINS_ENV=${GUBBINS_ENV:-"gubbins"}
TP_WGS_ENV=${TP_WGS_ENV:-"tp_wgs"}
SEQ_LENGTH=${SEQ_LENGTH:-1083903}
CLOCK_FILTER=${CLOCK_FILTER:-3}
FORCE=${FORCE:-0}

cd "${PROJECT_ROOT}"

MSA_DIR="${OUTBASE}/01_msa"
MASK_DIR="${OUTBASE}/02_masked"
GUBBINS_DIR="${OUTBASE}/03_gubbins"
IQTREE_DIR="${OUTBASE}/04_iqtree"
TEMPORAL_DIR="${OUTBASE}/05_temporal"
TMP_DIR="${PROJECT_ROOT}/tmp_nichols_check"

MSA_FILE="${MSA_DIR}/TPA_nicholsref_msa.fasta"
MSA_SOURCE_TABLE="${MSA_DIR}/TPA_nicholsref_msa_sources.tsv"
MASKED_FASTA="${MASK_DIR}/tpa_nichols_masked.fasta"
QC_SUMMARY="${MASK_DIR}/qc_summary.tsv"
SNP_ALN="${GUBBINS_DIR}/tpa_nichols.filtered_polymorphic_sites.fasta"
UNROOTED_TREE="${IQTREE_DIR}/tpa_nichols.treefile"
ROOTED_TREE="${IQTREE_DIR}/tpa_nichols_rooted.treefile"
FCONST_FILE="${IQTREE_DIR}/fconst.txt"

mkdir -p "${MSA_DIR}" "${MASK_DIR}" "${GUBBINS_DIR}" "${IQTREE_DIR}" "${TEMPORAL_DIR}" "${TMP_DIR}" logs

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; exit 1; }

activate_env() {
    local env_name="$1"
    [[ -z "${env_name}" || "${env_name}" == "none" ]] && return 0
    command -v conda >/dev/null 2>&1 || die "conda not found"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

[[ -f "${REFERENCE}" ]] || die "Nichols reference not found: ${REFERENCE}"
[[ -f "${PREMASK_BED}" ]] || die "Nichols premask BED not found: ${PREMASK_BED}. Generate it first (liftover from SS14 coordinates)."
[[ -f "${META_FILE}" ]] || die "Metadata file not found: ${META_FILE}"
[[ -f "${SAMPLE_DATES_FILE}" ]] || die "Sample dates file not found: ${SAMPLE_DATES_FILE}"

# ============================================================
# Step 1: Build the Nichols-coordinate MSA from all four sources,
# mirroring the existing SS14-coordinate MSA construction.
# ============================================================
if [[ "${FORCE}" == "1" || ! -s "${MSA_FILE}" ]]; then
    log "Step 1: Building Nichols-coordinate MSA from all sample sources."

    export REFERENCE NICHOLS_CONSENSUS_DIR SS14_A2N_CONSENSUS_DIR COMPLETE_GENOME_VCF_DIR
    export MSA_FILE MSA_SOURCE_TABLE TMP_DIR

    python3 <<'PYEOF'
import os
import subprocess
from pathlib import Path
from Bio import SeqIO

reference = Path(os.environ["REFERENCE"])
nichols_dir = Path(os.environ["NICHOLS_CONSENSUS_DIR"])
ss14_a2n_dir = Path(os.environ["SS14_A2N_CONSENSUS_DIR"])
cg_vcf_dir = Path(os.environ["COMPLETE_GENOME_VCF_DIR"])
msa_file = Path(os.environ["MSA_FILE"])
source_table = Path(os.environ["MSA_SOURCE_TABLE"])
tmp_dir = Path(os.environ["TMP_DIR"])
tmp_dir.mkdir(parents=True, exist_ok=True)

ref_rec = next(SeqIO.parse(str(reference), "fasta"))
ref_seq = str(ref_rec.seq).upper()
ref_len = len(ref_seq)
print(f"[REF] {reference}: length={ref_len}")

rows = []
seen_ids = set()


def read_one_record(path):
    recs = list(SeqIO.parse(str(path), "fasta"))
    if len(recs) != 1:
        raise SystemExit(f"Expected 1 record in {path}, found {len(recs)}")
    return recs[0]


def add_record(sample_id, seq, source_type, src_path):
    status = "PASS"
    if sample_id in seen_ids:
        status = "SKIP_DUPLICATE_SAMPLE_ID"
    elif len(seq) != ref_len:
        status = "FAIL_LENGTH"
    else:
        seen_ids.add(sample_id)
    rows.append({
        "sample_id": sample_id, "source_type": source_type, "fasta": str(src_path),
        "length": len(seq), "expected_length": ref_len, "status": status, "sequence": seq,
    })


# Source 1: Nichols reference itself
add_record("NC_021490.2", ref_seq, "nichols_reference", str(reference))

# Source 2: Nichols-lineage short reads aligned directly to Nichols
n_src2 = 0
for f in sorted(nichols_dir.glob("*.nichols.snp_only.consensus.fasta")):
    sample_id = f.name.split(".")[0]
    rec = read_one_record(f)
    add_record(sample_id, str(rec.seq).upper(), "nichols_wgs", str(f))
    n_src2 += 1
print(f"[SRC] nichols (direct): {n_src2} files")

# Source 3: SS14-lineage short reads cross-aligned to Nichols
n_src3 = 0
for f in sorted(ss14_a2n_dir.glob("*.ss14_a2n.snp_only.consensus.fasta")):
    sample_id = f.name.split(".")[0]
    rec = read_one_record(f)
    add_record(sample_id, str(rec.seq).upper(), "ss14_align2_nichols", str(f))
    n_src3 += 1
print(f"[SRC] ss14_align2_nichols (cross): {n_src3} files")

# Source 4: complete genomes vs Nichols -- generate SNP-only consensus from
# the existing normalized VCFs (no realignment needed). Samples with zero
# variants (e.g. a genome identical or nearly identical to the Nichols
# reference) are expected and are handled by falling back to the reference
# sequence itself, since bcftools consensus on a record-free VCF should
# in principle return the reference unchanged, but is not relied upon to
# do so silently here.
n_src4 = 0
n_src4_failed = 0
failed_log = tmp_dir / "complete_genome_consensus_failures.tsv"
with failed_log.open("w") as flog:
    flog.write("sample_id\tvcf\treason\n")
    for vcf in sorted(cg_vcf_dir.glob("*_vs_nichols.pass.norm.vcf.gz")):
        sample_id = vcf.name.replace("_vs_nichols.pass.norm.vcf.gz", "")
        out_fa = tmp_dir / f"{sample_id}.cg_vs_nichols.consensus.fasta"

        if not vcf.exists() or vcf.stat().st_size == 0:
            flog.write(f"{sample_id}\t{vcf}\tmissing_or_empty_file\n")
            n_src4_failed += 1
            continue

        try:
            n_records = int(subprocess.run(
                f"bcftools view -H {vcf} | wc -l", shell=True, check=True,
                capture_output=True, text=True,
            ).stdout.strip())
        except Exception as e:
            flog.write(f"{sample_id}\t{vcf}\tcould_not_count_records:{e}\n")
            n_src4_failed += 1
            continue

        if n_records == 0:
            # No variants relative to Nichols reference: consensus equals
            # the reference sequence itself.
            add_record(sample_id, ref_seq, "complete_genome_vs_nichols_novar", str(vcf))
            n_src4 += 1
            continue

        try:
            subprocess.run(
                f"bcftools consensus -f {reference} {vcf} | sed 's/^>.*/>{sample_id}/' > {out_fa}",
                shell=True, check=True,
            )
            rec = read_one_record(out_fa)
            add_record(sample_id, str(rec.seq).upper(), "complete_genome_vs_nichols", str(vcf))
            n_src4 += 1
        except Exception as e:
            flog.write(f"{sample_id}\t{vcf}\tbcftools_consensus_failed:{e}\n")
            n_src4_failed += 1

print(f"[SRC] complete_genome_vs_nichols: {n_src4} ok, {n_src4_failed} failed/skipped (see {failed_log})")

failures = [r for r in rows if r["status"] == "FAIL_LENGTH"]
if failures:
    for r in failures:
        print(f"[LENGTH FAIL] {r['sample_id']}: {r['length']} (expected {ref_len}) from {r['fasta']}")

with msa_file.open("w") as out:
    for r in rows:
        if r["status"] != "PASS":
            continue
        out.write(f">{r['sample_id']}\n")
        seq = r["sequence"]
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

with source_table.open("w") as out:
    out.write("sample_id\tsource_type\tfasta\tlength\texpected_length\tstatus\n")
    for r in rows:
        out.write(f"{r['sample_id']}\t{r['source_type']}\t{r['fasta']}\t{r['length']}\t{r['expected_length']}\t{r['status']}\n")

n_pass = sum(1 for r in rows if r["status"] == "PASS")
print(f"[OUT] MSA records written: {n_pass}")
print(f"[OUT] MSA: {msa_file}")
print(f"[OUT] Source table: {source_table}")
PYEOF
else
    log "Step 1: Nichols-coordinate MSA exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${MSA_FILE}" ]] || die "MSA construction failed: ${MSA_FILE}"

# ============================================================
# Step 2: QC filter (>=75% coverage) and mask the lifted-over
# hypervariable regions (tpr paralogs, arp, etc.) in Nichols coordinates.
# ============================================================
if [[ "${FORCE}" == "1" || ! -s "${MASKED_FASTA}" ]]; then
    log "Step 2: QC filter and mask hypervariable regions (Nichols coordinates)."
    export MSA_FILE REFERENCE QC_SUMMARY MASKED_FASTA PREMASK_BED MIN_COVERAGE

    python3 <<'PYEOF'
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

msa_file = Path(os.environ["MSA_FILE"])
ref_file = Path(os.environ["REFERENCE"])
qc_summary = Path(os.environ["QC_SUMMARY"])
masked_fasta = Path(os.environ["MASKED_FASTA"])
premask_bed = Path(os.environ["PREMASK_BED"])
min_coverage = float(os.environ.get("MIN_COVERAGE", "0.75"))

regions = []
with premask_bed.open() as f:
    for line in f:
        chrom, start, end = line.strip().split("\t")
        regions.append((int(start), int(end)))  # BED 0-based half-open

ref_rec = next(SeqIO.parse(str(ref_file), "fasta"))
ref_seq = str(ref_rec.seq).upper()
ref_len = len(ref_seq) - ref_seq.count("N") - ref_seq.count("-")

rows = []
masked_records = []
skipped = 0
for rec in SeqIO.parse(str(msa_file), "fasta"):
    seq_str = str(rec.seq).upper()
    n_count = seq_str.count("N") + seq_str.count("-")
    valid_bp = len(seq_str) - n_count
    coverage = valid_bp / ref_len
    pass_qc = coverage >= min_coverage
    rows.append((rec.id, len(seq_str), valid_bp, n_count, round(coverage, 6), pass_qc))
    if not pass_qc:
        skipped += 1
        continue
    seq = list(seq_str)
    for start, end in regions:
        for i in range(start, min(end, len(seq))):
            seq[i] = "N"
    masked_records.append(SeqRecord(Seq("".join(seq)), id=rec.id, description=""))

with qc_summary.open("w") as handle:
    handle.write("sample\ttotal_bp\tvalid_bp\tN_count\tcoverage\tpass_QC\n")
    for row in rows:
        handle.write("\t".join(map(str, row)) + "\n")

if not masked_records:
    raise SystemExit("No sequences passed QC")

SeqIO.write(masked_records, str(masked_fasta), "fasta")
print(f"[QC] total={len(rows)} passed={len(masked_records)} skipped={skipped} min_coverage={min_coverage}")
print(f"[MASK] regions masked: {len(regions)}")
print(f"[OUT] {masked_fasta}")
PYEOF
else
    log "Step 2: masked alignment exists; skipping. Use FORCE=1 to rerun."
fi

# ============================================================
# Step 3: Gubbins recombination filtering (Nichols coordinates)
# ============================================================
if [[ "${FORCE}" == "1" || ! -s "${SNP_ALN}" ]]; then
    log "Step 3: Gubbins recombination filtering."
    (
        activate_env "${GUBBINS_ENV}"
        export TMPDIR="${TMP_DIR}"
        cd "${GUBBINS_DIR}"
        run_gubbins.py \
            --prefix tpa_nichols \
            --first-tree-builder rapidnj \
            --first-model JC \
            --tree-builder raxml \
            --model GTRGAMMA \
            --mar \
            --converge-method recombination \
            --iterations 10 \
            --threads "${THREADS}" \
            --remove-identical-sequences \
            --verbose \
            "${PROJECT_ROOT}/${MASKED_FASTA}" \
            2>&1 | tee gubbins.log
    )
else
    log "Step 3: Gubbins SNP alignment exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${SNP_ALN}" ]] || die "Missing Gubbins SNP alignment: ${SNP_ALN}"

# ============================================================
# Step 4: IQ-TREE ML phylogeny + midpoint root (Nichols coordinates)
# ============================================================
if [[ "${FORCE}" == "1" || ! -s "${FCONST_FILE}" ]]; then
    log "Step 4a: Constant-site counts for IQ-TREE (Nichols, masked)."
    export REFERENCE PREMASK_BED FCONST_FILE
    python3 <<'PYEOF'
import os
from pathlib import Path
from Bio import SeqIO

ref_file = Path(os.environ["REFERENCE"])
premask_bed = Path(os.environ["PREMASK_BED"])
fconst_file = Path(os.environ["FCONST_FILE"])

regions = []
with premask_bed.open() as f:
    for line in f:
        chrom, start, end = line.strip().split("\t")
        regions.append((int(start), int(end)))

rec = next(SeqIO.parse(str(ref_file), "fasta"))
seq = list(str(rec.seq).upper())
for start, end in regions:
    for i in range(start, min(end, len(seq))):
        seq[i] = "N"
seq = "".join(seq)
counts = {base: seq.count(base) for base in "ACGT"}
fconst = f"{counts['A']},{counts['C']},{counts['G']},{counts['T']}"
fconst_file.parent.mkdir(parents=True, exist_ok=True)
fconst_file.write_text(fconst + "\n")
print(f"[FCONST] {fconst}")
PYEOF
fi

if [[ "${FORCE}" == "1" || ! -s "${UNROOTED_TREE}" ]]; then
    log "Step 4b: IQ-TREE ML phylogeny (Nichols coordinates)."
    (
        activate_env "${TP_WGS_ENV}"
        FCONST=$(tr -d '[:space:]' < "${FCONST_FILE}")
        iqtree \
            -s "${SNP_ALN}" \
            -m GTR+G4 \
            -fconst "${FCONST}" \
            -bb 1000 \
            -nt "${THREADS}" \
            -pre "${IQTREE_DIR}/tpa_nichols" \
            -redo \
            > "${IQTREE_DIR}/iqtree.log" 2>&1
    )
else
    log "Step 4b: IQ-TREE tree exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${UNROOTED_TREE}" ]] || die "Missing IQ-TREE treefile: ${UNROOTED_TREE}"

if [[ "${FORCE}" == "1" || ! -s "${ROOTED_TREE}" ]]; then
    log "Step 4c: Midpoint-root the tree."
    export UNROOTED_TREE ROOTED_TREE
    (
        activate_env "${TP_WGS_ENV}"
        Rscript --vanilla - <<'REOF'
suppressPackageStartupMessages({library(ape); library(phytools)})
infile <- Sys.getenv("UNROOTED_TREE")
outfile <- Sys.getenv("ROOTED_TREE")
tree <- read.tree(infile)
rooted <- midpoint.root(tree)
write.tree(rooted, file = outfile)
cat(sprintf("[ROOT] wrote %s (%d tips)\n", outfile, length(rooted$tip.label)))
REOF
    )
else
    log "Step 4c: rooted tree exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${ROOTED_TREE}" ]] || die "Missing rooted tree: ${ROOTED_TREE}"

# ============================================================
# Step 5: Re-test temporal signal (full dataset + Nichols-lineage subset)
# under Nichols coordinates, reusing the verified get_rtd_year /
# date_randomization logic from p4_01_temporal_signal_screening.R.
# ============================================================
log "Step 5: Re-testing temporal signal under Nichols coordinates."
export ROOTED_TREE META_FILE SAMPLE_DATES_FILE TEMPORAL_DIR
(
    activate_env "${TP_WGS_ENV}"
    Rscript --vanilla - <<'REOF'
suppressPackageStartupMessages({library(ape); library(phytools); library(readxl)})

MIN_N_TEMPORAL_TEST <- 10
info <- function(...) cat(sprintf(...), "\n")

tree_file <- Sys.getenv("ROOTED_TREE")
meta_file <- Sys.getenv("META_FILE")
dates_file <- Sys.getenv("SAMPLE_DATES_FILE")
outdir <- Sys.getenv("TEMPORAL_DIR")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(tree_file)

# Sampling years: headerless two-column TSV (sample_id, year), confirmed
# against sample_dates.tsv -- not the metadata Excel file, which does not
# reliably expose this column under a known name.
dates_tbl <- read.table(dates_file, sep = "\t", header = FALSE,
                         col.names = c("sample_id", "year"),
                         stringsAsFactors = FALSE)
dates <- setNames(suppressWarnings(as.numeric(dates_tbl$year)), dates_tbl$sample_id)

# Lineage assignment still comes from the metadata Excel file.
meta <- read_excel(meta_file)
sample_col <- "SRA_ID (Key)"
lineage_col <- "Clade (Lineage)"
lineages <- setNames(meta[[lineage_col]], meta[[sample_col]])

get_rtd_year <- function(tree, dates) {
    tips_with_dates <- intersect(tree$tip.label, names(dates))
    tips_with_dates <- tips_with_dates[!is.na(dates[tips_with_dates])]
    if (length(tips_with_dates) < MIN_N_TEMPORAL_TEST) return(NULL)
    if (length(tree$tip.label) > length(tips_with_dates)) {
        drop <- setdiff(tree$tip.label, tips_with_dates)
        tree <- drop.tip(tree, drop)
    }
    rtd_raw <- node.depth.edgelength(tree)
    rtd_vec <- rtd_raw[seq_along(tree$tip.label)]
    names(rtd_vec) <- tree$tip.label
    yr_vec <- as.numeric(dates[tree$tip.label])
    ok <- is.finite(rtd_vec) & is.finite(yr_vec)
    if (sum(ok) < MIN_N_TEMPORAL_TEST) return(NULL)
    list(rtd = rtd_vec[ok], year = yr_vec[ok], tips = tree$tip.label[ok])
}

calc_regression <- function(rtd, year) {
    fit <- lm(rtd ~ year)
    sm <- summary(fit)
    list(r2 = sm$r.squared, slope = coef(fit)[2])
}

date_randomization <- function(rtd, year, n_perm = 1000, seed = 42) {
    set.seed(seed)
    real <- calc_regression(rtd, year)
    perm_r2 <- numeric(n_perm)
    for (i in seq_len(n_perm)) perm_r2[i] <- calc_regression(rtd, sample(year))$r2
    list(n_samples = length(rtd), real_r2 = real$r2, slope = real$slope,
         perm_95pct = quantile(perm_r2, 0.95), p_value = mean(perm_r2 >= real$r2))
}

run_one <- function(dataset, label, tree, dates) {
    info("Temporal test: %s", label)
    data <- get_rtd_year(tree, dates)
    if (is.null(data)) { info("  Insufficient dated tips"); return(NULL) }

    rtd_points <- data.frame(sample = data$tips, year = data$year, rtd = data$rtd)
    write.table(rtd_points, file.path(outdir, sprintf("rtd_points_%s.tsv", dataset)),
                sep = "\t", row.names = FALSE, quote = FALSE)

    res <- date_randomization(data$rtd, data$year)
    info("  n=%d r2=%.4f slope=%.3e p=%.4f", res$n_samples, res$real_r2, res$slope, res$p_value)
    data.frame(dataset = dataset, label = label, n_samples = res$n_samples,
               real_r2 = res$real_r2, slope = res$slope,
               perm_95pct = res$perm_95pct, p_value = res$p_value)
}

results <- list()

# Full dataset under Nichols coordinates
results[[length(results) + 1]] <- run_one("full_nicholsref", "Full dataset (Nichols reference)", tree, dates)

# Nichols-lineage subset under Nichols coordinates
nichols_ids <- names(lineages)[lineages == "Nichols" & !is.na(lineages)]
nichols_tips <- intersect(tree$tip.label, nichols_ids)
if (length(nichols_tips) >= MIN_N_TEMPORAL_TEST) {
    tree_nichols_sub <- keep.tip(tree, nichols_tips)
    results[[length(results) + 1]] <- run_one(
        "nichols_lineage_nicholsref", "Nichols-lineage subset (Nichols reference)",
        tree_nichols_sub, dates
    )
}

summary_df <- do.call(rbind, results[!sapply(results, is.null)])
write.table(summary_df, file.path(outdir, "temporal_signal_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
info("\n[OUT] %s", file.path(outdir, "temporal_signal_summary.tsv"))
print(summary_df)
REOF
)

log "Stage 1 done: ${OUTBASE}/05_temporal/temporal_signal_summary.tsv"

# ============================================================
# Stage 2: Nichols-lineage subsampling strategies (from p4_98)
# Reuses the Stage 1 Nichols-reference tree; no realignment or rebuilding.
# ============================================================
log "Stage 2: Testing whether stratified subsampling rescues Nichols temporal signal."

PROFILE_DIR_SUBSAMPLE="${OUTBASE_SUBSAMPLE}/profiles"
mkdir -p "${PROFILE_DIR_SUBSAMPLE}"

export ROOTED_TREE META_FILE SAMPLE_DATES_FILE
export PROFILE_DIR_SUBSAMPLE OUTBASE_SUBSAMPLE
(
    activate_env "${TP_WGS_ENV}"
    Rscript --vanilla - <<'REOF'
suppressPackageStartupMessages({library(ape); library(readxl)})

MIN_N_TEMPORAL_TEST <- 10
info <- function(...) cat(sprintf(...), "\n")
set.seed(42)

tree_file   <- Sys.getenv("ROOTED_TREE")
dates_file  <- Sys.getenv("SAMPLE_DATES_FILE")
meta_file   <- Sys.getenv("META_FILE")
profile_dir <- Sys.getenv("PROFILE_DIR_SUBSAMPLE")
outdir      <- Sys.getenv("OUTBASE_SUBSAMPLE")

tree <- read.tree(tree_file)

dates_tbl <- read.table(dates_file, sep = "\t", header = FALSE,
                         col.names = c("sample_id", "year"), stringsAsFactors = FALSE)
dates <- setNames(suppressWarnings(as.numeric(dates_tbl$year)), dates_tbl$sample_id)

meta <- read_excel(meta_file)
sample_col  <- "SRA_ID (Key)"
lineage_col <- "Clade (Lineage)"
lineages <- setNames(meta[[lineage_col]], meta[[sample_col]])

nichols_ids <- names(lineages)[lineages == "Nichols" & !is.na(lineages)]
nichols_ids <- intersect(nichols_ids, tree$tip.label)
nichols_ids <- nichols_ids[!is.na(dates[nichols_ids])]
info("Nichols-lineage samples with tree placement + date: %d", length(nichols_ids))

nichols_years <- dates[nichols_ids]

get_rtd_year <- function(tree, dates) {
    tips_with_dates <- intersect(tree$tip.label, names(dates))
    tips_with_dates <- tips_with_dates[!is.na(dates[tips_with_dates])]
    if (length(tips_with_dates) < MIN_N_TEMPORAL_TEST) return(NULL)
    if (length(tree$tip.label) > length(tips_with_dates)) {
        tree <- keep.tip(tree, tips_with_dates)
    }
    rtd_raw <- node.depth.edgelength(tree)
    rtd_vec <- rtd_raw[seq_along(tree$tip.label)]
    names(rtd_vec) <- tree$tip.label
    yr_vec <- as.numeric(dates[tree$tip.label])
    ok <- is.finite(rtd_vec) & is.finite(yr_vec)
    if (sum(ok) < MIN_N_TEMPORAL_TEST) return(NULL)
    list(rtd = rtd_vec[ok], year = yr_vec[ok], tips = tree$tip.label[ok])
}

calc_regression <- function(rtd, year) {
    fit <- lm(rtd ~ year)
    sm <- summary(fit)
    list(r2 = sm$r.squared, slope = coef(fit)[2])
}

date_randomization <- function(rtd, year, n_perm = 1000, seed = 42) {
    set.seed(seed)
    real <- calc_regression(rtd, year)
    perm_r2 <- numeric(n_perm)
    for (i in seq_len(n_perm)) perm_r2[i] <- calc_regression(rtd, sample(year))$r2
    list(n_samples = length(rtd), real_r2 = real$r2, slope = real$slope,
         perm_95pct = quantile(perm_r2, 0.95), p_value = mean(perm_r2 >= real$r2))
}

run_one <- function(profile_name, label, tip_ids) {
    info("--- Profile: %s (%s) ---", profile_name, label)
    if (length(tip_ids) < MIN_N_TEMPORAL_TEST) {
        info("  Too few samples (n=%d), skipping.", length(tip_ids)); return(NULL)
    }
    sub_tree <- keep.tip(tree, tip_ids)
    writeLines(tip_ids, file.path(profile_dir, sprintf("%s_sample_ids.txt", profile_name)))
    write.tree(sub_tree, file.path(profile_dir, sprintf("%s_tree.nwk", profile_name)))

    data <- get_rtd_year(sub_tree, dates)
    if (is.null(data)) { info("  Insufficient dated tips after pruning."); return(NULL) }

    rtd_points <- data.frame(sample = data$tips, year = data$year, rtd = data$rtd)
    write.table(rtd_points, file.path(profile_dir, sprintf("rtd_points_%s.tsv", profile_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)

    res <- date_randomization(data$rtd, data$year)
    info("  n=%d  r2=%.4f  slope=%.3e  perm_95pct=%.4f  p=%.4f",
         res$n_samples, res$real_r2, res$slope, res$perm_95pct, res$p_value)

    data.frame(profile = profile_name, label = label, n_samples = res$n_samples,
               real_r2 = res$real_r2, slope = res$slope,
               perm_95pct = res$perm_95pct, p_value = res$p_value)
}

stratified_subsample <- function(ids, years, year_breaks, per_bin_cap, early_cutoff = 1980) {
    df <- data.frame(id = ids, year = years, stringsAsFactors = FALSE)
    early <- df[df$year < early_cutoff, ]
    late  <- df[df$year >= early_cutoff, ]

    late$bin <- cut(late$year, breaks = year_breaks, include.lowest = TRUE)
    kept_late <- do.call(rbind, lapply(split(late, late$bin), function(g) {
        if (nrow(g) == 0) return(g)
        if (nrow(g) <= per_bin_cap) return(g)
        g[sample(seq_len(nrow(g)), per_bin_cap), ]
    }))

    c(early$id, kept_late$id)
}

results <- list()
results[[length(results) + 1]] <- run_one(
    "profile_0_unsampled", "Full Nichols-lineage subset (no subsampling)", nichols_ids
)

ids_A <- stratified_subsample(
    nichols_ids, nichols_years,
    year_breaks = seq(1980, 2030, by = 1), per_bin_cap = 15
)
results[[length(results) + 1]] <- run_one(
    "profile_A_light_downsample", "Light downsampling (cap 15/year post-1980)", ids_A
)

ids_B <- stratified_subsample(
    nichols_ids, nichols_years,
    year_breaks = seq(1980, 2030, by = 1), per_bin_cap = 5
)
results[[length(results) + 1]] <- run_one(
    "profile_B_aggressive_downsample", "Aggressive downsampling (cap 5/year post-1980)", ids_B
)

ids_C <- stratified_subsample(
    nichols_ids, nichols_years,
    year_breaks = c(1980, 2000, 2010, 2020, 2030), per_bin_cap = 30
)
results[[length(results) + 1]] <- run_one(
    "profile_C_era_stratified", "Era-stratified (cap 30/era: 1980-2000/2000-10/2010-20/2020+)", ids_C
)

summary_df <- do.call(rbind, results[!sapply(results, is.null)])
out_path <- file.path(outdir, "temporal_signal_summary.tsv")
write.table(summary_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
info("\n[OUT] %s", out_path)
print(summary_df)
REOF
)
log "Stage 2 done: ${OUTBASE_SUBSAMPLE}/temporal_signal_summary.tsv"

# ============================================================
# Stage 3: Mixed SS14 + Nichols subsampling (from p4_97)
# Uses SS14 coordinates throughout, pruning the existing SS14-reference
# tree -- no realignment or rebuilding.
# ============================================================
log "Stage 3: Testing whether mixing the validated SS14 subset with Nichols genomes helps or hurts."

[[ -f "${SS14_TREE}" ]] || die "SS14 tree not found: ${SS14_TREE}"
[[ -f "${SS14_SUBSET_FILE}" ]] || die "${SS14_SUBSET_FILE} not found"

PROFILE_DIR_MIXED="${OUTBASE_MIXED}/profiles"
mkdir -p "${PROFILE_DIR_MIXED}"

export SS14_TREE SS14_SUBSET_FILE META_FILE SAMPLE_DATES_FILE
export PROFILE_DIR_MIXED OUTBASE_MIXED
(
    activate_env "${TP_WGS_ENV}"
    Rscript --vanilla - <<'REOF'
suppressPackageStartupMessages({library(ape); library(readxl)})

MIN_N_TEMPORAL_TEST <- 10
info <- function(...) cat(sprintf(...), "\n")
set.seed(42)

ss14_tree_file   <- Sys.getenv("SS14_TREE")
ss14_subset_file <- Sys.getenv("SS14_SUBSET_FILE")
dates_file       <- Sys.getenv("SAMPLE_DATES_FILE")
meta_file        <- Sys.getenv("META_FILE")
profile_dir      <- Sys.getenv("PROFILE_DIR_MIXED")
outdir           <- Sys.getenv("OUTBASE_MIXED")

full_tree <- read.tree(ss14_tree_file)

dates_tbl <- read.table(dates_file, sep = "\t", header = FALSE,
                         col.names = c("sample_id", "year"), stringsAsFactors = FALSE)
dates <- setNames(suppressWarnings(as.numeric(dates_tbl$year)), dates_tbl$sample_id)

meta <- read_excel(meta_file)
sample_col  <- "SRA_ID (Key)"
lineage_col <- "Clade (Lineage)"
lineages <- setNames(meta[[lineage_col]], meta[[sample_col]])

ss14_subset <- read.table(ss14_subset_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ss14_ids <- intersect(ss14_subset$sample, full_tree$tip.label)
info("Validated SS14 stratified subset: %d samples (of %d in file)", length(ss14_ids), nrow(ss14_subset))

nichols_ids_all <- names(lineages)[lineages == "Nichols" & !is.na(lineages)]
nichols_ids_all <- intersect(nichols_ids_all, full_tree$tip.label)
nichols_ids_all <- nichols_ids_all[!is.na(dates[nichols_ids_all])]
nichols_years_all <- dates[nichols_ids_all]
info("Nichols-lineage samples available (SS14-coordinate tree, with date): %d", length(nichols_ids_all))

get_rtd_year <- function(tree, dates) {
    tips_with_dates <- intersect(tree$tip.label, names(dates))
    tips_with_dates <- tips_with_dates[!is.na(dates[tips_with_dates])]
    if (length(tips_with_dates) < MIN_N_TEMPORAL_TEST) return(NULL)
    if (length(tree$tip.label) > length(tips_with_dates)) {
        tree <- keep.tip(tree, tips_with_dates)
    }
    rtd_raw <- node.depth.edgelength(tree)
    rtd_vec <- rtd_raw[seq_along(tree$tip.label)]
    names(rtd_vec) <- tree$tip.label
    yr_vec <- as.numeric(dates[tree$tip.label])
    ok <- is.finite(rtd_vec) & is.finite(yr_vec)
    if (sum(ok) < MIN_N_TEMPORAL_TEST) return(NULL)
    list(rtd = rtd_vec[ok], year = yr_vec[ok], tips = tree$tip.label[ok])
}

calc_regression <- function(rtd, year) {
    fit <- lm(rtd ~ year)
    sm <- summary(fit)
    list(r2 = sm$r.squared, slope = coef(fit)[2])
}

date_randomization <- function(rtd, year, n_perm = 1000, seed = 42) {
    set.seed(seed)
    real <- calc_regression(rtd, year)
    perm_r2 <- numeric(n_perm)
    for (i in seq_len(n_perm)) perm_r2[i] <- calc_regression(rtd, sample(year))$r2
    list(n_samples = length(rtd), real_r2 = real$r2, slope = real$slope,
         perm_95pct = quantile(perm_r2, 0.95), p_value = mean(perm_r2 >= real$r2))
}

run_one <- function(profile_name, label, tip_ids) {
    info("--- Profile: %s (%s) ---", profile_name, label)
    tip_ids <- intersect(tip_ids, full_tree$tip.label)
    if (length(tip_ids) < MIN_N_TEMPORAL_TEST) {
        info("  Too few samples (n=%d), skipping.", length(tip_ids)); return(NULL)
    }
    sub_tree <- keep.tip(full_tree, tip_ids)
    writeLines(tip_ids, file.path(profile_dir, sprintf("%s_sample_ids.txt", profile_name)))
    write.tree(sub_tree, file.path(profile_dir, sprintf("%s_tree.nwk", profile_name)))

    data <- get_rtd_year(sub_tree, dates)
    if (is.null(data)) { info("  Insufficient dated tips after pruning."); return(NULL) }

    rtd_points <- data.frame(sample = data$tips, year = data$year, rtd = data$rtd)
    write.table(rtd_points, file.path(profile_dir, sprintf("rtd_points_%s.tsv", profile_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)

    res <- date_randomization(data$rtd, data$year)
    info("  n=%d  r2=%.4f  slope=%.3e  perm_95pct=%.4f  p=%.4f",
         res$n_samples, res$real_r2, res$slope, res$perm_95pct, res$p_value)

    data.frame(profile = profile_name, label = label, n_samples = res$n_samples,
               real_r2 = res$real_r2, slope = res$slope,
               perm_95pct = res$perm_95pct, p_value = res$p_value)
}

stratified_nichols_subsample <- function(ids, years, per_bin_cap, early_cutoff = 1980) {
    df <- data.frame(id = ids, year = years, stringsAsFactors = FALSE)
    early <- df[df$year < early_cutoff, ]
    late  <- df[df$year >= early_cutoff, ]
    late$bin <- cut(late$year, breaks = seq(early_cutoff, 2030, by = 1), include.lowest = TRUE)
    kept_late <- do.call(rbind, lapply(split(late, late$bin), function(g) {
        if (nrow(g) == 0) return(g)
        if (nrow(g) <= per_bin_cap) return(g)
        g[sample(seq_len(nrow(g)), per_bin_cap), ]
    }))
    c(early$id, kept_late$id)
}

results <- list()

results[[length(results) + 1]] <- run_one(
    "ref1_ss14_alone", "SS14 stratified subset alone (n=555, validated)", ss14_ids
)

nichols_light <- stratified_nichols_subsample(nichols_ids_all, nichols_years_all, per_bin_cap = 15)
results[[length(results) + 1]] <- run_one(
    "ref2_nichols_alone_light", "Nichols alone, light downsample (cap 15/yr)", nichols_light
)

nichols_light_2 <- stratified_nichols_subsample(nichols_ids_all, nichols_years_all, per_bin_cap = 15)
results[[length(results) + 1]] <- run_one(
    "mixed_A_ss14_plus_nichols_light",
    "SS14 (555) + Nichols light downsample (cap 15/yr)",
    union(ss14_ids, nichols_light_2)
)

nichols_aggressive <- stratified_nichols_subsample(nichols_ids_all, nichols_years_all, per_bin_cap = 5)
results[[length(results) + 1]] <- run_one(
    "mixed_B_ss14_plus_nichols_aggressive",
    "SS14 (555) + Nichols aggressive downsample (cap 5/yr)",
    union(ss14_ids, nichols_aggressive)
)

nichols_early_only <- nichols_ids_all[nichols_years_all < 1980]
results[[length(results) + 1]] <- run_one(
    "mixed_C_ss14_plus_nichols_early_anchors_only",
    "SS14 (555) + Nichols pre-1980 samples only (n=13)",
    union(ss14_ids, nichols_early_only)
)

summary_df <- do.call(rbind, results[!sapply(results, is.null)])
out_path <- file.path(outdir, "temporal_signal_summary.tsv")
write.table(summary_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
info("\n[OUT] %s", out_path)
print(summary_df)
REOF
)
log "Stage 3 done: ${OUTBASE_MIXED}/temporal_signal_summary.tsv"

# ============================================================
# Stage 4: Backfill per-sample root-to-tip detail tables for the
# profiles needed by the visualization script (p6_02) but not already
# saved with per-sample distances (from p4_95). Root-to-tip distance is
# read directly from each tree's existing branch lengths -- no
# realignment or tree-rebuilding is performed or needed.
# ============================================================
log "Stage 4: Cross-reference validation -- SS14 subset and mixed matrix under Nichols coordinates."

[[ -f "${SS14_SUBSET_FILE}" ]] || die "${SS14_SUBSET_FILE} not found"

PROFILE_DIR_CROSSREF="${OUTBASE_CROSSREF}/profiles"
mkdir -p "${PROFILE_DIR_CROSSREF}"

export ROOTED_TREE META_FILE SAMPLE_DATES_FILE SS14_SUBSET_FILE
export PROFILE_DIR_CROSSREF OUTBASE_CROSSREF
(
    activate_env "${TP_WGS_ENV}"
    Rscript --vanilla - <<'REOF'
suppressPackageStartupMessages({library(ape); library(readxl)})

MIN_N_TEMPORAL_TEST <- 10
info <- function(...) cat(sprintf(...), "\n")
set.seed(42)

tree_file        <- Sys.getenv("ROOTED_TREE")
dates_file       <- Sys.getenv("SAMPLE_DATES_FILE")
meta_file        <- Sys.getenv("META_FILE")
ss14_subset_file <- Sys.getenv("SS14_SUBSET_FILE")
profile_dir      <- Sys.getenv("PROFILE_DIR_CROSSREF")
outdir           <- Sys.getenv("OUTBASE_CROSSREF")

# Nichols-reference tree built in Stage 1 -- reused here, not rebuilt.
full_tree <- read.tree(tree_file)

dates_tbl <- read.table(dates_file, sep = "\t", header = FALSE,
                         col.names = c("sample_id", "year"), stringsAsFactors = FALSE)
dates <- setNames(suppressWarnings(as.numeric(dates_tbl$year)), dates_tbl$sample_id)

meta <- read_excel(meta_file)
sample_col  <- "SRA_ID (Key)"
lineage_col <- "Clade (Lineage)"
lineages <- setNames(meta[[lineage_col]], meta[[sample_col]])

ss14_subset <- read.table(ss14_subset_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ss14_ids <- intersect(ss14_subset$sample, full_tree$tip.label)
info("Validated SS14 stratified subset placed on Nichols-reference tree: %d samples (of %d in file)",
     length(ss14_ids), nrow(ss14_subset))

nichols_ids_all <- names(lineages)[lineages == "Nichols" & !is.na(lineages)]
nichols_ids_all <- intersect(nichols_ids_all, full_tree$tip.label)
nichols_ids_all <- nichols_ids_all[!is.na(dates[nichols_ids_all])]
nichols_years_all <- dates[nichols_ids_all]
info("Nichols-lineage samples available (Nichols-reference tree, with date): %d", length(nichols_ids_all))

get_rtd_year <- function(tree, dates) {
    tips_with_dates <- intersect(tree$tip.label, names(dates))
    tips_with_dates <- tips_with_dates[!is.na(dates[tips_with_dates])]
    if (length(tips_with_dates) < MIN_N_TEMPORAL_TEST) return(NULL)
    if (length(tree$tip.label) > length(tips_with_dates)) {
        tree <- keep.tip(tree, tips_with_dates)
    }
    rtd_raw <- node.depth.edgelength(tree)
    rtd_vec <- rtd_raw[seq_along(tree$tip.label)]
    names(rtd_vec) <- tree$tip.label
    yr_vec <- as.numeric(dates[tree$tip.label])
    ok <- is.finite(rtd_vec) & is.finite(yr_vec)
    if (sum(ok) < MIN_N_TEMPORAL_TEST) return(NULL)
    list(rtd = rtd_vec[ok], year = yr_vec[ok], tips = tree$tip.label[ok])
}

calc_regression <- function(rtd, year) {
    fit <- lm(rtd ~ year)
    sm <- summary(fit)
    list(r2 = sm$r.squared, slope = coef(fit)[2])
}

date_randomization <- function(rtd, year, n_perm = 1000, seed = 42) {
    set.seed(seed)
    real <- calc_regression(rtd, year)
    perm_r2 <- numeric(n_perm)
    for (i in seq_len(n_perm)) perm_r2[i] <- calc_regression(rtd, sample(year))$r2
    list(n_samples = length(rtd), real_r2 = real$r2, slope = real$slope,
         perm_95pct = quantile(perm_r2, 0.95), p_value = mean(perm_r2 >= real$r2))
}

run_one <- function(profile_name, label, tip_ids) {
    info("--- Profile: %s (%s) ---", profile_name, label)
    tip_ids <- intersect(tip_ids, full_tree$tip.label)
    if (length(tip_ids) < MIN_N_TEMPORAL_TEST) {
        info("  Too few samples (n=%d), skipping.", length(tip_ids)); return(NULL)
    }
    sub_tree <- keep.tip(full_tree, tip_ids)
    writeLines(tip_ids, file.path(profile_dir, sprintf("%s_sample_ids.txt", profile_name)))
    write.tree(sub_tree, file.path(profile_dir, sprintf("%s_tree.nwk", profile_name)))

    data <- get_rtd_year(sub_tree, dates)
    if (is.null(data)) { info("  Insufficient dated tips after pruning."); return(NULL) }

    rtd_points <- data.frame(sample = data$tips, year = data$year, rtd = data$rtd)
    write.table(rtd_points, file.path(profile_dir, sprintf("rtd_points_%s.tsv", profile_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)

    res <- date_randomization(data$rtd, data$year)
    info("  n=%d  r2=%.4f  slope=%.3e  perm_95pct=%.4f  p=%.4f",
         res$n_samples, res$real_r2, res$slope, res$perm_95pct, res$p_value)

    data.frame(profile = profile_name, label = label, n_samples = res$n_samples,
               real_r2 = res$real_r2, slope = res$slope,
               perm_95pct = res$perm_95pct, p_value = res$p_value)
}

stratified_nichols_subsample <- function(ids, years, per_bin_cap, early_cutoff = 1980) {
    df <- data.frame(id = ids, year = years, stringsAsFactors = FALSE)
    early <- df[df$year < early_cutoff, ]
    late  <- df[df$year >= early_cutoff, ]
    late$bin <- cut(late$year, breaks = seq(early_cutoff, 2030, by = 1), include.lowest = TRUE)
    kept_late <- do.call(rbind, lapply(split(late, late$bin), function(g) {
        if (nrow(g) == 0) return(g)
        if (nrow(g) <= per_bin_cap) return(g)
        g[sample(seq_len(nrow(g)), per_bin_cap), ]
    }))
    c(early$id, kept_late$id)
}

results <- list()

# (a) Validated SS14 stratified subset, tested alone, under Nichols
# reference coordinates -- checks whether its temporal signal survives
# a reference-genome substitution on its own, independent of any
# mixing with Nichols genomes.
results[[length(results) + 1]] <- run_one(
    "crossref_a_ss14_alone_nicholsref",
    "SS14 stratified subset alone (n=555), Nichols reference",
    ss14_ids
)

# (b) Full Stage 3 mixed-lineage matrix, repeated under Nichols
# reference coordinates.
results[[length(results) + 1]] <- run_one(
    "crossref_b_nichols_alone_light_nicholsref",
    "Nichols alone, light downsample (cap 15/yr), Nichols reference",
    stratified_nichols_subsample(nichols_ids_all, nichols_years_all, per_bin_cap = 15)
)

results[[length(results) + 1]] <- run_one(
    "crossref_c_mixed_ss14_plus_nichols_light_nicholsref",
    "SS14 (555) + Nichols light downsample (cap 15/yr), Nichols reference",
    union(ss14_ids, stratified_nichols_subsample(nichols_ids_all, nichols_years_all, per_bin_cap = 15))
)

results[[length(results) + 1]] <- run_one(
    "crossref_d_mixed_ss14_plus_nichols_aggressive_nicholsref",
    "SS14 (555) + Nichols aggressive downsample (cap 5/yr), Nichols reference",
    union(ss14_ids, stratified_nichols_subsample(nichols_ids_all, nichols_years_all, per_bin_cap = 5))
)

nichols_early_only <- nichols_ids_all[nichols_years_all < 1980]
results[[length(results) + 1]] <- run_one(
    "crossref_e_mixed_ss14_plus_nichols_early_anchors_only_nicholsref",
    "SS14 (555) + Nichols pre-1980 samples only (n=13), Nichols reference",
    union(ss14_ids, nichols_early_only)
)

summary_df <- do.call(rbind, results[!sapply(results, is.null)])
out_path <- file.path(outdir, "temporal_signal_summary.tsv")
write.table(summary_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
info("\n[OUT] %s", out_path)
print(summary_df)
REOF
)
log "Stage 4 done: ${OUTBASE_CROSSREF}/temporal_signal_summary.tsv"

# ============================================================
# Stage 5: Consolidated console summary across all four stages.
# No new computation -- reads back the summary tables already written.
# ============================================================
log "Stage 5: Consolidated summary."
export OUTBASE OUTBASE_SUBSAMPLE OUTBASE_MIXED OUTBASE_CROSSREF
(
    activate_env "${TP_WGS_ENV}"
    Rscript --vanilla - <<'REOF'
info <- function(...) cat(sprintf(...), "\n")

read_summary <- function(path) {
    if (!file.exists(path)) return(NULL)
    read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

tables <- list(
    "Stage 1: reference substitution"         = read_summary(file.path(Sys.getenv("OUTBASE"), "05_temporal", "temporal_signal_summary.tsv")),
    "Stage 2: Nichols subsampling"             = read_summary(file.path(Sys.getenv("OUTBASE_SUBSAMPLE"), "temporal_signal_summary.tsv")),
    "Stage 3: mixed lineage (SS14 reference)"  = read_summary(file.path(Sys.getenv("OUTBASE_MIXED"), "temporal_signal_summary.tsv")),
    "Stage 4: cross-reference validation"      = read_summary(file.path(Sys.getenv("OUTBASE_CROSSREF"), "temporal_signal_summary.tsv"))
)

for (name in names(tables)) {
    info("\n=== %s ===", name)
    tbl <- tables[[name]]
    if (is.null(tbl)) { info("  (summary table not found)"); next }
    print(tbl)
}
REOF
)
log ""
log "Temporal-coordinate analyses complete. Generate current publication figures through the main figure workflow."