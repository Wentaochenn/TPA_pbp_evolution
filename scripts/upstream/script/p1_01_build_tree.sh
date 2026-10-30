#!/usr/bin/env bash
# p1_01_build_tree.sh
#
# Build the recombination-filtered maximum-likelihood phylogeny.
#
# Inputs:
#   msa.fasta or consensus/TPA_ss14ref_msa.fasta
#   refs/TPA_SS14.fasta
#
# Outputs:
#   results/01_masked_alignment/qc_summary.tsv
#   results/01_masked_alignment/tpa_ss14_masked.fasta
#   results/02_gubbins/tpa_ss14.filtered_polymorphic_sites.fasta
#   results/03_iqtree/tpa_ss14.treefile
#   results/03_iqtree/tpa_ss14_rooted.treefile
#
# Usage:
#   bash script/p1_01_build_tree.sh

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$(pwd)}
INPUT_MSA=${INPUT_MSA:-}
REFERENCE=${REFERENCE:-"refs/TPA_SS14.fasta"}
OUTDIR=${OUTDIR:-"results"}
THREADS=${THREADS:-24}
MIN_COVERAGE=${MIN_COVERAGE:-0.75}
GUBBINS_ENV=${GUBBINS_ENV:-"gubbins"}
TP_WGS_ENV=${TP_WGS_ENV:-"tp_wgs"}
FORCE=${FORCE:-0}

cd "${PROJECT_ROOT}"

MASK_DIR="${OUTDIR}/01_masked_alignment"
GUBBINS_DIR="${OUTDIR}/02_gubbins"
IQTREE_DIR="${OUTDIR}/03_iqtree"
TMP_DIR="${PROJECT_ROOT}/tmp"

QC_SUMMARY="${MASK_DIR}/qc_summary.tsv"
PASSED_IDS="${MASK_DIR}/passed_ids.txt"
FAILED_IDS="${MASK_DIR}/failed_ids.txt"
MASKED_FASTA="${MASK_DIR}/tpa_ss14_masked.fasta"
PREMASK_BED="${MASK_DIR}/premask_regions.bed"
PREMASK_REPORT="${MASK_DIR}/premask_report.txt"

SNP_ALN="${GUBBINS_DIR}/tpa_ss14.filtered_polymorphic_sites.fasta"
UNROOTED_TREE="${IQTREE_DIR}/tpa_ss14.treefile"
ROOTED_TREE="${IQTREE_DIR}/tpa_ss14_rooted.treefile"
FCONST_FILE="${IQTREE_DIR}/fconst.txt"
FCONST_REPORT="${IQTREE_DIR}/fconst_report.txt"

mkdir -p "${MASK_DIR}" "${GUBBINS_DIR}" "${IQTREE_DIR}" "${TMP_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; exit 1; }

activate_env() {
    local env_name="$1"
    if [[ -z "${env_name}" || "${env_name}" == "none" ]]; then
        return 0
    fi
    if ! command -v conda >/dev/null 2>&1; then
        die "conda not found; cannot activate environment: ${env_name}"
    fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

resolve_input_msa() {
    if [[ -n "${INPUT_MSA}" ]]; then
        [[ -f "${INPUT_MSA}" ]] || die "INPUT_MSA not found: ${INPUT_MSA}"
        echo "${INPUT_MSA}"
        return 0
    fi
    local candidates=(
        "msa.fasta"
        "consensus/TPA_ss14ref_msa.fasta"
        "consensus/msa.fasta"
        "results/00_msa/msa.fasta"
    )
    local f
    for f in "${candidates[@]}"; do
        if [[ -f "${f}" ]]; then
            echo "${f}"
            return 0
        fi
    done
    die "Input MSA not found. Set INPUT_MSA=/path/to/msa.fasta or place msa.fasta in the project root."
}

MSA_FILE=$(resolve_input_msa)
[[ -f "${REFERENCE}" ]] || die "Reference FASTA not found: ${REFERENCE}"

log "Input MSA: ${MSA_FILE}"
log "Reference FASTA: ${REFERENCE}"

# Step 1: Run sample-level QC and pre-mask hypervariable regions.
# Input:  whole-genome MSA and SS14 reference FASTA
# Output: QC tables, mask BED, pre-mask report, and masked alignment
if [[ "${FORCE}" == "1" || ! -s "${MASKED_FASTA}" ]]; then
    log "Step 1: QC and pre-mask alignment."
    export MSA_FILE REFERENCE QC_SUMMARY PASSED_IDS FAILED_IDS MASKED_FASTA PREMASK_BED PREMASK_REPORT MIN_COVERAGE
    python3 <<'PYEOF'
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

msa_file = Path(os.environ["MSA_FILE"])
ref_file = Path(os.environ["REFERENCE"])
qc_summary = Path(os.environ["QC_SUMMARY"])
passed_ids = Path(os.environ["PASSED_IDS"])
failed_ids = Path(os.environ["FAILED_IDS"])
masked_fasta = Path(os.environ["MASKED_FASTA"])
premask_bed = Path(os.environ["PREMASK_BED"])
premask_report = Path(os.environ["PREMASK_REPORT"])
min_coverage = float(os.environ.get("MIN_COVERAGE", "0.75"))

regions = [
    (8341, 10165, "tprA", "TPASS_RS00040"),
    (10397, 12379, "tprB", "TPASS_RS00045"),
    (134911, 136707, "tprC", "TPASS_RS00590"),
    (152370, 154160, "tprD", "TPASS_RS00675"),
    (329182, 331470, "tprE", "TPASS_RS01555"),
    (332337, 333531, "tprF", "TPASS_RS01565"),
    (333590, 335860, "tprG", "TPASS_RS01570"),
    (461728, 463542, "arp", "TPASS_RS02125"),
    (498768, 499709, "TPANIC_0470", "TPASS_RS02290"),
    (593141, 594457, "TPANIC_0548", "TPASS_RS02700"),
    (662728, 664809, "tprH", "TPASS_RS03020"),
    (672540, 674369, "tprI", "TPASS_RS03070"),
    (674428, 676704, "tprJ", "TPASS_RS03075"),
    (934605, 935789, "TPANIC_0856", "TPASS_RS04240"),
    (935891, 937117, "TPANIC_0858", "TPASS_RS04245"),
    (937199, 938683, "TPANIC_0859", "TPASS_RS04250"),
    (944847, 946289, "TPANIC_0865", "TPASS_RS04275"),
    (975797, 977377, "tprK", "TPASS_RS05385"),
    (1125912, 1127456, "tprL", "TPASS_RS05110"),
]

for path in [qc_summary.parent, masked_fasta.parent, premask_bed.parent, premask_report.parent]:
    path.mkdir(parents=True, exist_ok=True)

ref_rec = next(SeqIO.parse(str(ref_file), "fasta"))
ref_seq = str(ref_rec.seq).upper()
ref_len = len(ref_seq) - ref_seq.count("N") - ref_seq.count("-")
if ref_len <= 0:
    raise SystemExit("Reference effective length is zero")

passed = []
failed = []
rows = []
for rec in SeqIO.parse(str(msa_file), "fasta"):
    sid = rec.id
    seq = str(rec.seq).upper()
    n_count = seq.count("N") + seq.count("-")
    valid_bp = len(seq) - n_count
    coverage = valid_bp / ref_len
    pass_qc = coverage >= min_coverage
    rows.append((sid, len(seq), valid_bp, n_count, round(coverage, 6), pass_qc))
    if pass_qc:
        passed.append(sid)
    else:
        failed.append(sid)

with qc_summary.open("w") as handle:
    handle.write("sample\ttotal_bp\tvalid_bp\tN_count\tcoverage\tpass_QC\n")
    for row in rows:
        handle.write("\t".join(map(str, row)) + "\n")

passed_ids.write_text("\n".join(passed) + ("\n" if passed else ""))
failed_ids.write_text("\n".join(failed) + ("\n" if failed else ""))

with premask_bed.open("w") as handle:
    handle.write("# chrom\tstart_0based\tend\tgene\tlocus_tag\n")
    for start, end, gene, locus in regions:
        handle.write(f"NC_021508.1\t{start - 1}\t{end}\t{gene}\t{locus}\n")

passed_set = set(passed)
masked_records = []
skipped = 0
for rec in SeqIO.parse(str(msa_file), "fasta"):
    if rec.id not in passed_set:
        skipped += 1
        continue
    seq = list(str(rec.seq).upper())
    for start, end, gene, locus in regions:
        idx_start = start - 1
        idx_end = min(end, len(seq))
        for i in range(idx_start, idx_end):
            seq[i] = "N"
    masked_records.append(SeqRecord(Seq("".join(seq)), id=rec.id, description=""))

if not masked_records:
    raise SystemExit("No sequences passed QC")

SeqIO.write(masked_records, str(masked_fasta), "fasta")

genome_len = len(masked_records[0].seq)
total_masked_bp = sum(end - start + 1 for start, end, _, _ in regions)
with premask_report.open("w") as handle:
    handle.write("Pre-mask report\n")
    handle.write(f"Input MSA\t{msa_file}\n")
    handle.write(f"Reference\t{ref_file}\n")
    handle.write(f"Minimum coverage\t{min_coverage}\n")
    handle.write(f"Sequences in MSA\t{len(rows)}\n")
    handle.write(f"Sequences written\t{len(masked_records)}\n")
    handle.write(f"Sequences skipped\t{skipped}\n")
    handle.write(f"Genome length\t{genome_len}\n")
    handle.write(f"Total masked bp\t{total_masked_bp}\n")
    handle.write("\nRegion\tStart\tEnd\tSize_bp\tLocus_tag\n")
    for start, end, gene, locus in regions:
        handle.write(f"{gene}\t{start}\t{end}\t{end - start + 1}\t{locus}\n")

print(f"[QC] samples={len(rows)} passed={len(passed)} failed={len(failed)} min_coverage={min_coverage}")
print(f"[MASK] written={len(masked_records)} output={masked_fasta}")
PYEOF
else
    log "Step 1: masked alignment exists; skipping. Use FORCE=1 to rerun."
fi

# Step 2: Run Gubbins on the masked whole-genome alignment.
# Input:  results/01_masked_alignment/tpa_ss14_masked.fasta
# Output: results/02_gubbins/tpa_ss14.filtered_polymorphic_sites.fasta
if [[ "${FORCE}" == "1" || ! -s "${SNP_ALN}" ]]; then
    log "Step 2: Gubbins recombination filtering."
    (
        activate_env "${GUBBINS_ENV}"
        export TMPDIR="${TMP_DIR}"
        cd "${GUBBINS_DIR}"
        run_gubbins.py \
            --prefix tpa_ss14 \
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
    log "Step 2: Gubbins SNP alignment exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${SNP_ALN}" ]] || die "Missing Gubbins SNP alignment: ${SNP_ALN}"

# Step 3: Count constant sites for IQ-TREE ascertainment correction.
# Input:  refs/TPA_SS14.fasta
# Output: results/03_iqtree/fconst.txt
if [[ "${FORCE}" == "1" || ! -s "${FCONST_FILE}" ]]; then
    log "Step 3: Constant-site counts for IQ-TREE."
    export REFERENCE FCONST_FILE FCONST_REPORT
    python3 <<'PYEOF'
import os
from pathlib import Path
from Bio import SeqIO

ref_file = Path(os.environ["REFERENCE"])
fconst_file = Path(os.environ["FCONST_FILE"])
fconst_report = Path(os.environ["FCONST_REPORT"])
regions = [
    (8341, 10165, "tprA", "TPASS_RS00040"),
    (10397, 12379, "tprB", "TPASS_RS00045"),
    (134911, 136707, "tprC", "TPASS_RS00590"),
    (152370, 154160, "tprD", "TPASS_RS00675"),
    (329182, 331470, "tprE", "TPASS_RS01555"),
    (332337, 333531, "tprF", "TPASS_RS01565"),
    (333590, 335860, "tprG", "TPASS_RS01570"),
    (461728, 463542, "arp", "TPASS_RS02125"),
    (498768, 499709, "TPANIC_0470", "TPASS_RS02290"),
    (593141, 594457, "TPANIC_0548", "TPASS_RS02700"),
    (662728, 664809, "tprH", "TPASS_RS03020"),
    (672540, 674369, "tprI", "TPASS_RS03070"),
    (674428, 676704, "tprJ", "TPASS_RS03075"),
    (934605, 935789, "TPANIC_0856", "TPASS_RS04240"),
    (935891, 937117, "TPANIC_0858", "TPASS_RS04245"),
    (937199, 938683, "TPANIC_0859", "TPASS_RS04250"),
    (944847, 946289, "TPANIC_0865", "TPASS_RS04275"),
    (975797, 977377, "tprK", "TPASS_RS05385"),
    (1125912, 1127456, "tprL", "TPASS_RS05110"),
]
rec = next(SeqIO.parse(str(ref_file), "fasta"))
seq = list(str(rec.seq).upper())
genome_len = len(seq)
for start, end, gene, locus in regions:
    for i in range(start - 1, min(end, genome_len)):
        seq[i] = "N"
seq = "".join(seq)
counts = {base: seq.count(base) for base in "ACGT"}
fconst = f"{counts['A']},{counts['C']},{counts['G']},{counts['T']}"
fconst_file.parent.mkdir(parents=True, exist_ok=True)
fconst_file.write_text(fconst + "\n")
with fconst_report.open("w") as handle:
    handle.write("Constant-site counts for IQ-TREE\n")
    handle.write(f"Reference\t{ref_file}\n")
    handle.write(f"Genome length\t{genome_len}\n")
    handle.write(f"A\t{counts['A']}\n")
    handle.write(f"C\t{counts['C']}\n")
    handle.write(f"G\t{counts['G']}\n")
    handle.write(f"T\t{counts['T']}\n")
    handle.write(f"fconst\t{fconst}\n")
print(f"[FCONST] {fconst}")
PYEOF
else
    log "Step 3: fconst file exists; skipping. Use FORCE=1 to rerun."
fi

# Step 4: Infer the ML tree from the Gubbins SNP alignment.
# Input:  results/02_gubbins/tpa_ss14.filtered_polymorphic_sites.fasta
# Output: results/03_iqtree/tpa_ss14.treefile
if [[ "${FORCE}" == "1" || ! -s "${UNROOTED_TREE}" ]]; then
    log "Step 4: IQ-TREE ML phylogeny."
    (
        activate_env "${TP_WGS_ENV}"
        IQTREE_BIN=${IQTREE_BIN:-}
        if [[ -z "${IQTREE_BIN}" ]]; then
            if command -v iqtree2 >/dev/null 2>&1; then
                IQTREE_BIN="iqtree2"
            elif command -v iqtree >/dev/null 2>&1; then
                IQTREE_BIN="iqtree"
            else
                die "iqtree2 or iqtree not found"
            fi
        fi
        FCONST=$(tr -d '[:space:]' < "${PROJECT_ROOT}/${FCONST_FILE}")
        "${IQTREE_BIN}" \
            -s "${PROJECT_ROOT}/${SNP_ALN}" \
            -m GTR+G4 \
            -fconst "${FCONST}" \
            -bb 1000 \
            -nt "${THREADS}" \
            -pre "${PROJECT_ROOT}/${IQTREE_DIR}/tpa_ss14" \
            -redo \
            > "${PROJECT_ROOT}/${IQTREE_DIR}/iqtree.log" 2>&1
    )
else
    log "Step 4: IQ-TREE tree exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${UNROOTED_TREE}" ]] || die "Missing IQ-TREE treefile: ${UNROOTED_TREE}"

# Step 5: Midpoint-root the ML tree.
# Input:  results/03_iqtree/tpa_ss14.treefile
# Output: results/03_iqtree/tpa_ss14_rooted.treefile
if [[ "${FORCE}" == "1" || ! -s "${ROOTED_TREE}" ]]; then
    log "Step 5: Midpoint-root ML tree."
    export UNROOTED_TREE ROOTED_TREE
    (
        activate_env "${TP_WGS_ENV}"
        Rscript --vanilla - <<'REOF'
suppressPackageStartupMessages({
    library(ape)
    library(phytools)
})
infile <- Sys.getenv("UNROOTED_TREE")
outfile <- Sys.getenv("ROOTED_TREE")
tree <- read.tree(infile)
rooted <- midpoint.root(tree)
write.tree(rooted, file = outfile)
cat(sprintf("[ROOT] wrote %s\n", outfile))
REOF
    )
else
    log "Step 5: rooted tree exists; skipping. Use FORCE=1 to rerun."
fi
[[ -s "${ROOTED_TREE}" ]] || die "Missing rooted tree: ${ROOTED_TREE}"

log "p1 completed."
log "Masked alignment: ${MASKED_FASTA}"
log "Gubbins SNP alignment: ${SNP_ALN}"
log "Unrooted tree: ${UNROOTED_TREE}"
log "Rooted tree: ${ROOTED_TREE}"
