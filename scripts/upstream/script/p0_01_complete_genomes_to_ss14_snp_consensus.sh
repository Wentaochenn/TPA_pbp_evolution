#!/usr/bin/env bash
# p0_01_complete_genomes_to_ss14_snp_consensus.sh
#
# Generate SS14-coordinate SNP-only consensus FASTA files from complete genome FASTA inputs,
# then combine all SS14-coordinate SNP-only consensus FASTA files into one MSA.
#
# Inputs:
#   config/complete_genome_manifest.tsv
#   refs/TPA_SS14.fasta
#   consensus/ss14/*.ss14.snp_only.consensus.fasta
#   consensus/nichols_align2_ss14/*.nichols_a2s.snp_only.consensus.fasta
#
# Outputs:
#   consensus/complete_genome_align2_ss14/bam/<sample>_vs_ss14.bam
#   consensus/complete_genome_align2_ss14/vcf/<sample>_vs_ss14.vcf.gz
#   consensus/TPA_SS14_snp_only_consensus/<sample>_snp_consensus.fasta
#   consensus/TPA_SS14_snp_only_consensus/complete_genome_snp_consensus_for_msa.fasta
#   consensus/TPA_SS14_snp_only_consensus/complete_genome_snp_consensus_qc.tsv
#   consensus/TPA_ss14ref_msa.fasta
#   consensus/TPA_ss14ref_msa_sources.tsv
#
# Usage:
#   bash script/p0_01_complete_genomes_to_ss14_snp_consensus.sh
#   MANIFEST=config/complete_genome_manifest.tsv bash script/p0_01_complete_genomes_to_ss14_snp_consensus.sh
#   THREADS=8 FORCE=1 bash script/p0_01_complete_genomes_to_ss14_snp_consensus.sh

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$(pwd)}
REF=${REF:-"${PROJECT_ROOT}/refs/TPA_SS14.fasta"}
MANIFEST=${MANIFEST:-"${PROJECT_ROOT}/config/complete_genome_manifest.tsv"}
OUTBASE=${OUTBASE:-"${PROJECT_ROOT}/consensus/complete_genome_align2_ss14"}
CONSENSUS_DIR=${CONSENSUS_DIR:-"${PROJECT_ROOT}/consensus/TPA_SS14_snp_only_consensus"}
MSA_OUT=${MSA_OUT:-"${PROJECT_ROOT}/consensus/TPA_ss14ref_msa.fasta"}
MSA_SOURCE_TABLE=${MSA_SOURCE_TABLE:-"${PROJECT_ROOT}/consensus/TPA_ss14ref_msa_sources.tsv"}
REFERENCE_SAMPLE_ID=${REFERENCE_SAMPLE_ID:-"NC_021508.1"}
THREADS=${THREADS:-8}
MAPQ=${MAPQ:-20}
BASEQ=${BASEQ:-20}
FORCE=${FORCE:-0}
TMPDIR=${TMPDIR:-"${PROJECT_ROOT}/tmp"}
BACKUP_EXISTING_MSA=${BACKUP_EXISTING_MSA:-1}

BAM_DIR="${OUTBASE}/bam"
VCF_DIR="${OUTBASE}/vcf"
LOG_DIR="${OUTBASE}/logs"
QC_FILE="${CONSENSUS_DIR}/complete_genome_snp_consensus_qc.tsv"
COMPLETE_MSA_FASTA="${CONSENSUS_DIR}/complete_genome_snp_consensus_for_msa.fasta"
SIMPLE_REF="${TMPDIR}/TPA_SS14_simple.fasta"
SOURCE_LIST="${TMPDIR}/ss14_msa_source_files.tsv"
REFERENCE_ROW="${TMPDIR}/${REFERENCE_SAMPLE_ID}_ss14_reference_row.fasta"

mkdir -p "${BAM_DIR}" "${VCF_DIR}" "${LOG_DIR}" "${CONSENSUS_DIR}" "${TMPDIR}" "$(dirname "${MSA_OUT}")"

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

require_file() {
    local path="$1"
    local label="$2"
    [[ -f "${path}" ]] || die "Missing ${label}: ${path}"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Executable not found: $1"
}

fasta_len() {
    grep -v '^>' "$1" | tr -d '\n\r\t ' | wc -c | awk '{print $1}'
}

fasta_md5() {
    grep -v '^>' "$1" | tr -d '\n\r\t ' | tr '[:lower:]' '[:upper:]' | md5sum | awk '{print $1}'
}

write_qc_row() {
    local sample="$1"
    local fasta="$2"
    local role="$3"
    local bam="$4"
    local vcfgz="$5"
    local outfa="$6"
    local input_len="$7"
    local consensus_len="$8"
    local expected_len="$9"
    local n_variants="${10}"
    local status="${11}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sample}" "${fasta}" "${role}" "${bam}" "${vcfgz}" "${outfa}" \
        "${input_len}" "${consensus_len}" "${expected_len}" "${n_variants}" "${status}" \
        >> "${QC_FILE}"
}

write_reference_row() {
    sed 's/ .*//' "${REF}" | awk -v id="${REFERENCE_SAMPLE_ID}" '
        BEGIN { printed=0 }
        /^>/ { if (printed == 0) { print ">" id; printed=1 } ; next }
        { print toupper($0) }
    ' > "${REFERENCE_ROW}"
}

append_source_if_exists() {
    local source_type="$1"
    local file_path="$2"
    [[ -f "${file_path}" ]] || return 0
    printf '%s\t%s\n' "${source_type}" "${file_path}" >> "${SOURCE_LIST}"
}

append_sources_by_glob() {
    local source_type="$1"
    shift
    local pattern f
    for pattern in "$@"; do
        while IFS= read -r -d '' f; do
            append_source_if_exists "${source_type}" "${f}"
        done < <(find ${pattern%/*} -maxdepth 1 -type f -name "${pattern##*/}" -print0 2>/dev/null | sort -z)
    done
}

require_file "${REF}" "SS14 reference FASTA"
require_file "${MANIFEST}" "complete genome manifest"

for cmd in minimap2 samtools bcftools bgzip tabix awk sed grep tr wc md5sum python3; do
    require_cmd "${cmd}"
done

sed 's/ .*//' "${REF}" > "${SIMPLE_REF}"
REF_LEN=$(fasta_len "${SIMPLE_REF}")
REF_MD5=$(fasta_md5 "${SIMPLE_REF}")
[[ "${REF_LEN}" -gt 0 ]] || die "Reference length is zero: ${SIMPLE_REF}"

write_reference_row
REFERENCE_LEN=$(fasta_len "${REFERENCE_ROW}")
[[ "${REFERENCE_LEN}" == "${REF_LEN}" ]] || die "Reference row length mismatch: ${REFERENCE_ROW}"

log "Project root: ${PROJECT_ROOT}"
log "Reference FASTA: ${REF}"
log "Manifest: ${MANIFEST}"
log "Reference sample ID: ${REFERENCE_SAMPLE_ID}"
log "Reference length: ${REF_LEN}"
log "Consensus output: ${CONSENSUS_DIR}"
log "MSA output: ${MSA_OUT}"

printf 'sample_id\tinput_fasta\trole\tbam\tvcf_gz\tconsensus_fasta\tinput_length\tconsensus_length\texpected_length\tn_variants\tstatus\n' > "${QC_FILE}"
: > "${COMPLETE_MSA_FASTA}"
: > "${SOURCE_LIST}"

# Manifest columns:
#   sample_id<TAB>fasta<TAB>role
#
# The role column is optional.
# Accepted reference roles: reference, ref, ss14_reference.
# Accepted query roles: query, auto, include, complete_genome, blank.
# Lines beginning with # are ignored.
# A header line with sample_id and fasta is accepted.
awk -F'\t' 'BEGIN{OFS="\t"} NF>=2 && $1 !~ /^#/ {role=(NF>=3 ? $3 : "auto"); print $1, $2, role}' "${MANIFEST}" |
while IFS=$'\t' read -r sample fasta role; do
    [[ -n "${sample}" ]] || continue
    [[ "${sample}" == "sample_id" ]] && continue
    role=${role:-auto}
    role_lc=$(echo "${role}" | tr '[:upper:]' '[:lower:]')

    if [[ "${fasta}" != /* ]]; then
        fasta="${PROJECT_ROOT}/${fasta}"
    fi
    require_file "${fasta}" "complete genome FASTA for ${sample}"

    INPUT_LEN=$(fasta_len "${fasta}")
    INPUT_MD5=$(fasta_md5 "${fasta}")

    BAM="${BAM_DIR}/${sample}_vs_ss14.bam"
    VCF="${VCF_DIR}/${sample}_vs_ss14.vcf"
    VCFGZ="${VCF}.gz"
    OUTFA="${CONSENSUS_DIR}/${sample}_snp_consensus.fasta"
    SAMPLE_LOG="${LOG_DIR}/${sample}.log"

    log "Processing ${sample}"
    log "  Input length: ${INPUT_LEN}; MD5: ${INPUT_MD5}; role=${role_lc}"

    if [[ "${role_lc}" == "reference" || "${role_lc}" == "ref" || "${role_lc}" == "ss14_reference" ]]; then
        log "  Using ${sample} as the SS14 reference row in the MSA"
        REFERENCE_SAMPLE_ID="${sample}"
        sed 's/ .*//' "${REF}" | awk -v id="${sample}" '
            BEGIN { printed=0 }
            /^>/ { if (printed == 0) { print ">" id; printed=1 } ; next }
            { print toupper($0) }
        ' > "${REFERENCE_ROW}"
        write_qc_row "${sample}" "${fasta}" "${role_lc}" "NA" "NA" "${REFERENCE_ROW}" "${INPUT_LEN}" "${REF_LEN}" "${REF_LEN}" "0" "REFERENCE_ROW"
        continue
    fi

    if [[ "${INPUT_MD5}" == "${REF_MD5}" ]]; then
        log "  Input FASTA is identical to the SS14 reference; using ${sample} as the reference row"
        sed 's/ .*//' "${REF}" | awk -v id="${sample}" '
            BEGIN { printed=0 }
            /^>/ { if (printed == 0) { print ">" id; printed=1 } ; next }
            { print toupper($0) }
        ' > "${REFERENCE_ROW}"
        write_qc_row "${sample}" "${fasta}" "auto_reference" "NA" "NA" "${REFERENCE_ROW}" "${INPUT_LEN}" "${REF_LEN}" "${REF_LEN}" "0" "REFERENCE_ROW"
        continue
    fi

    if [[ ! -f "${BAM}" || "${FORCE}" == "1" ]]; then
        log "  Aligning complete genome FASTA to SS14 reference"
        minimap2 -a --cs -t "${THREADS}" "${REF}" "${fasta}" \
            2> "${SAMPLE_LOG}.minimap2" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
        samtools index "${BAM}"
    else
        log "  BAM exists; skipping alignment"
    fi

    if [[ ! -f "${VCFGZ}" || "${FORCE}" == "1" ]]; then
        log "  Calling SNPs"
        bcftools mpileup \
            -f "${REF}" \
            -q "${MAPQ}" \
            -Q "${BASEQ}" \
            "${BAM}" \
            | bcftools call -mv --ploidy 1 -o "${VCF}"
        bgzip -f "${VCF}"
        tabix -f -p vcf "${VCFGZ}"
    else
        log "  VCF exists; skipping SNP calling"
    fi

    N_VARIANTS=$(bcftools view -H "${VCFGZ}" | wc -l | awk '{print $1}')

    if [[ ! -f "${OUTFA}" || "${FORCE}" == "1" ]]; then
        log "  Generating SNP-only consensus FASTA"
        bcftools consensus \
            -f "${SIMPLE_REF}" \
            "${VCFGZ}" \
            | sed "s/^>.*/>${sample}/" \
            > "${OUTFA}"
    else
        log "  Consensus FASTA exists; skipping consensus generation"
    fi

    LEN=$(fasta_len "${OUTFA}")
    STATUS="PASS"
    if [[ "${LEN}" != "${REF_LEN}" ]]; then
        STATUS="FAIL_LENGTH"
    fi

    cat "${OUTFA}" >> "${COMPLETE_MSA_FASTA}"
    write_qc_row "${sample}" "${fasta}" "${role_lc}" "${BAM}" "${VCFGZ}" "${OUTFA}" "${INPUT_LEN}" "${LEN}" "${REF_LEN}" "${N_VARIANTS}" "${STATUS}"

    log "  Variants: ${N_VARIANTS}; length: ${LEN} (expected ${REF_LEN}); status=${STATUS}"
done

if grep -q $'\tFAIL_LENGTH$' "${QC_FILE}"; then
    die "One or more complete-genome consensus FASTA files failed the length check: ${QC_FILE}"
fi

# Step 2: Build the SS14-coordinate MSA from all SNP-only consensus FASTA files.
append_source_if_exists "ss14_reference" "${REFERENCE_ROW}"
append_sources_by_glob "ss14_wgs" "${PROJECT_ROOT}/consensus/ss14/*.ss14.snp_only.consensus.fasta"
append_sources_by_glob "nichols_align2_ss14" "${PROJECT_ROOT}/consensus/nichols_align2_ss14/*.nichols_a2s.snp_only.consensus.fasta"
append_sources_by_glob "complete_genome_align2_ss14" "${CONSENSUS_DIR}/*_snp_consensus.fasta"

if [[ ! -s "${SOURCE_LIST}" ]]; then
    die "No SS14-coordinate SNP-only consensus FASTA files were found for MSA construction"
fi

if [[ -f "${MSA_OUT}" && "${BACKUP_EXISTING_MSA}" == "1" ]]; then
    BACKUP="${MSA_OUT}.bak.$(date '+%Y%m%d_%H%M%S')"
    cp -p "${MSA_OUT}" "${BACKUP}"
    log "Existing MSA backed up: ${BACKUP}"
fi

python3 - "${SOURCE_LIST}" "${MSA_OUT}" "${MSA_SOURCE_TABLE}" "${REF_LEN}" <<'PY'
import sys
from pathlib import Path

source_list = Path(sys.argv[1])
msa_out = Path(sys.argv[2])
source_table = Path(sys.argv[3])
expected_len = int(sys.argv[4])

def read_fasta(path):
    path = Path(path)
    records = []
    header = None
    seq_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts).upper()))
    if len(records) != 1:
        raise SystemExit(f"Expected one FASTA record in {path}, found {len(records)}")
    return records[0]

rows = []
seen_ids = set()
with source_list.open() as handle:
    for raw in handle:
        raw = raw.rstrip("\n")
        if not raw:
            continue
        source_type, fasta = raw.split("\t", 1)
        sample_id, seq = read_fasta(fasta)
        status = "PASS"
        if sample_id in seen_ids:
            status = "SKIP_DUPLICATE_SAMPLE_ID"
        elif len(seq) != expected_len:
            status = "FAIL_LENGTH"
        else:
            seen_ids.add(sample_id)
        rows.append({
            "sample_id": sample_id,
            "source_type": source_type,
            "fasta": fasta,
            "length": len(seq),
            "expected_length": expected_len,
            "status": status,
            "sequence": seq,
        })

failures = [r for r in rows if r["status"] == "FAIL_LENGTH"]
if failures:
    for r in failures:
        print(f"Length failure: {r['sample_id']} {r['length']} expected {expected_len} from {r['fasta']}", file=sys.stderr)
    raise SystemExit("MSA construction stopped because one or more sequences failed length checks")

with msa_out.open("w") as out:
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
        out.write(
            f"{r['sample_id']}\t{r['source_type']}\t{r['fasta']}\t{r['length']}\t{r['expected_length']}\t{r['status']}\n"
        )

n_pass = sum(1 for r in rows if r["status"] == "PASS")
print(f"MSA records written: {n_pass}")
print(f"MSA output: {msa_out}")
print(f"Source table: {source_table}")
PY

log "Complete-genome QC table: ${QC_FILE}"
log "Complete-genome FASTA for MSA: ${COMPLETE_MSA_FASTA}"
log "MSA source table: ${MSA_SOURCE_TABLE}"
log "MSA output: ${MSA_OUT}"
log "Finished. Next step: bash master.sh p1"
