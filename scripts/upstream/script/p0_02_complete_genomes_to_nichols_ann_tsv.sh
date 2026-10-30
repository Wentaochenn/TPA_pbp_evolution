#!/usr/bin/env bash
# p0_02_complete_genomes_to_nichols_ann_tsv.sh
#
# Generate Nichols-coordinate snpEff annotation tables for complete genome FASTA files.
#
# Inputs:
#   config/complete_genome_manifest.tsv
#   refs/TPA_Nichols.fasta
#
# Outputs:
#   results/04_complete_genome_variant_calls/bam/<sample>_vs_nichols.bam
#   results/04_complete_genome_variant_calls/vcf/<sample>_vs_nichols.pass.norm.vcf.gz
#   results/04_complete_genome_variant_calls/ann_vcf/<sample>_vs_nichols.pass.ann.vcf.gz
#   results/04_complete_genome_variant_calls/ann_tsv/<sample>.ann.tsv
#   results/04_complete_genome_variant_calls/complete_genome_nichols_annotation_qc.tsv
#
# Usage:
#   bash script/p0_02_complete_genomes_to_nichols_ann_tsv.sh

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$(pwd)}
MANIFEST=${MANIFEST:-"config/complete_genome_manifest.tsv"}
REFERENCE=${REFERENCE:-"refs/TPA_Nichols.fasta"}
OUTDIR=${OUTDIR:-"results/04_complete_genome_variant_calls"}
THREADS=${THREADS:-24}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_BASEQ=${MIN_BASEQ:-20}
SNPEFF_DB=${SNPEFF_DB:-"TPA_Nichols"}
SNPEFF_ENV=${SNPEFF_ENV:-"tpa_main_env"}
FORCE=${FORCE:-0}

cd "${PROJECT_ROOT}"

BAM_DIR="${OUTDIR}/bam"
VCF_DIR="${OUTDIR}/vcf"
ANN_VCF_DIR="${OUTDIR}/ann_vcf"
ANN_TSV_DIR="${OUTDIR}/ann_tsv"
LOG_DIR="${OUTDIR}/logs"
QC_TABLE="${OUTDIR}/complete_genome_nichols_annotation_qc.tsv"
TMP_DIR="${PROJECT_ROOT}/tmp"

mkdir -p "${BAM_DIR}" "${VCF_DIR}" "${ANN_VCF_DIR}" "${ANN_TSV_DIR}" "${LOG_DIR}" "${TMP_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; exit 1; }

activate_env() {
    local env_name="$1"
    if [[ -z "${env_name}" || "${env_name}" == "none" ]]; then
        return 0
    fi
    command -v conda >/dev/null 2>&1 || die "conda not found; cannot activate environment: ${env_name}"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    hash -r
}

run_snpeff() {
    local norm_vcf="$1"
    local ann_vcf="$2"
    local html="$3"
    local csv="$4"
    (
        activate_env "${SNPEFF_ENV}"
        snpEff -Xmx8g \
            -v "${SNPEFF_DB}" \
            -stats "${html}" \
            -csvStats "${csv}" \
            "${norm_vcf}" \
        | "${BGZIP_BIN}" -c > "${ann_vcf}"
    )
}

run_snpsift() {
    local ann_vcf="$1"
    local ann_tsv="$2"
    (
        activate_env "${SNPEFF_ENV}"
        SnpSift extractFields -s ',' -e '.' \
            "${ann_vcf}" \
            CHROM POS REF ALT \
            'ANN[*].GENE' \
            'ANN[*].GENEID' \
            'ANN[*].EFFECT' \
            'ANN[*].IMPACT' \
            'ANN[*].HGVS_C' \
            'ANN[*].HGVS_P' \
            'ANN[*].RANK' \
            GEN[0].DP GEN[0].GT \
        > "${ann_tsv}"
    )
}

fasta_length() {
    awk 'BEGIN{n=0} /^>/{next} {gsub(/[[:space:]]/, ""); n += length($0)} END{print n+0}' "$1"
}

variant_count() {
    local vcf="$1"
    if [[ ! -s "${vcf}" ]]; then
        echo 0
        return 0
    fi
    bcftools view -H "${vcf}" 2>/dev/null | wc -l | awk '{print $1}'
}

write_empty_ann_tsv() {
    local ann_tsv="$1"
    printf "CHROM\tPOS\tREF\tALT\tANN[*].GENE\tANN[*].GENEID\tANN[*].EFFECT\tANN[*].IMPACT\tANN[*].HGVS_C\tANN[*].HGVS_P\tANN[*].RANK\tGEN[0].DP\tGEN[0].GT\n" > "${ann_tsv}"
}

[[ -f "${MANIFEST}" ]] || die "Complete genome manifest not found: ${MANIFEST}"
[[ -f "${REFERENCE}" ]] || die "Nichols reference FASTA not found: ${REFERENCE}"

for exe in minimap2 samtools bcftools bgzip tabix; do
    command -v "${exe}" >/dev/null 2>&1 || die "Required executable not found in current PATH: ${exe}"
done

BGZIP_BIN=$(command -v bgzip)
TABIX_BIN=$(command -v tabix)

if [[ ! -f "${REFERENCE}.fai" ]]; then
    log "Indexing Nichols reference with samtools faidx."
    samtools faidx "${REFERENCE}"
fi

if [[ "${FORCE}" == "1" || ! -f "${QC_TABLE}" ]]; then
    printf "sample_id\tfasta\trole\tfasta_length\tbam\traw_vcf\tnorm_vcf\tann_vcf\tann_tsv\traw_variant_count\tnorm_variant_count\tann_tsv_records\tstatus\n" > "${QC_TABLE}"
fi

log "Project root: ${PROJECT_ROOT}"
log "Manifest: ${MANIFEST}"
log "Nichols reference: ${REFERENCE}"
log "Output directory: ${OUTDIR}"
log "snpEff database: ${SNPEFF_DB}"
log "snpEff environment: ${SNPEFF_ENV}"

processed=0
skipped=0

while IFS=$'\t' read -r sample_id fasta role rest; do
    [[ -z "${sample_id}" ]] && continue
    [[ "${sample_id}" =~ ^# ]] && continue
    if [[ "${sample_id}" == "sample_id" ]]; then
        continue
    fi
    role=${role:-auto}
    role_lc=$(echo "${role}" | tr '[:upper:]' '[:lower:]')
    if [[ "${role_lc}" == "skip" || "${role_lc}" == "exclude" || "${role_lc}" == "ignore" ]]; then
        log "Skipping ${sample_id}: role=${role}"
        skipped=$((skipped + 1))
        continue
    fi

    [[ -f "${fasta}" ]] || die "Input FASTA not found for ${sample_id}: ${fasta}"

    sample_log="${LOG_DIR}/${sample_id}.log"
    bam="${BAM_DIR}/${sample_id}_vs_nichols.bam"
    raw_vcf="${VCF_DIR}/${sample_id}_vs_nichols.raw.vcf.gz"
    norm_vcf="${VCF_DIR}/${sample_id}_vs_nichols.pass.norm.vcf.gz"
    ann_vcf="${ANN_VCF_DIR}/${sample_id}_vs_nichols.pass.ann.vcf.gz"
    html="${ANN_VCF_DIR}/${sample_id}_snpEff_summary.html"
    csv="${ANN_VCF_DIR}/${sample_id}_snpEff_summary.csv"
    ann_tsv="${ANN_TSV_DIR}/${sample_id}.ann.tsv"

    if [[ "${FORCE}" != "1" && -s "${ann_tsv}" ]]; then
        log "${sample_id}: ann.tsv exists; skipping. Use FORCE=1 to rerun."
        skipped=$((skipped + 1))
        continue
    fi

    log "Processing ${sample_id}: ${fasta}"
    {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] sample_id=${sample_id}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] fasta=${fasta}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] reference=${REFERENCE}"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: minimap2 alignment to Nichols reference."
        minimap2 -a --cs -t "${THREADS}" "${REFERENCE}" "${fasta}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
        samtools index "${bam}"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: SNP calling."
        bcftools mpileup \
            -f "${REFERENCE}" \
            -q "${MIN_MAPQ}" \
            -Q "${MIN_BASEQ}" \
            --annotate FORMAT/AD,FORMAT/DP \
            "${bam}" \
        | bcftools call \
            -mv --ploidy 1 \
            -Oz -o "${raw_vcf}"
        bcftools index -t "${raw_vcf}"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: VCF normalization."
        bcftools norm \
            -f "${REFERENCE}" \
            -m -any \
            "${raw_vcf}" \
            -Oz -o "${norm_vcf}"
        bcftools index -t "${norm_vcf}"

        raw_n=$(variant_count "${raw_vcf}")
        norm_n=$(variant_count "${norm_vcf}")
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Raw variants: ${raw_n}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Normalized variants: ${norm_n}"

        if [[ "${norm_n}" -eq 0 ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 4: no variants; writing header-only ann.tsv."
            write_empty_ann_tsv "${ann_tsv}"
            : > "${ann_vcf}"
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 4: snpEff annotation."
            run_snpeff "${norm_vcf}" "${ann_vcf}" "${html}" "${csv}"
            "${TABIX_BIN}" -p vcf "${ann_vcf}"

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 5: extract annotation fields to TSV."
            run_snpsift "${ann_vcf}" "${ann_tsv}"
        fi

        ann_records=$(( $(wc -l < "${ann_tsv}") - 1 ))
        if [[ "${ann_records}" -lt 0 ]]; then ann_records=0; fi
        fa_len=$(fasta_length "${fasta}")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n" \
            "${sample_id}" "${fasta}" "${role}" "${fa_len}" "${bam}" "${raw_vcf}" "${norm_vcf}" "${ann_vcf}" "${ann_tsv}" "${raw_n}" "${norm_n}" "${ann_records}" \
            >> "${QC_TABLE}"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output ann.tsv: ${ann_tsv}"
    } 2>&1 | tee "${sample_log}"

    processed=$((processed + 1))
done < "${MANIFEST}"

log "p0_02 completed. processed=${processed} skipped=${skipped}"
log "Annotation TSV directory: ${ANN_TSV_DIR}"
log "QC table: ${QC_TABLE}"
