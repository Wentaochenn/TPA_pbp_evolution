#!/usr/bin/env bash
# master.sh
#
# Run the documented TPA upstream analysis stages from the project root.
# Each step explicitly changes to PROJECT_ROOT before executing commands.
#
# Usage:
#   bash master.sh check
#   bash master.sh all_until_beauti
#   bash master.sh all_after_beauti
#   bash master.sh all
#   bash master.sh figures

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-"$(pwd)"}
SCRIPT_DIR=${SCRIPT_DIR:-"${PROJECT_ROOT}/script"}
MASTER_LOG_DIR=${MASTER_LOG_DIR:-"${PROJECT_ROOT}/logs/master"}
THREADS=${THREADS:-24}
JOBS=${JOBS:-8}
BIG_MEM_JOBS=${BIG_MEM_JOBS:-2}
RESTART_TIMES=${RESTART_TIMES:-2}
SNAKEMAKE_BIN=${SNAKEMAKE_BIN:-snakemake}
PBP_ANALYSIS_PY=${PBP_ANALYSIS_PY:-python3}
PBP_MATRIX_FILE=${PBP_MATRIX_FILE:-results/05_extracted_variants/callability/sample_locus_ternary_matrix.tsv}
export PBP_MATRIX_FILE
INPUT_MSA=${INPUT_MSA:-""}

TP_WGS_ENV=${TP_WGS_ENV:-"tp_wgs"}
GUBBINS_ENV=${GUBBINS_ENV:-"gubbins"}
ESM2_ENV=${ESM2_ENV:-"evo2"}
SNPEFF_ENV=${SNPEFF_ENV:-"tpa_main_env"}
PBP_BBMAP_ENV=${PBP_BBMAP_ENV:-"pbp_bbmap_runtime"}
PBP_BBMAP_STARTUP_TIMEOUT=${PBP_BBMAP_STARTUP_TIMEOUT:-900}
export PBP_BBMAP_ENV PBP_BBMAP_STARTUP_TIMEOUT

TREETIME_MODE=${TREETIME_MODE:-"all_subsets"}
ESM2_MODEL_PATH=${ESM2_MODEL_PATH:-""}
ESM2_ALLOWED_HOST=${ESM2_ALLOWED_HOST:-"compute-0-3"}
N_PERMUTATIONS=${N_PERMUTATIONS:-10000}
BOOTSTRAP_N=${BOOTSTRAP_N:-2000}
RANDOM_SEED=${RANDOM_SEED:-42}
RUN_RELAXED_TREETIME=${RUN_RELAXED_TREETIME:-false}
FORCE=${FORCE:-0}
FORCE_MSA=${FORCE_MSA:-0}
FORCE_COMPLETE_ANN=${FORCE_COMPLETE_ANN:-0}

mkdir -p "${MASTER_LOG_DIR}"


PBP_PUBLICATION_PACKAGE=${PBP_PUBLICATION_PACKAGE:-"$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"}
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*" >&2; }

die() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; exit 1; }


usage() {
    cat <<USAGE
Usage:
  bash master.sh <command>

Commands:
  check                 Check scripts, conda environments, key inputs, and parser syntax.
  init                  Create expected output directories.

  wgs                   Reproduce both FASTQ reference branches and cross-reference variants through SGE.
  p0                    Convert complete genomes to SS14-coordinate SNP-only consensus and build consensus/TPA_ss14ref_msa.fasta.
  p0_02                 Generate Nichols-coordinate snpEff ann.tsv files for complete genomes.
  p1                    Run QC/pre-mask, Gubbins, IQ-TREE, and midpoint rooting.
  p2                    Extract codon-corrected variants and audit full-codon callability (0/1/NA).
  selection             Run HyPhy FUBAR selection tests.
  temporal              Run root-to-tip and date-randomization temporal screening.
  treetime              Run TreeTime on candidate subsets and select one subset.
  beast_prepare         Prepare BEAUti2/BEAST2 input files for the selected subset.
  beast_run             Run BEAST2 using BEAUti2-generated XML files.
  pbp_temporal          Compare PBP acquisition dates from both completed clock analyses.
  esm2                  Recompute four-protein scores in an allocated GPU job.
  structural            Generate current publication outputs, including reference structural views.
  figures               Generate manuscript figures.

  all_until_beauti      Run FASTQ through BEAST2 input preparation, stopping for BEAUti2 XML.
  all_after_beauti      Run BEAST2, PBP temporal dating, structural analysis, ESM-2, and figures.
  all                   Run the full pipeline; stops at BEAUti2 if XML is missing.
  all_from_p0           Run p0 explicitly, write the SS14-coordinate MSA, then stop for source-table review.

Environment variables:
  PROJECT_ROOT          ${PROJECT_ROOT}
  SCRIPT_DIR            ${SCRIPT_DIR}
  TP_WGS_ENV            ${TP_WGS_ENV}
  GUBBINS_ENV           ${GUBBINS_ENV}
  ESM2_ENV              ${ESM2_ENV}
  SNPEFF_ENV            ${SNPEFF_ENV}
  PBP_BBMAP_ENV         Dedicated BBMap 39.81/OpenJDK 17 environment. Default: ${PBP_BBMAP_ENV}
  PBP_BBMAP_STARTUP_TIMEOUT  Maximum seconds allowed for BBMap index loading. Default: ${PBP_BBMAP_STARTUP_TIMEOUT}
  ESM2_MODEL_PATH       ${ESM2_MODEL_PATH}
  TREETIME_MODE         recommended or all_subsets. Default: ${TREETIME_MODE}
  INPUT_MSA             Manually assembled SS14-coordinate SNP-only MSA. Default: consensus/TPA_ss14ref_msa.fasta.
  FORCE                 1 to rerun existing p1 outputs. Default: ${FORCE}
  FORCE_MSA             1 to rebuild consensus/TPA_ss14ref_msa.fasta during all/all_until_beauti. Default: ${FORCE_MSA}
  FORCE_COMPLETE_ANN    1 to regenerate complete-genome Nichols ann.tsv files before p2. Default: ${FORCE_COMPLETE_ANN}
  JOBS                  Maximum concurrent SGE jobs for each Snakemake stage. Default: ${JOBS}
  BIG_MEM_JOBS          Maximum concurrent 64-GB BBMap host-filtering jobs. Default: ${BIG_MEM_JOBS}
  RESTART_TIMES         Automatic retries for transient SGE/BBMap failures. Default: ${RESTART_TIMES}
  PBP_ANALYSIS_PY       Python with pandas and pysam for codon callability.
  PBP_MATRIX_FILE       The audited 0/1/NA matrix used downstream.
USAGE
}


cd_project() {
    [[ -d "${PROJECT_ROOT}" ]] || die "PROJECT_ROOT does not exist: ${PROJECT_ROOT}"
    cd "${PROJECT_ROOT}"
}


script_path() { echo "${SCRIPT_DIR}/$1"; }

require_script() { [[ -f "$1" ]] || die "Missing script: $1"; }

require_file() { [[ -f "$1" ]] || die "Missing $2: $1"; }


require_any_file() {
    local label="$1"; shift
    local f
    for f in "$@"; do
        if [[ -f "${f}" ]]; then
            log "Found ${label}: ${f}"
            return 0
        fi
    done
    die "Missing ${label}. Tried: $*"
}


conda_prefix() {
    command -v conda >/dev/null 2>&1 || die "conda not found"
    # Force the activated environment bin directory to the front of PATH.
    # Some cluster login shells prepend older R environments after conda activation.
    echo "source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate '$1' && export PATH=\"\${CONDA_PREFIX}/bin:\${PATH}\" && hash -r &&"
}


quote_for_bash() {
    printf "%q" "$1"
}


run_cmd() {
    local name="$1"
    local cmd="$2"
    local log_file="${MASTER_LOG_DIR}/${name}.log"
    local project_root_q
    project_root_q=$(quote_for_bash "${PROJECT_ROOT}")

    log "START ${name}"
    log "Command: cd ${PROJECT_ROOT} && ${cmd}"

    set +e
    (
        cd "${PROJECT_ROOT}" || exit 1
        bash -lc "cd ${project_root_q} && echo \"[MASTER_CWD] \$(pwd)\" && ${cmd}"
    ) 2>&1 | tee "${log_file}"
    local status=${PIPESTATUS[0]}
    set -e

    if [[ ${status} -ne 0 ]]; then
        die "Step failed: ${name}. Log: ${log_file}"
    fi
    log "DONE ${name}. Log: ${log_file}"
}


run_cmd_env() {
    local name="$1"
    local env_name="$2"
    local cmd="$3"
    run_cmd "${name}" "$(conda_prefix "${env_name}") cd \"${PROJECT_ROOT}\" && ${cmd}"
}


require_output() {
    local f="$1"
    local label="$2"
    local path="${f}"
    [[ "${path}" = /* ]] || path="${PROJECT_ROOT}/${path}"
    [[ -s "${path}" || -e "${path}" ]] || die "Expected output missing after step: ${label}: ${f}"
    log "Output OK: ${f}"
}


input_msa_path() {
    cd_project
    if [[ -n "${INPUT_MSA}" ]]; then
        [[ -f "${INPUT_MSA}" ]] || return 1
        echo "${INPUT_MSA}"
        return 0
    fi
    local candidates=(
        "consensus/TPA_ss14ref_msa.fasta"
        "msa.fasta"
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
    return 1
}


masked_alignment_path() {
    cd_project
    if [[ -f "results/01_masked_alignment/tpa_ss14_masked.fasta" ]]; then
        echo "results/01_masked_alignment/tpa_ss14_masked.fasta"
    elif [[ -f "results/01_masked/tpa_ss14_masked.fasta" ]]; then
        echo "results/01_masked/tpa_ss14_masked.fasta"
    else
        return 1
    fi
}


init_dirs() {
    cd_project
    mkdir -p \
        script config data refs VCF logs tmp \
        results/01_masked_alignment \
        results/02_gubbins \
        results/03_iqtree \
        results/04_complete_genome_variant_calls/bam \
        results/04_complete_genome_variant_calls/vcf \
        results/04_complete_genome_variant_calls/ann_vcf \
        results/04_complete_genome_variant_calls/ann_tsv \
        results/04_complete_genome_variant_calls/wgs_ann_tsv \
        results/04_complete_genome_variant_calls/logs \
        results/05_extracted_variants \
        results/06_asr \
        results/07_homoplasy \
        results/08_model_sensitivity \
        results/10_selection \
        results/11_temporal_signal \
        results/12_treetime \
        results/13_beast2_selected_subset/input \
        results/13_beast2_selected_subset/run \
        results/13_beast2_selected_subset/diagnostics \
        results/14_esm2_integration \
        results/15_structural
    log "Directory initialization finished."
}


check_scripts() {
    local scripts=(
        p0_00_bbmap_node_cache.sh
        p0_00_run_bbmap_guarded.sh
        p0_01_complete_genomes_to_ss14_snp_consensus.sh
        p0_02_complete_genomes_to_nichols_ann_tsv.sh
        p1_01_build_tree.sh
        p2_01_extract_codon_corrected_variants.py
        p3_01_asr_simmap_binary_loci.R
        p3_02_homoplasy_ci_ri.R
        p3_03_er_ard_model_sensitivity.R
        p3_05_hyphy_selection_tests.py
        p4_01_temporal_signal_screening.R
        p4_02_treetime_clock_screening.sh
        p4_03_prepare_beast2_inputs.sh
        p4_04_run_beast2_and_summarize.sh
        p2_02_audit_codon_callability.py
        p5_01_extract_reference_pbp_proteins.py
    )
    local s
    for s in "${scripts[@]}"; do require_script "$(script_path "${s}")"; done
    for s in p0_01_ss14_wgs.smk p0_02_nichols_wgs.smk p0_03_ss14_align2nichols.smk p0_04_nichols_align2ss14.smk; do
        require_script "${PROJECT_ROOT}/SC/${s}"
    done
    require_script "${PBP_PUBLICATION_PACKAGE}/master.sh"
    require_script "${PBP_PUBLICATION_PACKAGE}/analysis/clock_histories/run.py"
    require_script "${PBP_PUBLICATION_PACKAGE}/analysis/clock_histories/payload/date_origins.R"
    require_script "${SCRIPT_DIR}/pbp_mk_missing_helpers.R"
}

check_syntax() {
    cd_project
    local f
    for f in "${SCRIPT_DIR}"/*.sh; do
        [[ -f "${f}" ]] || continue
        bash -n "${f}"
        log "Bash syntax OK: ${f}"
    done
    for f in "${SCRIPT_DIR}"/*.py; do
        [[ -f "${f}" ]] || continue
        python3 -m py_compile "${f}"
        log "Python syntax OK: ${f}"
    done
    if command -v Rscript >/dev/null 2>&1; then
        for f in "${SCRIPT_DIR}"/*.R; do
            [[ -f "${f}" ]] || continue
            Rscript -e "parse(file='${f}')" >/dev/null
            log "R parse OK: ${f}"
        done
    else
        warn "Rscript not found; R parse checks skipped."
    fi
}


check_tool_current() {
    local t="$1"
    if command -v "${t}" >/dev/null 2>&1; then
        log "Executable found in current environment: ${t}"
    else
        warn "Executable not found in current environment: ${t}"
    fi
}


check_tool_env() {
    local env_name="$1"
    local tool="$2"
    local found_path=""

    if [[ "${CONDA_DEFAULT_ENV:-}" == "${env_name}" ]]; then
        found_path=$(command -v "${tool}" 2>/dev/null || true)
        if [[ -n "${found_path}" ]]; then
            log "Executable found in active ${env_name}: ${tool} (${found_path})"
            return 0
        fi
    fi

    found_path=$(bash -lc "source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate \"${env_name}\" && export PATH=\"\${CONDA_PREFIX}/bin:\${PATH}\" && hash -r && command -v \"${tool}\"" 2>/dev/null || true)
    if [[ -n "${found_path}" ]]; then
        log "Executable found in ${env_name}: ${tool} (${found_path})"
    else
        warn "Executable not found in ${env_name}: ${tool}"
    fi
}


check_envs_and_tools() {
    command -v conda >/dev/null 2>&1 || die "conda not found"
    conda env list | awk '{print $1}' | grep -qx "${TP_WGS_ENV}" || warn "Conda environment not listed: ${TP_WGS_ENV}"
    conda env list | awk '{print $1}' | grep -qx "${GUBBINS_ENV}" || warn "Conda environment not listed: ${GUBBINS_ENV}"
    conda env list | awk '{print $1}' | grep -qx "${ESM2_ENV}" || warn "Conda environment not listed: ${ESM2_ENV}"

    local t
    for t in python3 bash awk sed flock conda; do check_tool_current "${t}"; done
    for t in Rscript bwa samtools bcftools bgzip tabix iqtree treetime beast loganalyser treeannotator hyphy; do check_tool_env "${TP_WGS_ENV}" "${t}"; done
    check_tool_env "${GUBBINS_ENV}" run_gubbins.py
    check_tool_env "${ESM2_ENV}" python3
    check_tool_env "${PBP_BBMAP_ENV}" bbmap.sh
    check_tool_env "${PBP_BBMAP_ENV}" java
    command -v "${SNAKEMAKE_BIN}" >/dev/null 2>&1 || warn "Snakemake controller not found: ${SNAKEMAKE_BIN}"
    command -v qsub >/dev/null 2>&1 || warn "SGE qsub not found"

    local r_pkg_check
    r_pkg_check=$(bash -lc "source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate \"${TP_WGS_ENV}\" && export PATH=\"\${CONDA_PREFIX}/bin:\${PATH}\" && hash -r && Rscript -e 'pkgs <- c(\"ape\",\"phytools\",\"phangorn\",\"castor\",\"readxl\",\"tidyverse\"); print(sapply(pkgs, requireNamespace, quietly=TRUE))'" 2>&1 || true)
    log "R package check in ${TP_WGS_ENV}:"
    echo "${r_pkg_check}"
}


check_inputs() {
    cd_project
    require_file "refs/TPA_Nichols.fasta" "Nichols reference FASTA"
    require_any_file "SS14 reference FASTA" "refs/TPA_SS14.fasta" "refs/NC_021508.1.fasta"
    require_file "TPA_master_metadata.xlsx" "metadata file"
    require_file "sample_dates.tsv" "sample dates file"
    if [[ -f "config/complete_genome_manifest.tsv" ]]; then
        log "Found complete genome manifest: config/complete_genome_manifest.tsv"
    else
        warn "Complete genome manifest not found; p0 will be unavailable until config/complete_genome_manifest.tsv is created."
    fi
    local msa
    if msa=$(input_msa_path); then
        log "Found SS14-coordinate MSA: ${msa}"
    else
        warn "SS14-coordinate MSA not found. Run bash master.sh p0 to build consensus/TPA_ss14ref_msa.fasta, or set INPUT_MSA=/path/to/msa.fasta."
    fi
    log "Key input check finished."
}


check_all() {
    cd_project
    check_scripts
    check_syntax
    check_envs_and_tools
    check_inputs
    log "Check finished."
}


prepare_wgs() {
    cd_project
    command -v qsub >/dev/null 2>&1 || die "SGE qsub is required for FASTQ processing"
    command -v "${SNAKEMAKE_BIN}" >/dev/null 2>&1 || die "Snakemake executable not found: ${SNAKEMAKE_BIN}"
    conda run -n "${PBP_BBMAP_ENV}" which bbmap.sh >/dev/null 2>&1 || die "BBMap is missing from ${PBP_BBMAP_ENV}. Run bin/setup_bbmap_runtime.sh ${PBP_BBMAP_ENV}"
    local bbmap_java
    bbmap_java=$(conda run -n "${PBP_BBMAP_ENV}" java -version 2>&1 | head -n 1)
    [[ "${bbmap_java}" == *'"17.'* ]] || die "${PBP_BBMAP_ENV} must use OpenJDK 17; found: ${bbmap_java}"
    [[ "${PBP_BBMAP_STARTUP_TIMEOUT}" =~ ^[1-9][0-9]*$ ]] || die "PBP_BBMAP_STARTUP_TIMEOUT must be a positive integer"
    require_file "refs/TPA_SS14.fasta" "SS14 reference FASTA"
    require_file "refs/TPA_Nichols.fasta" "Nichols reference FASTA"
    run_cmd_env "wgs_reference_indexes" "${TP_WGS_ENV}" "for ref in refs/TPA_SS14.fasta refs/TPA_Nichols.fasta; do [[ -s \"\${ref}.fai\" ]] || samtools faidx \"\${ref}\"; [[ -s \"\${ref}.bwt\" ]] || bwa index \"\${ref}\"; dict=\"\${ref%.fasta}.dict\"; [[ -s \"\${dict}\" ]] || samtools dict -o \"\${dict}\" \"\${ref}\"; done"
}


run_snakemake_stage() {
    local stage="$1"; shift
    [[ "${BIG_MEM_JOBS}" =~ ^[1-9][0-9]*$ ]] || die "BIG_MEM_JOBS must be a positive integer: ${BIG_MEM_JOBS}"
    [[ "${RESTART_TIMES}" =~ ^[0-9]+$ ]] || die "RESTART_TIMES must be a non-negative integer: ${RESTART_TIMES}"
    local submit="${PBP_ANALYSIS_PY} ${SCRIPT_DIR}/../bin/submit_sge.py --project ${PROJECT_ROOT}"
    local targets=""
    local target
    for target in "$@"; do targets+=" '${target}'"; done
    run_cmd "wgs_${stage}" "'${SNAKEMAKE_BIN}' -s '${PROJECT_ROOT}/SC/${stage}.smk' --directory '${PROJECT_ROOT}' --jobs '${JOBS}' --local-cores 1 --resources 'big_mem_slots=${BIG_MEM_JOBS}' --restart-times '${RESTART_TIMES}' --rerun-incomplete --latency-wait 120 --cluster-sync '${submit}'${targets}"
}


run_wgs() {
    prepare_wgs
    local stage
    for stage in p0_01_ss14_wgs p0_02_nichols_wgs p0_03_ss14_align2nichols p0_04_nichols_align2ss14; do
        run_snakemake_stage "${stage}"
    done
}


ensure_rooted_tree() {
    cd_project
    if [[ -f "results/03_iqtree/tpa_ss14_rooted.treefile" ]]; then
        log "Rooted tree exists: results/03_iqtree/tpa_ss14_rooted.treefile"
        return 0
    fi
    require_file "results/03_iqtree/tpa_ss14.treefile" "unrooted IQ-TREE tree"
    run_cmd_env "p1_midpoint_root_tree" "${TP_WGS_ENV}" "Rscript - <<'RS'
suppressPackageStartupMessages({ library(ape); library(phytools) })
tree <- read.tree('results/03_iqtree/tpa_ss14.treefile')
tree_rooted <- phytools::midpoint.root(tree)
write.tree(tree_rooted, 'results/03_iqtree/tpa_ss14_rooted.treefile')
RS"
    require_output "results/03_iqtree/tpa_ss14_rooted.treefile" "rooted ML tree"
}


run_p0() {
    cd_project
    run_cmd_env "p0_01_complete_genomes_to_ss14_snp_consensus" "${TP_WGS_ENV}" "PROJECT_ROOT='${PROJECT_ROOT}' THREADS='${THREADS}' FORCE='${FORCE}' bash '${SCRIPT_DIR}/p0_01_complete_genomes_to_ss14_snp_consensus.sh'"
    require_output "consensus/TPA_SS14_snp_only_consensus/complete_genome_snp_consensus_qc.tsv" "complete genome SNP-consensus QC table"
    require_output "consensus/TPA_SS14_snp_only_consensus/complete_genome_snp_consensus_for_msa.fasta" "complete genome FASTA for MSA construction"
    require_output "consensus/TPA_ss14ref_msa_sources.tsv" "SS14-coordinate MSA source table"
    require_output "consensus/TPA_ss14ref_msa.fasta" "SS14-coordinate MSA"
    cat <<MSG

p0 finished.
Review:
  ${PROJECT_ROOT}/consensus/TPA_SS14_snp_only_consensus/complete_genome_snp_consensus_qc.tsv
  ${PROJECT_ROOT}/consensus/TPA_ss14ref_msa_sources.tsv

MSA written to:
  ${PROJECT_ROOT}/consensus/TPA_ss14ref_msa.fasta

Then run:
  bash master.sh p1

MSG
}


run_p0_02() {
    cd_project
    run_cmd_env "p0_02_complete_genomes_to_nichols_ann_tsv" "${TP_WGS_ENV}" "PROJECT_ROOT='${PROJECT_ROOT}' THREADS='${THREADS}' SNPEFF_ENV='${SNPEFF_ENV}' FORCE='${FORCE}' bash '${SCRIPT_DIR}/p0_02_complete_genomes_to_nichols_ann_tsv.sh'"
    require_output "results/04_complete_genome_variant_calls/complete_genome_nichols_annotation_qc.tsv" "complete-genome Nichols annotation QC table"
    require_output "results/04_complete_genome_variant_calls/ann_tsv" "complete-genome Nichols ann.tsv directory"
}


ensure_complete_genome_ann_for_p2() {
    cd_project
    local ann_dir="results/04_complete_genome_variant_calls/ann_tsv"
    if [[ "${FORCE_COMPLETE_ANN}" != "1" ]] && compgen -G "${ann_dir}/*.ann.tsv" >/dev/null; then
        log "Complete-genome Nichols ann.tsv files already exist; skipping p0_02: ${ann_dir}"
        log "To regenerate them, run: FORCE_COMPLETE_ANN=1 bash master.sh all_until_beauti"
        return 0
    fi
    if [[ "${FORCE_COMPLETE_ANN}" == "1" ]]; then
        log "FORCE_COMPLETE_ANN=1; regenerating complete-genome Nichols ann.tsv files before p2."
    else
        log "Complete-genome Nichols ann.tsv files not found; running p0_02 before p2."
    fi
    run_p0_02
}


ensure_msa_for_pipeline() {
    cd_project
    local msa="consensus/TPA_ss14ref_msa.fasta"
    if [[ -s "${msa}" && "${FORCE_MSA}" != "1" ]]; then
        log "SS14-coordinate MSA already exists; skipping p0: ${msa}"
        log "To rebuild the MSA, run: FORCE_MSA=1 bash master.sh all_until_beauti"
        return 0
    fi
    if [[ -s "${msa}" && "${FORCE_MSA}" == "1" ]]; then
        log "FORCE_MSA=1; rebuilding SS14-coordinate MSA before p1."
    else
        log "SS14-coordinate MSA not found; running p0 before p1."
    fi
    run_p0
}


run_p1() {
    cd_project
    run_cmd_env "p1_01_build_tree" "${TP_WGS_ENV}" "PROJECT_ROOT='${PROJECT_ROOT}' INPUT_MSA='${INPUT_MSA}' THREADS='${THREADS}' TP_WGS_ENV='${TP_WGS_ENV}' GUBBINS_ENV='${GUBBINS_ENV}' FORCE='${FORCE}' bash '${SCRIPT_DIR}/p1_01_build_tree.sh'"
    require_output "results/02_gubbins/tpa_ss14.filtered_polymorphic_sites.fasta" "Gubbins filtered SNP alignment"
    require_output "results/03_iqtree/tpa_ss14.treefile" "IQ-TREE ML tree"
    ensure_rooted_tree
}


prepare_p2_ann_inputs() {
    cd_project
    local combined_dir="results/04_complete_genome_variant_calls/wgs_ann_tsv"
    local manifest="results/04_complete_genome_variant_calls/wgs_ann_tsv_manifest.tsv"
    local ss14_dir="VCF/ss14_align2_nichols"
    local nichols_dir="VCF/nichols"
    local complete_dir="results/04_complete_genome_variant_calls/ann_tsv"

    [[ -d "${ss14_dir}" ]] || die "Missing SS14-to-Nichols ann.tsv directory: ${ss14_dir}"
    [[ -d "${nichols_dir}" ]] || die "Missing Nichols ann.tsv directory: ${nichols_dir}"
    [[ -d "${complete_dir}" ]] || die "Missing complete-genome Nichols ann.tsv directory: ${complete_dir}. Run: bash master.sh p0_02"

    mkdir -p "${combined_dir}"
    rm -f "${combined_dir}"/*.ann.tsv
    printf "sample_id	source_group	source_file	linked_file	status
" > "${manifest}"

    local d group f base target sample_id
    for spec in "${ss14_dir}:ss14_align2_nichols" "${nichols_dir}:nichols"; do
        d="${spec%%:*}"
        group="${spec##*:}"
        while IFS= read -r f; do
            [[ -n "${f}" ]] || continue
            base=$(basename "${f}")
            sample_id="${base%.ann.tsv}"
            target="${combined_dir}/${base}"
            if [[ -e "${target}" ]]; then
                printf "%s	%s	%s	%s	DUPLICATE_SKIPPED
" "${sample_id}" "${group}" "${f}" "${target}" >> "${manifest}"
                warn "Duplicate WGS ann.tsv basename skipped: ${base} from ${f}"
            else
                ln -s "$(realpath "${f}")" "${target}"
                printf "%s	%s	%s	%s	LINKED
" "${sample_id}" "${group}" "${f}" "${target}" >> "${manifest}"
            fi
        done < <(find "${d}" -maxdepth 1 -type f -name "*.ann.tsv" | sort)
    done

    local n_ss14 n_nichols n_combined n_complete
    n_ss14=$(find "${ss14_dir}" -maxdepth 1 -type f -name "*.ann.tsv" | wc -l)
    n_nichols=$(find "${nichols_dir}" -maxdepth 1 -type f -name "*.ann.tsv" | wc -l)
    n_combined=$(find "${combined_dir}" -maxdepth 1 -type l -name "*.ann.tsv" | wc -l)
    n_complete=$(find "${complete_dir}" -maxdepth 1 -type f -name "*.ann.tsv" | wc -l)

    log "Prepared p2 WGS ann.tsv directory: ${combined_dir}"
    log "  SS14-to-Nichols ann.tsv: ${n_ss14}"
    log "  Nichols ann.tsv: ${n_nichols}"
    log "  Combined WGS ann.tsv links: ${n_combined}"
    log "  Complete-genome ann.tsv: ${n_complete}"
    log "  Manifest: ${manifest}"

    [[ "${n_combined}" -gt 0 ]] || die "No WGS ann.tsv links were created in ${combined_dir}"
    [[ "${n_complete}" -gt 0 ]] || die "No complete-genome ann.tsv files found in ${complete_dir}. Run: bash master.sh p0_02"
}


run_p2() {
    cd_project
    run_cmd_env "p2_01_extract_codon_corrected_variants" "${TP_WGS_ENV}" "python3 '${SCRIPT_DIR}/p2_01_extract_codon_corrected_variants.py' --sample_ann_dir 'VCF/ss14_align2_nichols' --sample_ann_dir 'VCF/nichols' --complete_genome_ann_dir 'results/04_complete_genome_variant_calls/ann_tsv' --outdir 'results/05_extracted_variants' --ref_fasta 'refs/TPA_Nichols.fasta'"
    require_output "results/05_extracted_variants/sample_locus_binary_matrix.tsv" "codon-corrected binary matrix"
    run_cmd "p2_build_nichols_bam_paths" "'${PBP_ANALYSIS_PY}' '${SCRIPT_DIR}/../bin/build_nichols_bam_paths.py' --project '${PROJECT_ROOT}' --output 'results/05_extracted_variants/nichols_bam_paths.txt'"
    run_cmd "p2_02_audit_codon_callability" "'${PBP_ANALYSIS_PY}' '${SCRIPT_DIR}/p2_02_audit_codon_callability.py' --p2-script '${SCRIPT_DIR}/p2_01_extract_codon_corrected_variants.py' --matrix 'results/05_extracted_variants/sample_locus_binary_matrix.tsv' --bam-paths 'results/05_extracted_variants/nichols_bam_paths.txt' --ref-fasta 'refs/TPA_Nichols.fasta' --exclude-samples 'remove_samples.txt' --outdir 'results/05_extracted_variants/callability'"
    require_output "${PBP_MATRIX_FILE}" "audited 0/1/NA matrix"
    require_output "results/05_extracted_variants/callability/callability_summary.json" "codon audit completion summary"
}


run_p3() {
    cd_project
    run_cmd_env "p3_r_environment_check" "${TP_WGS_ENV}" "echo '[R_CHECK] CONDA_PREFIX='\"\${CONDA_PREFIX}\"; echo '[R_CHECK] which Rscript='\"\$(which Rscript)\"; Rscript -e 'cat(\"Rscript:\", Sys.which(\"Rscript\"), \"\n\"); print(.libPaths()); pkgs <- c(\"ape\",\"phytools\",\"phangorn\",\"castor\",\"readxl\"); print(sapply(pkgs, requireNamespace, quietly=TRUE))'"
    run_cmd_env "p3_01_asr_simmap_binary_loci" "${TP_WGS_ENV}" "Rscript '${SCRIPT_DIR}/p3_01_asr_simmap_binary_loci.R'"
    require_output "results/06_asr/asr_evolutionary_rates.tsv" "ASR/SIMMAP rates"
    run_cmd_env "p3_02_homoplasy_ci_ri" "${TP_WGS_ENV}" "Rscript '${SCRIPT_DIR}/p3_02_homoplasy_ci_ri.R'"
    require_output "results/07_homoplasy/homoplasy_ci_ri.tsv" "homoplasy CI/RI table"
    run_cmd_env "p3_03_er_ard_model_sensitivity" "${TP_WGS_ENV}" "Rscript '${SCRIPT_DIR}/p3_03_er_ard_model_sensitivity.R'"
    require_output "results/08_model_sensitivity/er_ard_model_sensitivity.tsv" "ER/ARD sensitivity table"
}


run_selection() {
    cd_project
    local masked
    masked=$(masked_alignment_path) || die "Masked alignment not found."
    run_cmd_env "p3_05_hyphy_selection_tests" "${TP_WGS_ENV}" "python3 '${SCRIPT_DIR}/p3_05_hyphy_selection_tests.py' --alignment '${masked}' --tree 'results/03_iqtree/tpa_ss14_rooted.treefile' --outdir 'results/10_selection'"
    require_output "results/10_selection" "HyPhy selection output directory"
}


run_temporal() {
    cd_project
    run_cmd_env "p4_01_temporal_signal_screening" "${TP_WGS_ENV}" "Rscript '${SCRIPT_DIR}/p4_01_temporal_signal_screening.R'"
    require_output "results/11_temporal_signal/temporal_signal_summary.tsv" "temporal signal summary"
    require_output "results/11_temporal_signal/recommended_temporal_subset.tsv" "recommended temporal subset"
    require_output "results/11_temporal_signal/recommended_temporal_tree.nwk" "recommended temporal tree"
}


run_treetime() {
    cd_project
    run_cmd_env "p4_02_treetime_clock_screening" "${TP_WGS_ENV}" "TREETIME_MODE='${TREETIME_MODE}' RUN_RELAXED='${RUN_RELAXED_TREETIME}' bash '${SCRIPT_DIR}/p4_02_treetime_clock_screening.sh'"
    require_output "results/12_treetime/treetime_subset_summary.tsv" "TreeTime subset summary"
    require_output "results/12_treetime/selected_subset/treetime_input_samples.tsv" "selected TreeTime sample table"
}


run_beast_prepare() {
    cd_project
    run_cmd_env "p4_03_prepare_beast2_inputs" "${TP_WGS_ENV}" "bash '${SCRIPT_DIR}/p4_03_prepare_beast2_inputs.sh'"
    require_output "results/13_beast2_selected_subset/input/final_samples.tsv" "final BEAST2 sample table"
    require_output "results/13_beast2_selected_subset/input/beast2_subset_snps_variable_only.fasta" "BEAST2 variable SNP alignment"
    require_output "results/13_beast2_selected_subset/input/beast2_input_manifest.tsv" "BEAST2 input manifest"
    require_output "results/13_beast2_selected_subset/input/beauti2_configuration_guide.txt" "BEAUti2 configuration guide"
    if ! has_beast_xml; then
        cat <<MSG

BEAUti2 manual step required.
Input files are ready in:
  ${PROJECT_ROOT}/results/13_beast2_selected_subset/input/

Create at least the strict-clock XML in BEAUti2 and save it as:
  tpa_ss14_selected_subset_n<final_n>_strict_<sites>sites_ascertained.xml

Then run:
  bash master.sh beast_run
  bash master.sh all_after_beauti

MSG
    fi
}


has_beast_xml() {
    cd "${PROJECT_ROOT}"
    compgen -G "results/13_beast2_selected_subset/input/tpa_ss14_selected_subset_n*_strict_*sites_ascertained.xml" >/dev/null
}


run_beast_run() {
    cd_project
    has_beast_xml || die "Strict-clock BEAST2 XML not found. Run beast_prepare, generate XML in BEAUti2, then run beast_run."
    run_cmd_env "p4_04_run_beast2_and_summarize" "${TP_WGS_ENV}" "bash '${SCRIPT_DIR}/p4_04_run_beast2_and_summarize.sh'"
    require_output "results/13_beast2_selected_subset/run/strict/mcc_strict.tree" "strict-clock MCC tree"
}


run_pbp_temporal() {
    cd_project
    local config="${PBP_CLOCK_HISTORY_CONFIG:-${PROJECT_ROOT}/config/clock_histories.json}"
    [[ -s "${config}" ]] || die "Provide the strict and relaxed clock input paths in ${config}, or set PBP_CLOCK_HISTORY_CONFIG. Use analysis/clock_histories/config.example.json as the template."
    local output="${PBP_CLOCK_HISTORY_OUTDIR:-${PROJECT_ROOT}/results/clock_histories/run_$(date +%Y%m%d_%H%M%S)_$$}"
    [[ ! -e "${output}" ]] || die "Use a new clock-history output directory: ${output}"
    run_cmd_env "clock_histories" "${TP_WGS_ENV}" "python3 '${PBP_PUBLICATION_PACKAGE}/analysis/clock_histories/run.py' --root '${PROJECT_ROOT}' --config '${config}' --outdir '${output}'"
    require_output "${output}/clock_sensitivity/origin_clock_comparison.tsv" "strict/relaxed PBP event comparison"
    require_output "${output}/RUN_STATUS.json" "clock-history completion status"
}


run_esm2() {
    cd_project
    [[ -n "${ESM2_MODEL_PATH}" ]] || die "Set ESM2_MODEL_PATH to the local esm2_t36_3B_UR50D model directory."
    [[ "$(hostname -s)" == "${ESM2_ALLOWED_HOST}" ]] || die "Run this stage in an allocated GPU job on ${ESM2_ALLOWED_HOST}."
    local output="${PROJECT_ROOT}/results/esm2_recomputed/run_$(date +%Y%m%d_%H%M%S)_$$"
    run_cmd_env "esm2_four_proteins" "${ESM2_ENV}" "python3 '${PBP_PUBLICATION_PACKAGE}/analysis/esm2/esm2_stage.py' recompute --project '${PROJECT_ROOT}' --model '${ESM2_MODEL_PATH}' --out '${output}'"
}

run_structural() {
    run_figures
}

run_figures() {
    cd_project
    require_script "${PBP_PUBLICATION_PACKAGE}/master.sh"
    run_cmd_env "publication_outputs" "${PBP_FIGURE_ENV:-pbp_figures_release}" "bash '${PBP_PUBLICATION_PACKAGE}/master.sh' publish --project '${PROJECT_ROOT}'"
}

run_all_until_beauti() {
    init_dirs
    run_wgs
    ensure_msa_for_pipeline
    run_p1
    ensure_complete_genome_ann_for_p2
    run_p2
    run_p3
    run_selection
    run_temporal
    run_treetime
    run_beast_prepare
}


run_all_after_beauti() {
    run_beast_run
    run_pbp_temporal
    run_figures
}


run_all() {
    run_all_until_beauti
    if has_beast_xml; then
        run_all_after_beauti
    else
        log "BEAST2 XML is not present. Pipeline stopped after beast_prepare."
    fi
}


run_all_from_p0() {
    init_dirs
    run_p0
    log "Pipeline stopped after p0 for MSA source-table review. Continue with: bash master.sh p1"
}

cmd=${1:-""}
case "${cmd}" in
    check) check_all ;;
    init) init_dirs ;;
    wgs) run_wgs ;;
    p0) run_p0 ;;
    p0_02) run_p0_02 ;;
    p1) run_p1 ;;
    p2) run_p2 ;;
    p3) run_p3 ;;
    selection) run_selection ;;
    temporal) run_temporal ;;
    treetime) run_treetime ;;
    beast_prepare) run_beast_prepare ;;
    beast_run) run_beast_run ;;
    pbp_temporal) run_pbp_temporal ;;
    esm2) run_esm2 ;;
    structural|figures) run_figures ;;
    all_until_beauti) run_all_until_beauti ;;
    all_after_beauti) run_all_after_beauti ;;
    all) run_all ;;
    all_from_p0) run_all_from_p0 ;;
    -h|--help|help|"") usage ;;
    *) usage; die "Unknown command: ${cmd}" ;;
esac
