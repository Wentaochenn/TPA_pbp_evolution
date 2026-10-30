#!/usr/bin/env bash
# p4_02_treetime_clock_screening.sh
#
# Run TreeTime clock screening for temporal-signal candidate sets.
#
# Inputs:
#   results/11_temporal_signal/recommended_temporal_tree.nwk
#   results/11_temporal_signal/recommended_temporal_subset.tsv
#   results/11_temporal_signal/subsampled_tree_<profile>.nwk
#   results/11_temporal_signal/subsampled_phylo_informed_<profile>.tsv
#
# Outputs:
#   results/12_treetime/<subset_id>/strict_clock/
#   results/12_treetime/selected_subset/
#   results/12_treetime/treetime_subset_summary.tsv
#
# Usage:
#   bash SC/p4_02_treetime_clock_screening.sh
#   TREETIME_MODE=all_subsets bash SC/p4_02_treetime_clock_screening.sh

set -euo pipefail

TEMPORAL_DIR=${TEMPORAL_DIR:-"results/11_temporal_signal"}
OUTDIR=${OUTDIR:-"results/12_treetime"}
TREETIME_MODE=${TREETIME_MODE:-"recommended"}
SEQ_LENGTH=${SEQ_LENGTH:-1083903}
CLOCK_FILTER=${CLOCK_FILTER:-3}
RUN_RELAXED=${RUN_RELAXED:-false}
SELECTED_SUBSET_ID=${SELECTED_SUBSET_ID:-"selected_subset"}

mkdir -p "${OUTDIR}" logs
SUMMARY="${OUTDIR}/treetime_subset_summary.tsv"
LOG_FILE="${OUTDIR}/p4_02_treetime_clock_screening.log"

die() { echo "[ERROR] $*" | tee -a "${LOG_FILE}" >&2; exit 1; }
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"; }

command -v treetime >/dev/null 2>&1 || die "TreeTime executable not found."

write_dates_csv() {
    local tsv="$1"
    local csv="$2"
    echo "name,date" > "${csv}"
    awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++){if($i=="sample") s=i; if($i=="year") y=i}} NR>1 && s>0 && y>0 {print $s","$y}' "${tsv}" >> "${csv}"
}

run_subset() {
    local subset_id="$1"
    local tree="$2"
    local samples_tsv="$3"

    [[ -f "${tree}" ]] || { log "[SKIP] ${subset_id}: tree not found (${tree})"; return 0; }
    [[ -f "${samples_tsv}" ]] || { log "[SKIP] ${subset_id}: sample table not found (${samples_tsv})"; return 0; }

    local subset_dir="${OUTDIR}/${subset_id}"
    local strict_dir="${subset_dir}/strict_clock"
    local dates_csv="${subset_dir}/treetime_dates.csv"
    mkdir -p "${subset_dir}"

    cp "${tree}" "${subset_dir}/treetime_input_tree.nwk"
    cp "${samples_tsv}" "${subset_dir}/treetime_input_samples.tsv"
    write_dates_csv "${samples_tsv}" "${dates_csv}"

    local input_n
    input_n=$(awk 'NR>1' "${samples_tsv}" | wc -l | tr -d ' ')

    log "Running TreeTime subset=${subset_id}; samples=${input_n}"
    treetime --tree "${tree}" --dates "${dates_csv}" --sequence-length "${SEQ_LENGTH}" --outdir "${strict_dir}" --clock-filter "${CLOCK_FILTER}" --reroot best --coalescent opt --verbose 2 > "${subset_dir}/treetime_stdout.log" 2> "${subset_dir}/treetime_stderr.log" || true

    local outliers=0
    if [[ -f "${strict_dir}/outliers.tsv" ]]; then
        outliers=$(awk 'NR>1' "${strict_dir}/outliers.tsv" | wc -l | tr -d ' ')
    fi
    local final_n=$(( input_n - outliers ))
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${subset_id}" "${tree}" "${samples_tsv}" "${input_n}" "${outliers}" "${final_n}" >> "${SUMMARY}"

    if [[ "${RUN_RELAXED}" == "true" ]]; then
        treetime --tree "${tree}" --dates "${dates_csv}" --sequence-length "${SEQ_LENGTH}" --outdir "${subset_dir}/relaxed_clock" --clock-filter "${CLOCK_FILTER}" --reroot best --coalescent skyline --relax 0.1 0 --verbose 2 >> "${subset_dir}/treetime_stdout.log" 2>> "${subset_dir}/treetime_stderr.log" || true
    fi
}

printf "subset_id\ttree_file\tsample_file\tinput_samples\toutliers_removed\tfinal_samples_after_outlier_filter\n" > "${SUMMARY}"

if [[ "${TREETIME_MODE}" == "recommended" ]]; then
    run_subset "recommended_temporal_subset" "${TEMPORAL_DIR}/recommended_temporal_tree.nwk" "${TEMPORAL_DIR}/recommended_temporal_subset.tsv"
else
    for tree in "${TEMPORAL_DIR}"/subsampled_tree_*.nwk; do
        [[ -e "${tree}" ]] || continue
        base=$(basename "${tree}" .nwk)
        profile=${base#subsampled_tree_}
        run_subset "${profile}" "${tree}" "${TEMPORAL_DIR}/subsampled_phylo_informed_${profile}.tsv"
    done
    run_subset "recommended_temporal_subset" "${TEMPORAL_DIR}/recommended_temporal_tree.nwk" "${TEMPORAL_DIR}/recommended_temporal_subset.tsv"
fi

selected_source="recommended_temporal_subset"
if [[ "${TREETIME_MODE}" == "all_subsets" ]]; then
    selected_source=$(awk -F'\t' 'NR>1 {print $1"\t"$6}' "${SUMMARY}" | sort -k2,2nr | head -1 | cut -f1)
fi

[[ -n "${selected_source}" ]] || die "No TreeTime subset was selected."
rm -rf "${OUTDIR}/${SELECTED_SUBSET_ID}"
cp -r "${OUTDIR}/${selected_source}" "${OUTDIR}/${SELECTED_SUBSET_ID}"
echo "${selected_source}" > "${OUTDIR}/selected_subset_id.txt"

log "TreeTime screening finished."
log "Summary: ${SUMMARY}"
log "Selected subset: ${selected_source} -> ${OUTDIR}/${SELECTED_SUBSET_ID}"
