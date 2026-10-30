#!/usr/bin/env bash
# p4_04_run_beast2_and_summarize.sh
#
# Run BEAST2 using XML files generated in BEAUti2.
#
# Inputs:
#   results/13_beast2_selected_subset/input/tpa_ss14_selected_subset_n*_strict_*sites_ascertained.xml
#   results/13_beast2_selected_subset/input/tpa_ss14_selected_subset_n*_relaxed_*sites_ascertained.xml
#
# Outputs:
#   results/13_beast2_selected_subset/run/strict/
#   results/13_beast2_selected_subset/run/relaxed/
#   results/13_beast2_selected_subset/diagnostics/
#
# Usage:
#   bash SC/p4_04_run_beast2_and_summarize.sh

set -euo pipefail

INDIR=${INDIR:-"results/13_beast2_selected_subset/input"}
OUTDIR=${OUTDIR:-"results/13_beast2_selected_subset/run"}
DIAGDIR=${DIAGDIR:-"results/13_beast2_selected_subset/diagnostics"}
MANIFEST=${MANIFEST:-"${INDIR}/beast2_input_manifest.tsv"}

BEAST=${BEAST:-beast}
LOGANALYSER=${LOGANALYSER:-loganalyser}
TREEANNOTATOR=${TREEANNOTATOR:-treeannotator}
TREEANNOTATOR_XMX=${TREEANNOTATOR_XMX:-"64g"}
THREADS=${THREADS:-8}
SEED=${SEED:-42}
BURNIN_PCT=${BURNIN_PCT:-10}

mkdir -p "${OUTDIR}" "${DIAGDIR}" logs
LOG_FILE="${OUTDIR}/p4_04_run_beast2_and_summarize.log"

die() { echo "[ERROR] $*" | tee -a "${LOG_FILE}" >&2; exit 1; }
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"; }

resolve_exec() {
    local x="$1"
    if [[ "${x}" == */* ]]; then
        [[ -x "${x}" ]] || return 1
        echo "${x}"
    else
        command -v "${x}"
    fi
}

BEAST_BIN=$(resolve_exec "${BEAST}") || die "BEAST executable not found: ${BEAST}"
LOGANALYSER_BIN=$(resolve_exec "${LOGANALYSER}") || die "LogAnalyser executable not found: ${LOGANALYSER}"
TREEANNOTATOR_BIN=$(resolve_exec "${TREEANNOTATOR}") || die "TreeAnnotator executable not found: ${TREEANNOTATOR}"

manifest_value() {
    local key="$1"
    [[ -f "${MANIFEST}" ]] || return 1
    awk -F'\t' -v k="${key}" 'NR>1 && $1==k {print $2; exit}' "${MANIFEST}"
}

XML_STRICT=${XML_STRICT:-$(manifest_value strict_xml || true)}
XML_RELAXED=${XML_RELAXED:-$(manifest_value relaxed_xml || true)}

if [[ -z "${XML_STRICT}" || ! -f "${XML_STRICT}" ]]; then
    XML_STRICT=$(ls "${INDIR}"/tpa_ss14_selected_subset_n*_strict_*sites_ascertained.xml 2>/dev/null | head -1 || true)
fi
if [[ -z "${XML_RELAXED}" || ! -f "${XML_RELAXED}" ]]; then
    XML_RELAXED=$(ls "${INDIR}"/tpa_ss14_selected_subset_n*_relaxed_*sites_ascertained.xml 2>/dev/null | head -1 || true)
fi

[[ -f "${XML_STRICT}" ]] || die "Strict-clock XML not found in ${INDIR}. Generate it in BEAUti2 and save it to the input directory."

check_xml_basic() {
    local xml="$1"
    python3 - "$xml" <<'PYEOF'
import re, sys
from pathlib import Path
xml = Path(sys.argv[1])
txt = xml.read_text()
if "beast2_subset_snps_variable_only" not in txt:
    raise SystemExit("[ERROR] XML does not reference beast2_subset_snps_variable_only")
seqs = re.findall(r'<sequence\b[^>]*\bvalue="([^"]*)"', txt)
print(f"XML: {xml}")
print(f"Embedded sequences: {len(seqs)}")
if not seqs:
    raise SystemExit("[ERROR] No embedded sequence records found in XML")
PYEOF
}

run_one() {
    local xml="$1"
    local tag="$2"
    local required="$3"
    local rundir="${OUTDIR}/${tag}"
    mkdir -p "${rundir}"

    if [[ ! -f "${xml}" ]]; then
        if [[ "${required}" == "yes" ]]; then
            die "Required XML not found: ${xml}"
        else
            log "[SKIP] ${tag}: XML not found (${xml})"
            return 0
        fi
    fi

    local xml_abs
    xml_abs=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${xml}")

    log "Running BEAST2 ${tag}: ${xml_abs}"
    (cd "${rundir}" && "${BEAST_BIN}" -threads "${THREADS}" -seed "${SEED}" -overwrite "${xml_abs}")

    local logfile treefile
    logfile=$(ls "${rundir}"/*.log 2>/dev/null | head -1 || true)
    treefile=$(ls "${rundir}"/*.trees 2>/dev/null | head -1 || true)
    [[ -n "${logfile}" && -n "${treefile}" ]] || die "${tag}: .log or .trees output not found in ${rundir}"

    "${LOGANALYSER_BIN}" -b "${BURNIN_PCT}" "${logfile}" > "${rundir}/beast_${tag}.loganalyser.txt" 2>&1 || true
    cp "${rundir}/beast_${tag}.loganalyser.txt" "${DIAGDIR}/${tag}_loganalyser.txt" || true

    log "${tag}: running TreeAnnotator (Xmx=${TREEANNOTATOR_XMX})"
    _JAVA_OPTIONS="-Xmx${TREEANNOTATOR_XMX}" "${TREEANNOTATOR_BIN}" -burnin "${BURNIN_PCT}" -height median "${treefile}" "${rundir}/mcc_${tag}.tree" > "${rundir}/treeannotator_${tag}.log" 2>&1 || true
    if [[ -s "${rundir}/mcc_${tag}.tree" ]]; then
        log "${tag} MCC tree: ${rundir}/mcc_${tag}.tree"
    else
        log "[WARNING] ${tag}: TreeAnnotator did not produce an MCC tree. Check ${rundir}/treeannotator_${tag}.log (common cause: insufficient heap memory for large .trees files -- increase TREEANNOTATOR_XMX)."
    fi
}

log "Input directory: ${INDIR}"
log "Output directory: ${OUTDIR}"
log "Strict XML: ${XML_STRICT}"
log "Relaxed XML: ${XML_RELAXED:-NA}"

check_xml_basic "${XML_STRICT}"
run_one "${XML_STRICT}" "strict" "yes"
if [[ -n "${XML_RELAXED}" ]]; then
    check_xml_basic "${XML_RELAXED}"
    run_one "${XML_RELAXED}" "relaxed" "no"
else
    log "Relaxed-clock XML not provided."
fi

{
    echo -e "field\tvalue"
    echo -e "input_dir\t${INDIR}"
    echo -e "output_dir\t${OUTDIR}"
    echo -e "diagnostics_dir\t${DIAGDIR}"
    echo -e "strict_xml\t${XML_STRICT}"
    echo -e "relaxed_xml\t${XML_RELAXED:-NA}"
    echo -e "burnin_pct\t${BURNIN_PCT}"
    echo -e "seed\t${SEED}"
} > "results/13_beast2_selected_subset/beast2_run_manifest.tsv"

log "BEAST2 run finished."
log "Run manifest: results/13_beast2_selected_subset/beast2_run_manifest.tsv"