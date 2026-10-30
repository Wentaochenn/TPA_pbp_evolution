#!/usr/bin/env bash
set -euo pipefail
: "${ENSEMBLE_ROOT:?}" "${ENSEMBLE_PACKAGE:?}" "${ENSEMBLE_RUN:?}" "${ENSEMBLE_ACTION:?}"
mkdir -p "$ENSEMBLE_RUN/scheduler"
cd "$ENSEMBLE_ROOT"
case "$ENSEMBLE_ACTION" in
  prepare) LOG="$ENSEMBLE_RUN/prepare.console.log" ;;
  fit) : "${SGE_TASK_ID:?}"; LOG="$ENSEMBLE_RUN/task_${SGE_TASK_ID}.console.log" ;;
  collect) LOG="$ENSEMBLE_RUN/collect.console.log" ;;
  *) exit 2 ;;
esac
exec >> "$LOG" 2>&1
finish() { code=$?; trap - EXIT; printf '[END] %s action=%s exit_code=%s\n' "$(date -Is)" "$ENSEMBLE_ACTION" "$code"; exit "$code"; }
trap finish EXIT
printf '[START] %s job=%s task=%s action=%s\n' "$(date -Is)" "${JOB_ID:-manual}" "${SGE_TASK_ID:-none}" "$ENSEMBLE_ACTION"
if [ -n "${PBP_CONDA_SH:-}" ]; then
  source "$PBP_CONDA_SH"
elif command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  printf 'Conda not available; set PBP_CONDA_SH or activate the analysis environment.\n'; exit 2
fi
conda activate "${ENSEMBLE_ENV:-tp_wgs}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS="${NSLOTS:-4}"
python3 - "$ENSEMBLE_RUN/config.json" "${NSLOTS:-4}" <<'PY'
import json,sys
c=json.load(open(sys.argv[1]));slots=int(sys.argv[2])
if slots<c['threads']:raise SystemExit('Allocated slots fewer than configured HyPhy threads')
PY
args=("$ENSEMBLE_ACTION" --root "$ENSEMBLE_ROOT" --run "$ENSEMBLE_RUN")
if [ "$ENSEMBLE_ACTION" = fit ]; then args+=(--task "$SGE_TASK_ID"); fi
python3 -u "$ENSEMBLE_PACKAGE/ensemble.py" "${args[@]}"
