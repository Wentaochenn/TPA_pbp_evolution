#!/usr/bin/env bash
set -euo pipefail

LOG_FILE=${1:?log file is required}
STARTUP_TIMEOUT=${2:?startup timeout is required}
shift 2

[[ "${STARTUP_TIMEOUT}" =~ ^[1-9][0-9]*$ ]] || {
    echo "Invalid BBMap startup timeout: ${STARTUP_TIMEOUT}" >&2
    exit 2
}
[[ "$#" -gt 0 ]] || {
    echo "BBMap command is required" >&2
    exit 2
}

mkdir -p "$(dirname "${LOG_FILE}")"
: > "${LOG_FILE}"

setsid "$@" > >(tee -a "${LOG_FILE}") 2> >(tee -a "${LOG_FILE}" >&2) &
COMMAND_PID=$!
DEADLINE=$(( $(date +%s) + STARTUP_TIMEOUT ))

terminate_group() {
    kill -TERM -- "-${COMMAND_PID}" 2>/dev/null || true
    sleep 2
    kill -KILL -- "-${COMMAND_PID}" 2>/dev/null || true
    wait "${COMMAND_PID}" 2>/dev/null || true
}

while kill -0 "${COMMAND_PID}" 2>/dev/null; do
    if grep -Eq 'Processing reads in|Started [0-9]+ mapping threads' "${LOG_FILE}"; then
        set +e
        wait "${COMMAND_PID}"
        STATUS=$?
        set -e
        exit "${STATUS}"
    fi
    if grep -Eq 'Exception in thread|OutOfMemoryError|Cannot read field "sites"|BBMAP CACHE ERROR' "${LOG_FILE}"; then
        echo "[BBMAP GUARD] fatal startup error detected" >&2
        terminate_group
        exit 70
    fi
    if (( $(date +%s) >= DEADLINE )); then
        echo "[BBMAP GUARD] index loading did not finish within ${STARTUP_TIMEOUT}s" >&2
        terminate_group
        exit 124
    fi
    sleep 10
done

set +e
wait "${COMMAND_PID}"
STATUS=$?
set -e
if grep -Eq 'Processing reads in|Started [0-9]+ mapping threads' "${LOG_FILE}"; then
    exit "${STATUS}"
fi
if [[ "${STATUS}" -eq 0 ]]; then
    echo "[BBMAP GUARD] BBMap exited before read processing was detected" >&2
    exit 71
fi
exit "${STATUS}"
