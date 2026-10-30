#!/usr/bin/env bash
set -euo pipefail

ACTION=${1:-}
SOURCE=${2:-}
CACHE_KEY=${3:-}
RUN_USER=${USER:-$(id -un)}
CACHE_ROOT=${PBP_BBMAP_CACHE_ROOT:-/dev/shm/${RUN_USER}/tp_wgs_bbmap_cache}
MIN_FREE_MB=${PBP_BBMAP_CACHE_MIN_FREE_MB:-20480}

die() {
    echo "[BBMAP CACHE ERROR] $*" >&2
    exit 1
}

validate_root() {
    case "${CACHE_ROOT}" in
        /dev/shm/${RUN_USER}/tp_wgs_bbmap_cache|/dev/shm/${RUN_USER}/tp_wgs_bbmap_cache/*) ;;
        *) die "Cache root must remain under /dev/shm/${RUN_USER}/tp_wgs_bbmap_cache: ${CACHE_ROOT}" ;;
    esac
}

validate_root

case "${ACTION}" in
    prepare)
        [[ -d "${SOURCE}" ]] || die "Index directory not found: ${SOURCE}"
        [[ "${CACHE_KEY}" =~ ^[A-Za-z0-9._-]+$ ]] || die "Invalid cache key: ${CACHE_KEY}"
        [[ "${MIN_FREE_MB}" =~ ^[0-9]+$ ]] || die "PBP_BBMAP_CACHE_MIN_FREE_MB must be a non-negative integer"
        command -v flock >/dev/null 2>&1 || die "flock is required on compute nodes"

        mkdir -p "${CACHE_ROOT}"
        DESTINATION="${CACHE_ROOT}/${CACHE_KEY}"
        LOCK_FILE="${CACHE_ROOT}/.${CACHE_KEY}.lock"
        exec 9>"${LOCK_FILE}"
        flock -x 9

        CACHE_VALID=0
        if [[ -f "${DESTINATION}/.complete" && -f "${DESTINATION}/.source_path" ]]; then
            [[ "$(cat "${DESTINATION}/.source_path")" == "${SOURCE}" ]] && CACHE_VALID=1
        fi
        if (( CACHE_VALID == 0 )); then
            SOURCE_KB=$(du -sk "${SOURCE}" | awk '{print $1}')
            AVAILABLE_KB=$(df -Pk "${CACHE_ROOT}" | awk 'NR==2 {print $4}')
            REQUIRED_KB=$((SOURCE_KB + MIN_FREE_MB * 1024))
            (( AVAILABLE_KB >= REQUIRED_KB )) || die "Insufficient local space: need ${REQUIRED_KB} KB, available ${AVAILABLE_KB} KB in ${CACHE_ROOT}"

            PARTIAL=$(mktemp -d "${CACHE_ROOT}/.${CACHE_KEY}.partial.XXXXXX")
            cleanup_partial() {
                [[ -n "${PARTIAL:-}" && -d "${PARTIAL}" ]] && rm -rf -- "${PARTIAL}"
            }
            trap cleanup_partial EXIT
            cp -a "${SOURCE}/." "${PARTIAL}/"
            printf '%s\n' "${SOURCE}" > "${PARTIAL}/.source_path"
            touch "${PARTIAL}/.complete"
            if [[ -e "${DESTINATION}" ]]; then
                INCOMPLETE="${CACHE_ROOT}/.${CACHE_KEY}.incomplete.$$"
                mv "${DESTINATION}" "${INCOMPLETE}"
                rm -rf -- "${INCOMPLETE}"
            fi
            mv "${PARTIAL}" "${DESTINATION}"
            PARTIAL=""
            trap - EXIT
            echo "[BBMAP CACHE] $(hostname): cached ${SOURCE} at ${DESTINATION}" >&2
        fi

        printf '%s\n' "${DESTINATION}"
        ;;
    status)
        hostname
        df -h "${CACHE_ROOT}" 2>/dev/null || df -h /dev/shm
        if [[ -d "${CACHE_ROOT}" ]]; then
            du -sh "${CACHE_ROOT}"
            find "${CACHE_ROOT}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
        else
            echo "cache_absent=${CACHE_ROOT}"
        fi
        ;;
    cleanup)
        if [[ -d "${CACHE_ROOT}" ]]; then
            rm -rf -- "${CACHE_ROOT}"
        fi
        echo "cache_removed=$(hostname):${CACHE_ROOT}"
        ;;
    *)
        die "Usage: $0 {prepare SOURCE KEY|status|cleanup}"
        ;;
esac
