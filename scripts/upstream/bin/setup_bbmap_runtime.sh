#!/usr/bin/env bash
set -euo pipefail

ENV_NAME=${1:-pbp_bbmap_runtime}

[[ "${ENV_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "Environment name must contain only letters, numbers, dot, underscore or dash: ${ENV_NAME}" >&2
    exit 2
}

command -v conda >/dev/null 2>&1 || {
    echo "conda is required" >&2
    exit 1
}

if conda run -n "${ENV_NAME}" true >/dev/null 2>&1; then
    conda install -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
        --strict-channel-priority bbmap=39.81 openjdk=17 pigz
else
    conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
        --strict-channel-priority bbmap=39.81 openjdk=17 pigz
fi

conda run -n "${ENV_NAME}" which bbmap.sh
conda run -n "${ENV_NAME}" which java
conda run -n "${ENV_NAME}" java -version
