#!/usr/bin/env bash
set -euo pipefail
PACKAGE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ENV_NAME=${PBP_FIGURE_ENV:-pbp_figures_release}
if ! command -v conda >/dev/null 2>&1; then
  if [[ -n "${PBP_CONDA_SH:-}" ]]; then
    source "$PBP_CONDA_SH"
  else
    printf "Conda is required. Activate a Conda installation first.\n" >&2
    exit 1
  fi
fi
env_list_file=$(mktemp "${TMPDIR:-/tmp}/pbp_fig_env.XXXXXX")
trap 'rm -f "$env_list_file"' EXIT
conda env list --json > "$env_list_file"
if python3 - "$env_list_file" "$ENV_NAME" <<'PY'
import json,sys
from pathlib import Path
names={Path(p).name for p in json.load(open(sys.argv[1]))['envs']}
raise SystemExit(0 if sys.argv[2] in names else 1)
PY
then
  conda env update -n "$ENV_NAME" -f "$PACKAGE/figure_environment.yml"
else
  conda env create -n "$ENV_NAME" -f "$PACKAGE/figure_environment.yml"
fi
