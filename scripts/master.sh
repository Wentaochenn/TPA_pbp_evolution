#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "$HERE/configure_locale.sh"
export PYTHONDONTWRITEBYTECODE=1
exec python3 "$HERE/cli.py" "$@"
