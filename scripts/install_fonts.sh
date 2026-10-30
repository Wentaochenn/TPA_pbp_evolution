#!/usr/bin/env bash
set -euo pipefail
PACKAGE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
DEST="${XDG_DATA_HOME:-$HOME/.local/share}/fonts/pbp-cluster-figures-20261027"
mkdir -p "$DEST"
cp "$PACKAGE"/assets/fonts/NimbusSans-*.otf "$DEST/"
fc-cache -f "$DEST"
