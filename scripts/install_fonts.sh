#!/usr/bin/env bash
set -euo pipefail
PACKAGE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ $# -gt 1 ]]; then
    echo "Usage: bash scripts/install_fonts.sh [/path/to/Arial]" >&2
    exit 2
fi
PBP_FONT_CHECK="$(mktemp -d)"
trap 'rm -rf -- "$PBP_FONT_CHECK"' EXIT
if [[ $# -eq 1 ]]; then
    python3 "$PACKAGE/figure_fonts.py" check --font-dir "$1" --out "$PBP_FONT_CHECK"
    PBP_FONT_DEST="${XDG_DATA_HOME:-$HOME/.local/share}/fonts/TPA_pbp_evolution/Arial"
    python3 - "$PBP_FONT_CHECK/FONT_PREFLIGHT.json" "$PBP_FONT_DEST" <<'PY'
import hashlib
import json
from pathlib import Path
import shutil
import sys

record = json.loads(Path(sys.argv[1]).read_text())
if record['state'] != 'PASS':
    raise SystemExit('The supplied fonts must contain all four genuine Arial faces')
destination = Path(sys.argv[2]); destination.mkdir(parents=True, exist_ok=True)
for source in {Path(face['file']) for face in record['faces'].values()}:
    target = destination/source.name
    if target.exists() and target.resolve() != source.resolve():
        if hashlib.sha256(target.read_bytes()).digest() != hashlib.sha256(source.read_bytes()).digest():
            raise SystemExit('A different font file already exists: '+str(target))
    elif target.resolve() != source.resolve():
        shutil.copy2(source, target)
PY
    fc-cache -f "$PBP_FONT_DEST"
    python3 "$PACKAGE/figure_fonts.py" check --out "$PBP_FONT_CHECK/installed"
else
    python3 "$PACKAGE/figure_fonts.py" check --out "$PBP_FONT_CHECK"
fi
echo "Arial regular, bold, italic and bold-italic fonts are available."
