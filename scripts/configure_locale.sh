#!/usr/bin/env bash
# Use an installed UTF-8 locale when available; keep numeric and sort conventions fixed.
unset LC_ALL
pbp_utf8_locale="$(LC_ALL=C locale -a 2>/dev/null | LC_ALL=C awk 'tolower($0) ~ /^(c|en_us)[.](utf8|utf-8)$/ && !found {print; found=1}' || true)"
export LANG="${pbp_utf8_locale:-C}" LC_CTYPE="${pbp_utf8_locale:-C}"
export LC_NUMERIC=C LC_COLLATE=C LC_TIME=C LC_MESSAGES=C
unset pbp_utf8_locale
