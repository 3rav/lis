#!/usr/bin/env bash
# Public deterministic 8x8 synthetic mode discovery and audit.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
cat > "$tmp/A.mtx" <<'A'
%%MatrixMarket matrix coordinate real general
8 8 8
1 1 10
2 2 12
3 3 0.01
4 4 0.02
5 5 0.03
6 6 4
7 7 5
8 8 6
A
cat > "$tmp/base.mtx" <<'B'
%%MatrixMarket matrix array real general
8 2
1
0
0
0
0
0
0
0
0
1
0
0
0
0
0
0
B
python3 "$here/mfamg_mode_auto_candidates.py" "$tmp/A.mtx" "$tmp/base.mtx" \
  --count 2 --steps 180 --seed 51 --output "$tmp/auto.mtx" > "$tmp/report.txt"
grep -q '^AUTO_MODE_EXTRA_COUNT=2$' "$tmp/report.txt"
grep -q '^8 4$' "$tmp/auto.mtx"
python3 "$here/mfamg_mode_audit.py" "$tmp/A.mtx" "$tmp/base.mtx" "$tmp/auto.mtx" \
  --csv-out "$tmp/audit.csv" > "$tmp/audit.txt"
grep -q '^MODE_AUDIT_EXTRA_NEW_RANK=2$' "$tmp/audit.txt"
if python3 "$here/mfamg_mode_auto_candidates.py" "$tmp/A.mtx" "$tmp/base.mtx" \
  --count 2 --steps 1 --seed 51 --output "$tmp/auto.mtx" > /dev/null 2>&1; then
  echo 'ERROR: automatic mode output overwritten' >&2
  exit 1
fi
echo AUTO_MODE_CANDIDATES_SMOKE_OK
