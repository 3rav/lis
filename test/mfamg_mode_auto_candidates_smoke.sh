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
# Signed curvature is tested on a known symmetric-indefinite operator.
cat > "$tmp/A_indef.mtx" <<'N'
%%MatrixMarket matrix coordinate real general
8 8 8
1 1 10
2 2 12
3 3 -3
4 4 -1
5 5 0.03
6 6 4
7 7 5
8 8 6
N
python3 "$here/mfamg_mode_auto_candidates.py" "$tmp/A_indef.mtx" "$tmp/base.mtx" \
  --target negative-curvature --count 2 --steps 180 --seed 51 \
  --output "$tmp/negative.mtx" > "$tmp/negative_report.txt"
grep -q '^AUTO_MODE_TARGET=negative-curvature' "$tmp/negative_report.txt"
grep -q '^AUTO_MODE_01_RAYLEIGH_AFTER=-' "$tmp/negative_report.txt"
grep -q '^AUTO_MODE_02_RAYLEIGH_AFTER=-' "$tmp/negative_report.txt"
python3 "$here/mfamg_mode_audit.py" "$tmp/A_indef.mtx" "$tmp/base.mtx" "$tmp/negative.mtx" \
  --csv-out "$tmp/negative_audit.csv" > "$tmp/negative_audit.txt"
grep -q '^MODE_AUDIT_EXTRA_NEW_RANK=2' "$tmp/negative_audit.txt"
grep -q '^MODE_AUDIT_EXTRA_SYMPART_NEGATIVE_COUNT=2' "$tmp/negative_audit.txt"
python3 "$here/mfamg_mode_auto_candidates.py" "$tmp/A_indef.mtx" "$tmp/base.mtx" \
  --target negative-curvature --count 2 --steps 180 --seed 51 \
  --output "$tmp/negative_copy.mtx" > /dev/null
cmp "$tmp/negative.mtx" "$tmp/negative_copy.mtx"
if python3 "$here/mfamg_mode_auto_candidates.py" "$tmp/A.mtx" "$tmp/base.mtx" \
  --target negative-curvature --count 1 --steps 180 --seed 51 \
  --output "$tmp/not_negative.mtx" > "$tmp/positive_attempt.txt" 2>&1; then
  echo 'ERROR: negative-curvature filter accepted a positive definite matrix' >&2
  exit 1
fi
grep -q 'has no verified negative curvature' "$tmp/positive_attempt.txt"
test ! -e "$tmp/not_negative.mtx"
echo AUTO_MODE_CANDIDATES_SMOKE_OK
