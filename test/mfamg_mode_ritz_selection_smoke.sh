#!/usr/bin/env bash
# Public synthetic test of projected Ritz rotation, preserving baseline.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
cat > "$tmp/A.mtx" <<'A'
%%MatrixMarket matrix coordinate real symmetric
5 5 6
1 1 -3
2 2 2
3 3 -1
4 4 4
5 5 5
2 1 0.25
A
cat > "$tmp/base.mtx" <<'B'
%%MatrixMarket matrix array real general
5 2
0
1
0
0
0
0
0
0
1
0
B
cat > "$tmp/extra.mtx" <<'E'
%%MatrixMarket matrix array real general
5 5
0
1
0
0
0
0
0
0
1
0
1
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
1
E
mkdir "$tmp/out"
python3 "$here/mfamg_mode_ritz_selection.py" \
  "$tmp/A.mtx" "$tmp/base.mtx" "$tmp/extra.mtx" \
  --output-dir "$tmp/out" > "$tmp/result.txt"
grep -q '^RITZ_SELECTION_NEGATIVE_COUNT=2$' "$tmp/result.txt"
grep -q '^RITZ_SELECTION_POSITIVE_COUNT=1$' "$tmp/result.txt"
grep -q '^RITZ_SELECTION_PROJECTED_EIGENVALUES=-3,-1,5$' "$tmp/result.txt"
grep -q '^RITZ_SELECTION_NOTE=offline_rotation_not_online_mode_discovery$' "$tmp/result.txt"
grep -q '^5 4$' "$tmp/out/modes_negative_ritz.mtx"
grep -q '^5 5$' "$tmp/out/modes_all_ritz.mtx"
python3 "$here/mfamg_mode_audit.py" \
  "$tmp/A.mtx" "$tmp/base.mtx" "$tmp/out/modes_negative_ritz.mtx" \
  --csv-out "$tmp/neg.csv" > "$tmp/neg_check.txt"
grep -q '^MODE_AUDIT_EXTRA_SYMPART_NEGATIVE_COUNT=2$' "$tmp/neg_check.txt"
grep -q '^MODE_AUDIT_EXTRA_NEW_RANK=2$' "$tmp/neg_check.txt"
python3 "$here/mfamg_mode_audit.py" \
  "$tmp/A.mtx" "$tmp/base.mtx" "$tmp/out/modes_all_ritz.mtx" \
  --csv-out "$tmp/all.csv" > "$tmp/all_check.txt"
grep -q '^MODE_AUDIT_EXTRA_NEW_RANK=3$' "$tmp/all_check.txt"
if python3 "$here/mfamg_mode_ritz_selection.py" \
  "$tmp/A.mtx" "$tmp/base.mtx" "$tmp/extra.mtx" \
  --output-dir "$tmp/out" > "$tmp/retry.out" 2> "$tmp/retry.err"; then
  echo 'ERROR: overwrite of existing modes was permitted' >&2
  exit 1
fi
grep -q '^RITZ_SELECTION_ERROR=output files exist' "$tmp/retry.err"
echo RITZ_SELECTION_SMOKE_OK
