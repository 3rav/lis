#!/usr/bin/env bash
# Public deterministic algebra checks, no private CCX data or SciPy.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
cat > "$tmp/A.mtx" <<'A'
%%MatrixMarket matrix coordinate real symmetric
% Five-row indefinite matrix with a positive/negative coupling
5 5 6
1 1 -3
2 2 2
3 3 -1
4 4 4
5 5 5
2 1 0.25
A
cat > "$tmp/base.mtx" <<'BASE'
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
BASE
cat > "$tmp/enriched.mtx" <<'MORE'
%%MatrixMarket matrix array real general
5 4
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
MORE
python3 "$here/mfamg_mode_audit.py" "$tmp/A.mtx" "$tmp/base.mtx" \
  "$tmp/enriched.mtx" --csv-out "$tmp/report.csv" > "$tmp/report.txt"
grep -q '^MODE_AUDIT_N=5$' "$tmp/report.txt"
grep -q '^MODE_AUDIT_A_MM_SYMMETRY=symmetric$' "$tmp/report.txt"
grep -q '^MODE_AUDIT_BASE_RANK=2$' "$tmp/report.txt"
grep -q '^MODE_AUDIT_EXTRA_NEW_RANK=2$' "$tmp/report.txt"
grep -q '^MODE_AUDIT_EXTRA_MIN_DISTANCE_FROM_BASE=1$' "$tmp/report.txt"
grep -q '^MODE_AUDIT_EXTRA_SYMPART_NEGATIVE_COUNT=2$' "$tmp/report.txt"
grep -q '^MODE_AUDIT_EXTRA_SYMPART_EIGENVALUES=-3,-1$' "$tmp/report.txt"
test "$(wc -l < "$tmp/report.csv")" -eq 5
# An off-diagonal term makes e1 an imperfect eigenvector even though
# its Rayleigh quotient is -3.
awk -F, 'NR==4 && $4 == -3 && $5 > 0 && $5 < 1 {ok=1} END {exit !ok}' "$tmp/report.csv"
# Do not quietly accept changes to original baseline modes.
awk 'NR==4 {$0="0"} {print}' "$tmp/enriched.mtx" > "$tmp/corrupt.mtx"
if python3 "$here/mfamg_mode_audit.py" "$tmp/A.mtx" "$tmp/base.mtx" \
  "$tmp/corrupt.mtx" --csv-out "$tmp/invalid.csv" > "$tmp/invalid.out" 2> "$tmp/invalid.err"; then
  echo 'ERROR: altered baseline was accepted' >&2
  exit 1
fi
grep -q '^MODE_AUDIT_ERROR=enriched first 2 modes differ from baseline' "$tmp/invalid.err"
test ! -e "$tmp/invalid.csv"
cat > "$tmp/badA.mtx" <<'BAD'
%%MatrixMarket matrix coordinate complex general
5 5 1
1 1 2 3
BAD
if python3 "$here/mfamg_mode_audit.py" "$tmp/badA.mtx" "$tmp/base.mtx" \
  "$tmp/enriched.mtx" --csv-out "$tmp/bad.csv" > "$tmp/bad.out" 2> "$tmp/bad.err"; then
  echo 'ERROR: unsupported matrix was accepted' >&2
  exit 1
fi
grep -q 'MODE_AUDIT_ERROR=' "$tmp/bad.err"
test ! -e "$tmp/bad.csv"
echo MODE_AUDIT_SMOKE_OK
