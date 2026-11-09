#!/usr/bin/env bash
# Synthetic 240-row A1->A2 fixed-structure cross-operator regression.
# No private contact matrices. Uses only Bash, awk and the built replay EXE.
set -euo pipefail
if (( $# != 2 )); then
  echo "usage: $0 /absolute/path/mfamgreplay NEW-output-directory" >&2
  exit 2
fi
exe="$1"
out="$2"
test -s "$exe"
if [[ -e "$out" ]]; then
  echo "output exists: $out" >&2
  exit 2
fi
mkdir -m 700 "$out"
cd "$out"
for diagonal in 4 5; do
  path="A${diagonal}.mtx"
  {
    echo '%%MatrixMarket matrix coordinate real general'
    echo '240 240 718'
    for ((i=1; i<=240; i++)); do
      echo "$i $i $diagonal"
      if (( i > 1 )); then echo "$i $((i-1)) -1"; fi
      if (( i < 240 )); then echo "$i $((i+1)) -1"; fi
    done
  } > "$path"
done
# Same-CSR graph as A4 but a genuinely nonsymmetric off-diagonal A2.
{
  echo '%%MatrixMarket matrix coordinate real general'
  echo '240 240 718'
  for ((i=1; i<=240; i++)); do
    echo "$i $i 5"
    if (( i > 1 )); then echo "$i $((i-1)) -0.9"; fi
    if (( i < 240 )); then echo "$i $((i+1)) -1.1"; fi
  done
} > A_mixed.mtx
# Identical dimensions and nnz but different column graph in row 1.
{
  echo '%%MatrixMarket matrix coordinate real general'
  echo '240 240 718'
  for ((i=1; i<=240; i++)); do
    echo "$i $i 5"
    if (( i > 1 )); then echo "$i $((i-1)) -0.9"; fi
    if (( i < 240 )); then
      if (( i == 1 )); then
        echo '1 3 -1.1'
      else
        echo "$i $((i+1)) -1.1"
      fi
    fi
  done
} > A_bad_graph.mtx
for step in 1 2; do
  {
    echo '%%MatrixMarket matrix array real general'
    echo '240 1'
    for ((i=1; i<=240; i++)); do
      if (( step == 1 )); then
        echo "$(( (i*37)%23 - 11 ))"
      else
        echo "$(( (i*29)%19 - 9 ))"
      fi
    done
  } > "b${step}.mtx"
done
{
  echo '%%MatrixMarket matrix array real general'
  echo '240 2'
  for ((i=1; i<=240; i++)); do echo 1; done
  for ((i=1; i<=240; i++)); do
    awk -v i="$i" 'BEGIN {printf "%.17g\n",(i-120.5)/240}'
  done
} > modes.mtx

# Read-only preflight must compare exactly the same imported CSR arrays
# that the later value-replacement path uses. No solve, no state mutation.
"$exe" --compare-csr A4.mtx A_mixed.mtx > graph_ok.log 2>&1
grep -q '^REPLAY_CSR_COMPARE compatible=1 n=240 nnz=718 changed_nnz=718$' graph_ok.log
! grep -q '^REPLAY_RESULT' graph_ok.log
"$exe" --compare-csr A4.mtx A5.mtx > graph_diagonal.log 2>&1
grep -q '^REPLAY_CSR_COMPARE compatible=1 n=240 nnz=718 changed_nnz=240$' graph_diagonal.log
if "$exe" --compare-csr A4.mtx A_bad_graph.mtx > graph_bad.log 2>&1; then
  echo 'ERROR: CSR preflight accepted incompatible column graph' >&2
  exit 1
fi
grep -q '^REPLAY_CSR_COMPARE compatible=0$' graph_bad.log
grep -q 'replace-values CSR column graph mismatch' graph_bad.log
if "$exe" --compare-csr A4.mtx A4.mtx > graph_unchanged.log 2>&1; then
  echo 'ERROR: CSR preflight accepted unchanged operator as changed' >&2
  exit 1
fi
grep -q '^REPLAY_CSR_COMPARE compatible=0$' graph_unchanged.log
grep -q 'replace-values requires a genuinely changed operator' graph_unchanged.log
if "$exe" --compare-csr A4.mtx nonexistent.mtx > graph_missing.log 2>&1; then
  echo 'ERROR: CSR preflight accepted missing A2' >&2
  exit 1
fi

export LIS_MFAMG_EXPERIMENT_FINE_AGG=1
export LIS_MFAMG_FINE_AGG_FACTOR=2
export LIS_MFAMG_EXPERIMENT_DEEP_MERGE=1
export LIS_MFAMG_VCYCLE_SMOOTH_STEPS=1
export LIS_MFAMG_REPLAY_RESTART=2
export LIS_MFAMG_REPLAY_MAXITER=100
export LIS_FGMRES_PROJECTED_RELATIVE=1
export LIS_MFAMG_REPLAY_METHOD=fgcrodr
export LIS_MFAMG_REPLAY_RECYCLE_DIM=1
export LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1
export LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1
export LIS_MFAMG_REPLAY_NEXT_A=A5.mtx
export LIS_MFAMG_REPLAY_NEXT_B=b2.mtx

"$exe" A4.mtx b1.mtx modes.mtx > separate.log 2>&1
grep -q '^REPLAY_RESULT solve_rc=0 status=0 ' separate.log
grep -q '^REPLAY_NEWTON_MATRIX mode=separate_object' separate.log
grep -q '^REPLAY_NEWTON_RESULT step=2 solve_rc=0 status=0 ' separate.log
grep -q '^REPLAY_NEWTON_RECYCLE_CHECK rc=0 ' separate.log
grep -Eq '^REPLAY_NEWTON_POLICY cross_uses_delta=0 keeps_delta=0 transports_delta=0 rejects_delta=[1-9][0-9]*$' separate.log

# SAME matrix object, changed value generation: the old Z must carry.
LIS_MFAMG_REPLAY_NEXT_DIAG_SHIFT=-1 \
  "$exe" A4.mtx b1.mtx modes.mtx > inplace.log 2>&1
grep -q '^REPLAY_NEWTON_MATRIX mode=inplace_diag_shift shift=-1' inplace.log
grep -q '^REPLAY_NEWTON_RESULT step=2 solve_rc=0 status=0 ' inplace.log
grep -q '^REPLAY_NEWTON_RECYCLE_CHECK rc=0 ' inplace.log
grep -Eq '^REPLAY_NEWTON_POLICY cross_uses_delta=[1-9][0-9]* ' inplace.log

# Test a genuine same-graph arbitrary-value A2, including off-diagonals.
LIS_MFAMG_REPLAY_NEXT_A=A_mixed.mtx LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1 \
  "$exe" A4.mtx b1.mtx modes.mtx > replacement.log 2>&1
grep -Eq '^REPLAY_NEWTON_MATRIX mode=inplace_replace_values changed_nnz=[1-9][0-9]*$' replacement.log
grep -q '^REPLAY_NEWTON_RESULT step=2 solve_rc=0 status=0 ' replacement.log
grep -q '^REPLAY_NEWTON_RECYCLE_CHECK rc=0 ' replacement.log
grep -Eq '^REPLAY_NEWTON_POLICY cross_uses_delta=[1-9][0-9]* ' replacement.log

unset LIS_MFAMG_REPLAY_NEXT_A LIS_MFAMG_REPLAY_NEXT_B
"$exe" A5.mtx b2.mtx modes.mtx > fresh.log 2>&1
"$exe" A_mixed.mtx b2.mtx modes.mtx > fresh_mixed.log 2>&1
grep -q '^REPLAY_RESULT solve_rc=0 status=0 ' fresh_mixed.log
grep -q '^REPLAY_RESULT solve_rc=0 status=0 ' fresh.log

for item in "separate.log:REPLAY_NEWTON_RESULT" \
            "inplace.log:REPLAY_NEWTON_RESULT" \
            "replacement.log:REPLAY_NEWTON_RESULT" \
            "fresh_mixed.log:REPLAY_RESULT" \
            "fresh.log:REPLAY_RESULT"; do
  log="${item%%:*}"
  marker="${item#*:}"
  awk -v prefix="$marker" '
    $1==prefix {
      count++
      for (i=2;i<=NF;i++) if ($i ~ /^true_relative=/) {
        split($i,a,"=")
        if (a[2] !~ /^[0-9.eE+-]+$/ ||
            (a[2]+0)>1e-5 || (a[2]+0)<0) bad=1
        found++
      }
    }
    END {if (count!=1 || found!=1 || bad) exit 1}
  ' "$log"
done

# Wrong stated shift must not pass independent A2 true residual check.
export LIS_MFAMG_REPLAY_NEXT_A=A5.mtx
export LIS_MFAMG_REPLAY_NEXT_B=b2.mtx
if LIS_MFAMG_REPLAY_NEXT_DIAG_SHIFT=1 \
    "$exe" A4.mtx b1.mtx modes.mtx > wrong_shift.log 2>&1; then
  echo "ERROR: wrong shift accepted against independent A2" >&2
  exit 1
fi
grep -q '^REPLAY_NEWTON_MATRIX mode=inplace_diag_shift shift=1' wrong_shift.log
if LIS_MFAMG_REPLAY_NEXT_B=nonexistent.mtx \
    "$exe" A4.mtx b1.mtx modes.mtx > missing_b2.log 2>&1; then
  echo "ERROR: missing second RHS passed" >&2
  exit 1
fi
grep -q '^REPLAY_RESULT solve_rc=0 status=0 ' missing_b2.log
# Detect graph mismatch BEFORE touching any old value.
if LIS_MFAMG_REPLAY_NEXT_A=A_bad_graph.mtx \
   LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1 \
   "$exe" A4.mtx b1.mtx modes.mtx > bad_graph.log 2>&1; then
  echo "ERROR: incompatible CSR graph accepted" >&2
  exit 1
fi
grep -q 'replace-values CSR column graph mismatch' bad_graph.log
! grep -q '^REPLAY_NEWTON_STEP step=2' bad_graph.log

if LIS_MFAMG_REPLAY_NEXT_A=A4.mtx \
   LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1 \
   "$exe" A4.mtx b1.mtx modes.mtx > identical.log 2>&1; then
  echo "ERROR: unchanged operator accepted as changed" >&2
  exit 1
fi
grep -q 'replace-values requires a genuinely changed operator' identical.log

if LIS_MFAMG_REPLAY_NEXT_A=A_mixed.mtx \
   LIS_MFAMG_REPLAY_NEXT_DIAG_SHIFT=-1 \
   LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1 \
   "$exe" A4.mtx b1.mtx modes.mtx > conflicting.log 2>&1; then
  echo "ERROR: conflicting mutation modes accepted" >&2
  exit 1
fi
grep -q 'replace-values requires NEXT_A and NEXT_B' conflicting.log
echo MFAMG_V8_CROSS_OPERATOR_REJECT_AND_TRANSPORT_OK
