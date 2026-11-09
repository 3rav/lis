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

unset LIS_MFAMG_REPLAY_NEXT_A LIS_MFAMG_REPLAY_NEXT_B
"$exe" A5.mtx b2.mtx modes.mtx > fresh.log 2>&1
grep -q '^REPLAY_RESULT solve_rc=0 status=0 ' fresh.log

for item in "separate.log:REPLAY_NEWTON_RESULT" \
            "inplace.log:REPLAY_NEWTON_RESULT" \
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
echo MFAMG_V8_CROSS_OPERATOR_REJECT_AND_TRANSPORT_OK
