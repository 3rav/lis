#!/usr/bin/env bash
# MFAMG v8: difficult near-singular *indefinite* synthetic regression pilot.
# Two disjoint 240x240 tridiagonal blocks: +T and -T, where
# T has diagonal 2+delta and +/- 1 neighbors (delta=1e-5).
# Exact |lambda|min = 2+delta - 2*cos(pi/241), max = 2+delta+2*cos(pi/241).
# The matrix is symmetric, indefinite, nonsingular and ill-conditioned.
# Gating convergence and algebraic regression; does not assert speedup.
# Usage: bash test/mfamg_replay_illconditioned_ci.sh ABS_EXE ABS_OUTPUT_DIR
set -euo pipefail
if (( $# != 2 )); then
  echo "usage: $0 ABS_EXE ABS_OUTPUT_DIR" >&2
  exit 2
fi
exe="$1"
output="$2"
test -s "$exe"
mkdir -p "$output"
cd "$output"
n=480
block=240
{
  echo '%%MatrixMarket matrix coordinate real general'
  echo "$n $n $((3*n-4))"
  for ((i=1;i<=n;i++)); do
    diag=2.00001
    off=-1
    if (( i>block )); then
      diag=-2.00001
      off=1
    fi
    printf '%d %d %s\n' "$i" "$i" "$diag"
    if (( i>1 && i!=block+1 )); then
      printf '%d %d %d\n' "$i" "$((i-1))" "$off"
    fi
    if (( i<n && i!=block )); then
      printf '%d %d %d\n' "$i" "$((i+1))" "$off"
    fi
  done
} > A_ill.mtx
{
  echo '%%MatrixMarket matrix array real general'
  echo "$n 1"
  for ((i=1;i<=n;i++)); do echo "$(( (i*37)%29-14 ))"; done
} > b_ill.mtx
{
  echo '%%MatrixMarket matrix array real general'
  echo "$n 2"
  for ((i=1;i<=n;i++)); do echo 1; done
  awk -v n="$n" 'BEGIN {for(i=1;i<=n;i++) printf "%.17g\n",(i-(n+1)/2)/n}'
} > modes_ill.mtx

export LIS_MFAMG_EXPERIMENT_FINE_AGG=1
export LIS_MFAMG_FINE_AGG_FACTOR=2
export LIS_MFAMG_EXPERIMENT_DEEP_MERGE=1
export LIS_MFAMG_VCYCLE_SMOOTH_STEPS=1
export LIS_MFAMG_STAR_SCHUR=0
export LIS_FGMRES_PROJECTED_RELATIVE=1
export LIS_MFAMG_REPLAY_RESTART=4
export LIS_MFAMG_REPLAY_MAXITER=300

run_case() {
  local label="$1" rc log
  shift
  log="ill_$label.log"
  set +e
  env "$@" "$exe" A_ill.mtx b_ill.mtx modes_ill.mtx > "$log" 2>&1
  rc=$?
  set -e
  echo "ILL_CASE name=$label process_exit=$rc"
  grep -E '^REPLAY n=|^REPLAY_RECYCLE|^REPLAY_RESULT|^LIS_FGCRODR_INTRA_RECYCLE ' "$log" | tail -n 20 || true
  if ! grep -q '^REPLAY_RESULT ' "$log"; then
    echo "ERROR: no REPLAY_RESULT in $log" >&2
    tail -n 40 "$log"
    return 1
  fi
  if grep -qiE 'nan|[+-]inf' "$log"; then
    echo "ERROR: nonfinite text in $log" >&2
    return 1
  fi
  if (( rc != 0 )) || ! grep -q '^REPLAY_RESULT solve_rc=0 status=0 ' "$log"; then
    echo "ERROR: solver did not converge in $log" >&2
    return 1
  fi
  awk '/^REPLAY_RESULT / {
         for(i=1;i<=NF;i++) if($i ~ /^true_relative=/) {
           split($i,p,"="); if(p[2]+0>=0 && p[2]+0<1.0e-5) ok=1
         }
       } END {exit !ok}' "$log" || {
    echo "ERROR: true residual above tolerance in $log" >&2
    return 1
  }
}

run_case fgmres LIS_MFAMG_REPLAY_METHOD=fgmres LIS_FGCRODR_INTRA_SOLVE_RECYCLE=0
run_case k0 LIS_MFAMG_REPLAY_METHOD=fgcrodr LIS_MFAMG_REPLAY_RECYCLE_DIM=0 LIS_FGCRODR_INTRA_SOLVE_RECYCLE=0
run_case k1 LIS_MFAMG_REPLAY_METHOD=fgcrodr LIS_MFAMG_REPLAY_RECYCLE_DIM=1 LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1 LIS_FGCRODR_INTRA_TRACE=1 LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1
run_case k2 LIS_MFAMG_REPLAY_METHOD=fgcrodr LIS_MFAMG_REPLAY_RECYCLE_DIM=2 LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1 LIS_FGCRODR_INTRA_TRACE=1 LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1

signature() {
  grep '^REPLAY_RESULT ' "$1" |
    sed -E 's/^.*status=([^ ]+) iter=([^ ]+) reported_residual=([^ ]+) true_relative=([^ ]+) wall_seconds=.*$/\1 \2 \3 \4/'
}
test "$(signature ill_fgmres.log)" = "$(signature ill_k0.log)" || {
  echo "ERROR: k=0 differs from FGMRES control on ill-conditioned problem" >&2
  exit 1
}
# Ensure the hard baseline really needs multiple restarts.
awk '/^REPLAY_RESULT / {for(i=1;i<=NF;i++) if($i ~ /^iter=/) {
       split($i,p,"="); if(p[2]>8) ok=1
     }} END {exit !ok}' ill_fgmres.log || {
  echo "ERROR: ill-conditioned baseline did not require two restarts" >&2
  exit 1
}
for label in k1 k2; do
  log="ill_$label.log"
  count=$(grep -c '^LIS_FGCRODR_INTRA_RECYCLE iter=' "$log" || true)
  if (( count < 2 )); then
    echo "ERROR: less than two active recycling restarts in $log" >&2
    exit 1
  fi
  grep -q '^REPLAY_RECYCLE_CHECK rc=0 ' "$log" || {
    echo "ERROR: invalid recycle-space algebra in $log" >&2
    exit 1
  }
  awk '/^REPLAY_RECYCLE_STATS / {
         for(i=1;i<=NF;i++) if($i ~ /^updates=/) {
           split($i,p,"="); if(p[2]+0>=2) ok=1
         }
       } END {exit !ok}' "$log" || {
    echo "ERROR: fewer than two recycle updates in $log" >&2
    exit 1
  }
done
echo "ILL_REGRESSION_COMPLETE: validated convergence and multiple active restarts."
