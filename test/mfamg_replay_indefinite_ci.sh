#!/usr/bin/env bash
# MFAMG v8 pilot: deterministic real symmetric indefinite system.
# Strict row diagonal dominance makes the matrix nonsingular; positive
# and negative diagonal entries prove indefiniteness.
# Usage: bash test/mfamg_replay_indefinite_ci.sh ABS_EXE ABS_OUTPUT_DIR
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
{
  echo '%%MatrixMarket matrix coordinate real general'
  echo "$n $n $((3*n-2))"
  for ((i=1;i<=n;i++)); do
    diag=3.25
    if (( (i>=121 && i<=180) || (i>=301 && i<=360) )); then
      diag=-3.25
    fi
    printf '%d %d %s\n' "$i" "$i" "$diag"
    if (( i>1 )); then printf '%d %d -1\n' "$i" "$((i-1))"; fi
    if (( i<n )); then printf '%d %d -1\n' "$i" "$((i+1))"; fi
  done
} > A_indef.mtx
{
  echo '%%MatrixMarket matrix array real general'
  echo "$n 1"
  for ((i=1;i<=n;i++)); do echo "$(( (i*37)%29-14 ))"; done
} > b_indef.mtx
{
  echo '%%MatrixMarket matrix array real general'
  echo "$n 2"
  for ((i=1;i<=n;i++)); do echo 1; done
  awk -v n="$n" 'BEGIN {for(i=1;i<=n;i++) printf "%.17g\n",(i-(n+1)/2)/n}'
} > modes_indef.mtx
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
  log="indef_$label.log"
  set +e
  env "$@" "$exe" A_indef.mtx b_indef.mtx modes_indef.mtx > "$log" 2>&1
  rc=$?
  set -e
  echo "INDEFINITE_CASE name=$label process_exit=$rc"
  grep -E '^REPLAY n=|^REPLAY_RECYCLE|^REPLAY_RESULT|^LIS_FGCRODR_INTRA_RECYCLE ' "$log" || true
  grep -q '^REPLAY_RESULT ' "$log" || {
    tail -n 30 "$log"
    echo "ERROR: no solver result in $log" >&2
    return 1
  }
  if grep -qiE 'nan|[+-]inf' "$log"; then
    echo "ERROR: non-finite text in $log" >&2
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
# k=0 must remain a safe no-op even when intra-solve is requested.
run_case k0_optin LIS_MFAMG_REPLAY_METHOD=fgcrodr LIS_MFAMG_REPLAY_RECYCLE_DIM=0 LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1 LIS_FGCRODR_INTRA_TRACE=1
run_case k1 LIS_MFAMG_REPLAY_METHOD=fgcrodr LIS_MFAMG_REPLAY_RECYCLE_DIM=1 LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1 LIS_FGCRODR_INTRA_TRACE=1 LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1
run_case k2 LIS_MFAMG_REPLAY_METHOD=fgcrodr LIS_MFAMG_REPLAY_RECYCLE_DIM=2 LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1 LIS_FGCRODR_INTRA_TRACE=1 LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1
# k=0 must be numerically identical to FGMRES.
signature() {
  grep '^REPLAY_RESULT ' "$1" |
    sed -E 's/^.*status=([^ ]+) iter=([^ ]+) reported_residual=([^ ]+) true_relative=([^ ]+) wall_seconds=.*$/\1 \2 \3 \4/'
}
test "$(signature indef_fgmres.log)" = "$(signature indef_k0.log)" || {
  echo "ERROR: k=0 differs from FGMRES baseline" >&2
  exit 1
}
test "$(signature indef_fgmres.log)" = "$(signature indef_k0_optin.log)" || {
  echo "ERROR: opt-in k=0 differs from FGMRES baseline" >&2
  exit 1
}
awk '/^REPLAY_RESULT / {for(i=1;i<=NF;i++) if($i ~ /^iter=/) {
       split($i,p,"="); if(p[2]>4) ok=1
     }} END {exit !ok}' indef_fgmres.log || {
  echo "ERROR: baseline did not reach an Arnoldi restart" >&2
  exit 1
}
for label in k1 k2; do
  log="indef_$label.log"
  grep -q '^LIS_FGCRODR_INTRA_RECYCLE iter=' "$log" || {
    echo "ERROR: no within-solve reuse in $log" >&2
    exit 1
  }
  grep -q '^REPLAY_RECYCLE_CHECK rc=0 ' "$log" || {
    echo "ERROR: invalid recycle-space algebra in $log" >&2
    exit 1
  }
done
# Only convergence and algebra are gated; do NOT claim speedup.
echo 'INDEFINITE_PROBE_COMPLETE: examine status, iter, true_relative, reuse.'
