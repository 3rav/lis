#!/usr/bin/env bash
# Five new offline Ritz leave-one-out replays; case 06 already tested.
set -euo pipefail
data="${MFAMG_REPLAY_DATA:-$HOME/mfamg_replay}"
exe="${MFAMG_REPLAY_EXE:-$HOME/lis-v8/bin6000/mfamgreplay.exe}"
out="${MFAMG_RITZ_LOO_RESULTS_DIR:-$data/v8_mode_audit/loo/logs}"
dir="$data/v8_mode_audit/loo"
A="$data/challenge_31/A04.mtx"
b="$data/challenge_31/b04.mtx"

for f in "$exe" "$A" "$b"; do
  if [[ ! -s "$f" ]]; then
    echo "RITZ_LOO_REPLAY_ERROR=missing_file path=$f" >&2
    exit 2
  fi
done
mkdir -p "$out"
for i in 01 02 03 04 05; do
  if [[ ! -s "$dir/modes_omit_ritz_${i}.mtx" ]]; then
    echo "RITZ_LOO_REPLAY_ERROR=missing_modes case=$i" >&2
    exit 2
  fi
  if [[ -e "$out/omit_${i}_fgcrodr_3000.log" ]]; then
    echo "RITZ_LOO_REPLAY_ERROR=log_already_exists case=$i" >&2
    exit 2
  fi
done

export LIS_MFAMG_EXPERIMENT_FINE_AGG=1
export LIS_MFAMG_EXPERIMENT_DEEP_MERGE=1
export LIS_MFAMG_FINE_AGG_FACTOR=4
export LIS_MFAMG_STAR_SCHUR=0
export LIS_MFAMG_VCYCLE_SMOOTH_STEPS=3
export LIS_MFAMG_VCYCLE_DIAG=0
export LIS_MFAMG_VCYCLE_TRACE=0
export LIS_MFAMG_DIAG=1
export LIS_FGMRES_TRACE=1
export LIS_FGMRES_PROJECTED_RELATIVE=1
export LIS_MFAMG_REPLAY_METHOD=fgcrodr
export LIS_MFAMG_REPLAY_RECYCLE_DIM=8
export LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1
export LIS_FGCRODR_INTRA_TRACE=1
export LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1
export LIS_MFAMG_REPLAY_RESTART=180
export LIS_MFAMG_REPLAY_MAXITER=3000
unset LIS_MFAMG_MODES_EXPORT_PREFIX

for i in 01 02 03 04 05; do
  mode="$dir/modes_omit_ritz_${i}.mtx"
  log="$out/omit_${i}_fgcrodr_3000.log"
  echo "===== RITZ LOO omit $i ====="
  if "$exe" "$A" "$b" "$mode" > "$log" 2>&1; then
    rc=0
  else
    rc=$?
  fi
  echo "process_exit=$rc"
  grep -E '^REPLAY n=|^REPLAY_RECYCLE_STATS|^REPLAY_RECYCLE_CHECK|^REPLAY_RESULT' "$log" || true
  if ! grep -q '^REPLAY_RESULT ' "$log"; then
    echo "RITZ_LOO_REPLAY_ERROR=no_REPLAY_RESULT case=$i log=$log" >&2
    exit 1
  fi
  if [[ "$rc" -ne 0 && "$rc" -ne 1 ]]; then
    echo "RITZ_LOO_REPLAY_ERROR=unexpected_process_exit case=$i rc=$rc" >&2
    exit 1
  fi
done
echo RITZ_LOO_REPLAY_COMPLETE
