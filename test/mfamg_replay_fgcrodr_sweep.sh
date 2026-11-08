#!/usr/bin/env bash
# MFAMG v8: compare online harmonic-Ritz recycling against FGMRES.
# Optional: MFAMG_V8_FULL=1 bash mfamg_replay_fgcrodr_sweep.sh
# MFAMG_REPLAY_DATA points to challenge_31/ and control_116/.
set -uo pipefail

root="${MFAMG_REPLAY_DATA:-$HOME/mfamg_replay}"
exe="${MFAMG_REPLAY_EXE:-/ucrt64/bin/mfamgreplay.exe}"
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
cd -- "$root" || exit 2
root="$PWD"
export PATH="/ucrt64/bin:$PATH"

A31=challenge_31/A04.mtx
B31=challenge_31/b04.mtx
M12=challenge_31/modes04.mtx
M18=challenge_31/modes04_12plus6negative.mtx
A116=control_116/A02.mtx
B116=control_116/b02.mtx
M116=control_116/modes02.mtx

[[ -s "$exe" ]] || { echo "ERROR: missing EXE $exe" >&2; exit 2; }
for file in "$A31" "$B31" "$M12" "$M18" "$A116" "$B116" "$M116"; do
  [[ -s "$file" ]] || { echo "ERROR: missing data $root/$file" >&2; exit 2; }
done

probe=$(LIS_MFAMG_REPLAY_METHOD=invalid "$exe" "$A31" "$B31" "$M12" 2>&1 || true)
if [[ "$probe" != *"LIS_MFAMG_REPLAY_METHOD must be"* ]]; then
  echo "ERROR: replay EXE lacks v8 method selector" >&2
  exit 2
fi

export LIS_MFAMG_EXPERIMENT_FINE_AGG=1
export LIS_MFAMG_EXPERIMENT_DEEP_MERGE=1
export LIS_MFAMG_STAR_SCHUR=0
export LIS_MFAMG_VCYCLE_SMOOTH_STEPS=3
export LIS_MFAMG_VCYCLE_DIAG=0
export LIS_MFAMG_VCYCLE_TRACE=0
export LIS_MFAMG_DIAG=1
export LIS_FGMRES_TRACE=1
export LIS_FGMRES_PROJECTED_RELATIVE=1
export LIS_MFAMG_REPLAY_RESTART=180
export LIS_MFAMG_REPLAY_MAXITER=1000
unset LIS_MFAMG_MODES_EXPORT_PREFIX

# Keep private snapshot logs/results under the snapshot root (or local override).
results_dir="${MFAMG_V8_RESULTS_DIR:-$root/v8_results}"
mkdir -p -- "$results_dir" || exit 2
results_dir="$(cd -- "$results_dir" && pwd)"
case_logs=()
case_errors=0

run_case() {
  local id="$1" A="$2" b="$3" modes="$4" factor="$5" method="$6" k="$7"
  local rc snapshot=31 log="$results_dir/v8_${id}.log"
  [[ "$id" == 116_* ]] && snapshot=116
  case_logs+=("$log")
  export LIS_MFAMG_FINE_AGG_FACTOR="$factor"
  export LIS_MFAMG_REPLAY_METHOD="$method"
  if [[ "$method" == fgcrodr ]]; then
    export LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1
    export LIS_FGCRODR_INTRA_TRACE=1
  else
    unset LIS_FGCRODR_INTRA_SOLVE_RECYCLE LIS_FGCRODR_INTRA_TRACE
  fi
  if [[ "$method" == fgcrodr ]]; then
    export LIS_MFAMG_REPLAY_RECYCLE_DIM="$k"
  else
    unset LIS_MFAMG_REPLAY_RECYCLE_DIM
  fi
  echo
  echo "==== $id | method=$method | factor=$factor | k=$k ===="
  printf "REPLAY_SWEEP_META case=%s snapshot=%s factor=%s\n" "$id" "$snapshot" "$factor" > "$log"
  "$exe" "$A" "$b" "$modes" >>"$log" 2>&1
  rc=$?
  printf "REPLAY_SWEEP_EXIT process_exit=%d\n" "$rc" >>"$log"
  echo "process_exit=$rc log=$log bytes=$(wc -c < "$log")"
  grep -E '^REPLAY n=|^REPLAY_RECYCLE|^REPLAY_RESULT|^REPLAY LIS error|^REPLAY:|^LIS_FGCRODR_INTRA_RECYCLE ' "$log" | tail -n 10 || true
  grep -q '^REPLAY_RESULT ' "$log" || {
    echo "ERROR: no REPLAY_RESULT; inspect $log" >&2
    tail -n 30 "$log"
    return 2
  }
  if [[ "$method" == fgcrodr && "$k" != 0 ]]; then
    if ! grep -q '^REPLAY_RECYCLE ' "$log"; then
      echo "ERROR: no recycle-space diagnostics in $log" >&2
      return 2
    fi
  fi
  if ! grep -q "method=$method recycle_k=$k" "$log"; then
    echo "ERROR: method confirmation failed in $log" >&2
    return 2
  fi
  return 0
}

run_case 31_12m_f2_fgmres "$A31" "$B31" "$M12" 2 fgmres 0 || case_errors=1
run_case 31_12m_f2_fgcrodr_k0 "$A31" "$B31" "$M12" 2 fgcrodr 0 || case_errors=1
run_case 31_12m_f2_fgcrodr_k8 "$A31" "$B31" "$M12" 2 fgcrodr 8 || case_errors=1
run_case 31_18m_f4_fgcrodr_k8 "$A31" "$B31" "$M18" 4 fgcrodr 8 || case_errors=1
run_case 116_12m_f2_fgmres "$A116" "$B116" "$M116" 2 fgmres 0 || case_errors=1
run_case 116_12m_f2_fgcrodr_k8 "$A116" "$B116" "$M116" 2 fgcrodr 8 || case_errors=1

if [[ "${MFAMG_V8_FULL:-0}" == 1 ]]; then
  run_case 31_12m_f2_fgcrodr_k16 "$A31" "$B31" "$M12" 2 fgcrodr 16 || case_errors=1
  run_case 31_18m_f4_fgmres "$A31" "$B31" "$M18" 4 fgmres 0 || case_errors=1
  run_case 31_18m_f4_fgcrodr_k16 "$A31" "$B31" "$M18" 4 fgcrodr 16 || case_errors=1
fi
 summary_errors=0
awk -f "$script_dir/mfamg_replay_v8_report.awk" "${case_logs[@]}" >"$results_dir/v8_summary.csv" || summary_errors=1
echo "MFAMG v8 finished: CSV=$results_dir/v8_summary.csv"
cat "$results_dir/v8_summary.csv"
echo "Reminder: process_exit and solve_rc do NOT imply status=0 or convergence."

if (( case_errors || summary_errors )); then
  echo "ERROR: malformed case; logs and CSV retained for inspection." >&2
  exit 2
fi
