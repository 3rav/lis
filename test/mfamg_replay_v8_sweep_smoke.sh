#!/usr/bin/env bash
# End-to-end synthetic test: zero private snapshot content.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$tmp/data/challenge_31" "$tmp/data/control_116" "$tmp/data/results"
for f in challenge_31/A04.mtx challenge_31/b04.mtx challenge_31/modes04.mtx \
         challenge_31/modes04_12plus6negative.mtx control_116/A02.mtx \
         control_116/b02.mtx control_116/modes02.mtx; do
  echo placeholder > "$tmp/data/$f"
done
cat > "$tmp/fake_replay.sh" <<'FAKE'
#!/usr/bin/env bash
set -euo pipefail
method="$LIS_MFAMG_REPLAY_METHOD"
if [[ "$method" == invalid ]]; then
  echo "REPLAY: LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
k=0
if [[ "$method" == fgcrodr ]]; then k="$LIS_MFAMG_REPLAY_RECYCLE_DIM"; fi
n=60265; iter=1000; status=4; relative=1.832690841870e-03; rc=1
if [[ "$1" == control_116/* ]]; then
  n=95000; iter=592; status=0; relative=9.000000000000e-06; rc=0
  if [[ "$method" == fgcrodr ]]; then iter=580; fi
elif [[ "$method" == fgcrodr && "$k" != 0 ]]; then
  relative=1.100000000000e-03
fi
modes=12
if [[ "$3" == *12plus6negative* ]]; then modes=18; fi
echo "REPLAY n=$n nnz=123456 modes=$modes restart=180 maxiter=1000 tol=1e-5 method=$method recycle_k=$k"
if [[ "$method" == fgcrodr ]]; then
  echo "REPLAY_RECYCLE dim=$k action=2 quality=1 min_quality=1"
  echo "REPLAY_RECYCLE_STATS builds=1 updates=2 uses=1 rejects=0"
  if (( k > 0 )); then
    echo "LIS_FGCRODR_INTRA_RECYCLE iter=180 k=$k relative=0.1"
    echo "LIS_FGCRODR_INTRA_RECYCLE iter=360 k=$k relative=0.01"
    echo "REPLAY_RECYCLE_CHECK rc=0 orth_error=4.4e-16 relation_error=1.8e-13"
  fi
fi
echo "REPLAY_RESULT solve_rc=$status status=$status iter=$iter reported_residual=$relative true_relative=$relative wall_seconds=12.500"
exit "$rc"
FAKE
chmod +x "$tmp/fake_replay.sh"
MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_replay.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results" MFAMG_V8_FULL=1 \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/sweep.stdout"
csv="$tmp/data/results/v8_summary.csv"
test -s "$csv"
test "$(wc -l < "$csv")" -eq 10
grep -q '^31_12m_f2_fgmres,31,60265,12,2,fgmres,0,180,1000,1,4,4,1000,' "$csv"
grep -q '^31_12m_f2_fgcrodr_k0,31,60265,12,2,fgcrodr,0,180,1000,1,4,4,1000,' "$csv"
grep -q '^116_12m_f2_fgmres,116,95000,12,2,fgmres,0,180,1000,0,0,0,592,' "$csv"
grep -q '^116_12m_f2_fgcrodr_k8,116,95000,12,2,fgcrodr,8,180,1000,0,0,0,580,' "$csv"
grep -q ',8,1,2,1,0,2,NA,NA,NA,v8_116_12m_f2_fgcrodr_k8.log$' "$csv"
test "$(find "$tmp/data/results" -maxdepth 1 -name '*.log' | wc -l)" -eq 9
test -s "$tmp/data/results/v8_assessment.txt"
grep -q '^ASSESS_NUMERIC_GATE=PASS_CONTACT_REPLAY_ONLYcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
test -s "$tmp/data/results_failed/v8_assessment.txt"
grep -q '^ASSESS_NUMERIC_GATE=INCOMPLETE_OR_INVALID_DATAgrep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/results/v8_assessment.txt"
grep -q '^ASSESS_PRODUCTION_GATE=BLOCKED_NEEDS_REPEATED_TIMING_AND_NEWTON_SEQUENCEcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/results/v8_assessment.txt"
awk -F, 'BEGIN {OFS=","} /^116_12m_f2_fgcrodr_k8,/ {$13=700} {print}' "$csv" > "$tmp/data/numerically_worse.csv"
awk -f "$here/mfamg_replay_v8_assess.awk" "$tmp/data/numerically_worse.csv" > "$tmp/data/worse_assessment.txt"
grep -q '^ASSESS_116_K8_CONTROL=FAILcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/worse_assessment.txt"
grep -q '^ASSESS_NUMERIC_GATE=FAIL_OR_INCONCLUSIVEcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/worse_assessment.txt"
cat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/results_failed/v8_assessment.txt"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/results/v8_assessment.txt"
grep -q '^ASSESS_PRODUCTION_GATE=BLOCKED_NEEDS_REPEATED_TIMING_AND_NEWTON_SEQUENCEcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/results/v8_assessment.txt"
awk -F, 'BEGIN {OFS=","} /^116_12m_f2_fgcrodr_k8,/ {$13=700} {print}' "$csv" > "$tmp/data/numerically_worse.csv"
awk -f "$here/mfamg_replay_v8_assess.awk" "$tmp/data/numerically_worse.csv" > "$tmp/data/worse_assessment.txt"
grep -q '^ASSESS_116_K8_CONTROL=FAILcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/worse_assessment.txt"
grep -q '^ASSESS_NUMERIC_GATE=FAIL_OR_INCONCLUSIVEcat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
 "$tmp/data/worse_assessment.txt"
cat > "$tmp/fake_crash.sh" <<'CRASH'
#!/usr/bin/env bash
if [[ "$LIS_MFAMG_REPLAY_METHOD" == invalid ]]; then
  echo "LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr" >&2
  exit 1
fi
if [[ "$LIS_MFAMG_REPLAY_METHOD" == fgcrodr && "$LIS_MFAMG_REPLAY_RECYCLE_DIM" == 8 ]]; then
  echo "simulated crash before REPLAY_RESULT" >&2
  exit 17
fi
exec "$(dirname "$0")/fake_replay.sh" "$@"
CRASH
chmod +x "$tmp/fake_crash.sh"
mkdir "$tmp/data/results_failed"
if MFAMG_REPLAY_DATA="$tmp/data" MFAMG_REPLAY_EXE="$tmp/fake_crash.sh" \
  MFAMG_V8_RESULTS_DIR="$tmp/data/results_failed" \
  bash "$here/mfamg_replay_fgcrodr_sweep.sh" > "$tmp/failure.stdout" 2>"$tmp/failure.stderr"; then
  echo "ERROR: malformed case was not signaled" >&2
  exit 1
fi
test "$(find "$tmp/data/results_failed" -name '*.log' | wc -l)" -eq 6
test -s "$tmp/data/results_failed/v8_summary.csv"
grep -q '^116_12m_f2_fgmres,116,95000,' "$tmp/data/results_failed/v8_summary.csv"
grep -q 'NA,NA,NA,NA,NA,NA' "$tmp/data/results_failed/v8_summary.csv"
echo "REPLAY_SWEEP_END_TO_END_SMOKE_OK"
