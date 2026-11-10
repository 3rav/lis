#!/usr/bin/env bash
# LOCAL UCRT64 Newton triple: read-only preflight first, optional real solve.
# Never guesses genuine Newton chronology from numbered A/b filenames.
# Uses existing exact upstream-#214 winfix1 EXE, not old bin6000.
# No upload, no checkout/reset/clean/merge, no existing output overwrite.
set -euo pipefail

cd "$HOME/lis-v8"
ROOT="$HOME/mfamg_replay"
BRANCH="experiment/mfamg-fgcrodr-replay-v8-20261108"
REF="refs/remotes/origin/$BRANCH"
EXE="$ROOT/build-v8-upstream214-20261109-winfix1/test/.libs/mfamgreplay.exe"
SOLVER="$ROOT/src-v8-upstream214-20261109-winfix1/src/solver/lis_solver_fgcrodr.c"
EXPECTED_BLOB="815a73ada1d18d970ec9b6d24306d23ff631fb35"
CSV="$(printenv MFAMG_V8_NEWTON_STEPS_CSV || true)"
MODES="$(printenv MFAMG_V8_NEWTON_MODES || true)"
FACTOR="$(printenv MFAMG_V8_NEWTON_FINE_AGG_FACTOR || echo 4)"
EXECUTE="$(printenv MFAMG_V8_NEWTON_EXECUTE || echo 0)"
CONFIRMED="$(printenv MFAMG_V8_NEWTON_CONSECUTIVE_CONFIRMED || echo 0)"

if [[ "$(printenv MSYSTEM || true)" != "UCRT64" ]]; then
  echo "STOP: requires MSYS2 UCRT64"
  exit 2
fi
for cmd in git tar cygpath python3; do
  command -v "$cmd" >/dev/null || { echo "STOP: missing $cmd"; exit 2; }
done

# This is metadata discovery only, no filename-based chronology invention.
if [[ -z "$CSV" ]]; then
  count=0
  for candidate in \
    "$ROOT/newton_steps.csv" \
    "$ROOT/newton/newton_steps.csv" \
    "$ROOT/challenge_31/newton_steps.csv"; do
    if [[ -s "$candidate" ]]; then
      count=$((count+1))
      CSV="$candidate"
    fi
  done
  if (( count != 1 )); then
    echo "V8_NEWTON_GATE=NEEDS_EXPLICIT_CURATED_STEPS_CSV"
    echo "V8_NEWTON_CSV_CANDIDATES_FOUND=$count"
    echo "Set MFAMG_V8_NEWTON_STEPS_CSV to a local ordered step,matrix,rhs CSV."
    echo "Do NOT infer Newton chronology from A01/A02/A03 filenames."
    exit 3
  fi
fi

# Default 12 research near-nullspace modes for the established 31-spring
# context ONLY; shape and provenance are still checked at preflight.
if [[ -z "$MODES" ]]; then
  MODES="$ROOT/challenge_31/modes04.mtx"
fi
if [[ "$FACTOR" != 2 && "$FACTOR" != 4 ]]; then
  echo "STOP: MFAMG factor must be 2 or 4"
  exit 2
fi
if [[ "$EXECUTE" != 0 && "$EXECUTE" != 1 ]]; then
  echo "STOP: MFAMG_V8_NEWTON_EXECUTE must be 0 or 1"
  exit 2
fi
for file in "$EXE" "$SOLVER" "$CSV" "$MODES"; do
  if [[ ! -s "$file" ]]; then
    echo "V8_NEWTON_GATE=MISSING_REQUIRED_LOCAL_FILE"
    exit 2
  fi
done
if [[ "$EXECUTE" == 1 ]]; then
  PY_SRC="$ROOT/v8_newton_triplet_py_20261110"
else
  PY_SRC="$ROOT/v8_newton_preflight_py_20261110"
fi
OUT="$ROOT/v8_newton_triplet_upstream214_20261110"
if [[ -e "$PY_SRC" || -e "$OUT" ]]; then
  echo "STOP: new Python source or results directory already exists, no overwrite"
  exit 2
fi
if [[ "$(git hash-object "$SOLVER")" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: selected local solver source is not accepted upstream develop #214"
  exit 2
fi
git fetch --no-tags origin "refs/heads/$BRANCH:$REF"
if [[ "$(git show "$REF:src/solver/lis_solver_fgcrodr.c" | git hash-object --stdin)" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: remote solver differs from accepted upstream develop #214"
  exit 2
fi

# Archive only vetted public source. Original A/b/modes/manifest are not
# copied, transformed or written. Never materialize tracked COPYING symlink.
umask 077
mkdir "$PY_SRC"
git archive "$REF" \
  test/mfamg_mode_audit.py \
  test/mfamg_replay_v8_repeats.py \
  test/mfamg_v8_snapshot_sequence.py \
  test/mfamg_v8_triplet_windows.py \
  test/mfamg_v8_newton_preflight.py | tar -xf - -C "$PY_SRC"
cd "$PY_SRC"
echo "V8_NEWTON_PREFLIGHT_START"
python3 test/mfamg_v8_newton_preflight.py \
  --exe "$(cygpath -m "$EXE")" \
  --modes "$(cygpath -m "$MODES")" \
  --snapshots "$(cygpath -m "$CSV")" \
  --max-windows 1

if [[ "$EXECUTE" != 1 ]]; then
  echo "V8_NEWTON_GATE=PREFLIGHT_ONLY_NOT_A_REAL_SOLVE"
  echo "V8_NEWTON_GATE_NOTE=CSV_ORDER_IS_USER_DECLARED_NOT_PROVEN"
  exit 0
fi
if [[ "$CONFIRMED" != 1 ]]; then
  echo "V8_NEWTON_GATE=NEEDS_EXPLICIT_CHRONOLOGY_CONFIRMATION"
  echo "Confirm CSV contains three ACTUALLY consecutive Newton iterations."
  exit 3
fi

echo "V8_NEWTON_TRIPLET_SAME_SOLVER_START"
set +e
python3 test/mfamg_v8_triplet_windows.py \
  --exe "$(cygpath -m "$EXE")" \
  --modes "$(cygpath -m "$MODES")" \
  --snapshots "$(cygpath -m "$CSV")" \
  --fine-agg-factor "$FACTOR" --max-windows 1 --fresh-third \
  --out-dir "$(cygpath -m "$OUT")"
RC=$?
set -e
echo "V8_NEWTON_TRIPLET_EXIT=$RC"
if [[ -f "$OUT/summary.json" ]]; then
  echo "=== LOCAL AGGREGATE NEWTON TRIPLET SUMMARY ==="
  cat "$OUT/summary.json"
fi
exit "$RC"
