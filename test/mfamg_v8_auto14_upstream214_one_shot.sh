#!/usr/bin/env bash
# PRIVATE, optional 31-spring A04/b04 automatic-mode feasibility control.
# Uses the EXISTING accepted upstream develop #214 winfix1 replay EXE.
# Compares 12 baseline vs 12+2 SMALL-SINGULAR automatically filtered modes
# vs 18 OFFLINE Ritz modes; the generator never reads the 18-mode oracle.
# The offline control is passed ONLY to the independent comparison phase.
# No private A/b/modes/logs are uploaded, and no existing files overwritten.
set -euo pipefail
cd "$HOME/lis-v8"

ROOT="$HOME/mfamg_replay"
BRANCH="experiment/mfamg-fgcrodr-replay-v8-20261108"
REF="refs/remotes/origin/$BRANCH"
EXE="$ROOT/build-v8-upstream214-20261109-winfix1/test/.libs/mfamgreplay.exe"
SOLVER="$ROOT/src-v8-upstream214-20261109-winfix1/src/solver/lis_solver_fgcrodr.c"
SOLVER_BLOB="815a73ada1d18d970ec9b6d24306d23ff631fb35"
A="$ROOT/challenge_31/A04.mtx"
B="$ROOT/challenge_31/b04.mtx"
BASE="$ROOT/challenge_31/modes04.mtx"
OFFLINE="$ROOT/v8_mode_audit/ritz/modes_all_ritz.mtx"
PY_SRC="$ROOT/v8_auto14_scripts_upstream214_20261110"
AUTO="$ROOT/v8_auto14_small_singular_upstream214_20261110.mtx"
OUT="$ROOT/v8_auto14_compare_upstream214_20261110"

if [[ "$(printenv MSYSTEM || true)" != "UCRT64" ]]; then
  echo "STOP: requires MSYS2 UCRT64"
  exit 2
fi
for command in git tar python3 cygpath; do
  command -v "$command" >/dev/null || { echo "STOP: missing $command"; exit 2; }
done
for file in "$EXE" "$SOLVER" "$A" "$B" "$BASE" "$OFFLINE"; do
  if [[ ! -s "$file" ]]; then
    echo "STOP: missing/empty private input or upstream replay EXE"
    exit 2
  fi
done
for path in "$PY_SRC" "$AUTO" "$OUT"; do
  if [[ -e "$path" || -L "$path" ]]; then
    echo "STOP: new scripts/modes/results already exist; no overwrite"
    exit 2
  fi
done
if [[ "$(git hash-object "$SOLVER")" != "$SOLVER_BLOB" ]]; then
  echo "STOP: local solver source does not match accepted upstream develop #214"
  exit 2
fi
git fetch --no-tags origin "refs/heads/$BRANCH:$REF"
if [[ "$(git show "$REF:src/solver/lis_solver_fgcrodr.c" | git hash-object --stdin)" != "$SOLVER_BLOB" ]]; then
  echo "STOP: remote solver source diverged from accepted upstream develop #214"
  exit 2
fi

# Fetch ONLY public scripts. Generator gets no path to the offline 18-mode
# basis, ensuring the candidate is not an oracle vector copied from A04 Ritz.
umask 077
mkdir "$PY_SRC"
git archive "$REF" \
  test/mfamg_mode_audit.py \
  test/mfamg_mode_ritz_selection.py \
  test/mfamg_mode_auto_candidates.py \
  test/mfamg_replay_v8_repeats.py \
  test/mfamg_replay_v8_candidate_compare.py | tar -xf - -C "$PY_SRC"
cd "$PY_SRC"
# Fail before the potentially expensive generator when the offline
# control is not a true baseline-preserving 18-column reference.
# This comparison does not feed any offline vector into the generator.
python3 - "$(cygpath -m "$BASE")" "$(cygpath -m "$OFFLINE")" <<'PY'
import sys
sys.path.insert(0, "test")
from mfamg_mode_audit import read_modes
n0, base = read_modes(sys.argv[1])
n1, offline = read_modes(sys.argv[2])
if n0 != n1 or len(base) != 12 or len(offline) != 18 or offline[:12] != base:
    print("V8_AUTO14_PREFLIGHT=BLOCKED_INVALID_12_18_BASELINE")
    raise SystemExit(2)
print("V8_AUTO14_PREFLIGHT=BASELINE12_OFFLINE18_MATCH")
PY
echo "V8_AUTO14_GENERATOR_START"
set +e
python3 test/mfamg_mode_auto_candidates.py \
  "$(cygpath -m "$A")" "$(cygpath -m "$BASE")" \
  --output "$(cygpath -m "$AUTO")" \
  --target small-singular --count 2 --steps 60 --seed 20261109
GEN_RC=$?
set -e
echo "V8_AUTO14_GENERATOR_EXIT=$GEN_RC"
if [[ "$GEN_RC" != 0 ]]; then
  echo "V8_AUTO14_ASSESSMENT=INCONCLUSIVE_GENERATOR_FAILED"
  echo "No fabricated candidate modes; private inputs remain unchanged."
  exit "$GEN_RC"
fi
if [[ ! -s "$AUTO" ]]; then
  echo "V8_AUTO14_ASSESSMENT=INCONCLUSIVE_MISSING_MODES"
  exit 2
fi

echo "V8_AUTO14_MATCHED_CONTROL_START"
set +e
python3 test/mfamg_replay_v8_candidate_compare.py \
  --exe "$(cygpath -m "$EXE")" \
  --matrix "$(cygpath -m "$A")" \
  --rhs "$(cygpath -m "$B")" \
  --baseline "$(cygpath -m "$BASE")" \
  --generated "$(cygpath -m "$AUTO")" \
  --offline-control "$(cygpath -m "$OFFLINE")" \
  --out-dir "$(cygpath -m "$OUT")" \
  --repetitions 1 --require-offline-control-converged
RC=$?
set -e
echo "V8_AUTO14_COMPARE_EXIT=$RC"
if [[ -f "$OUT/assessment.json" ]]; then
  echo "=== LOCAL AGGREGATE 12 / 14 / 18 ASSESSMENT ==="
  cat "$OUT/assessment.json"
else
  echo "V8_AUTO14_ASSESSMENT=INCONCLUSIVE_NO_ASSESSMENT"
fi
exit "$RC"
