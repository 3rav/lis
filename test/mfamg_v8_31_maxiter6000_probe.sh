#!/usr/bin/env bash
# UCRT64 private one-pair convergence probe with the SAME freshly built
# upstream-#214 EXE: FGMRES vs FGCRODR, common MAXITER=6000.
# Does not overwrite any existing src/build/result path or touch git worktree.
# No private data or logs are uploaded.
set -euo pipefail
cd "$HOME/lis-v8"

ROOT="$HOME/mfamg_replay"
BRANCH="experiment/mfamg-fgcrodr-replay-v8-20261108"
REF="refs/remotes/origin/$BRANCH"
EXE="$ROOT/build-v8-upstream214-20261109-winfix1/test/.libs/mfamgreplay.exe"
LOCAL_SOLVER="$ROOT/src-v8-upstream214-20261109-winfix1/src/solver/lis_solver_fgcrodr.c"
EXPECTED_BLOB="815a73ada1d18d970ec9b6d24306d23ff631fb35"
SCRIPTS="$ROOT/v8_pair_6000_python_20261109"
OUT="$ROOT/v8_pair_6000_probe_31_20261109"
MATRIX="$ROOT/challenge_31/A04.mtx"
RHS="$ROOT/challenge_31/b04.mtx"
MODES="$ROOT/v8_mode_audit/ritz/modes_all_ritz.mtx"

if [[ "$(printenv MSYSTEM || true)" != "UCRT64" ]]; then
  echo "STOP: wymagany MSYS2 UCRT64"
  exit 2
fi
for cmd in git tar python3 cygpath; do
  command -v "$cmd" >/dev/null || { echo "STOP: brak $cmd"; exit 2; }
done
for file in "$EXE" "$LOCAL_SOLVER" "$MATRIX" "$RHS" "$MODES"; do
  if [[ ! -s "$file" ]]; then
    echo "STOP: brak prywatnego EXE/danych: $file"
    exit 2
  fi
done
if [[ -e "$SCRIPTS" || -e "$OUT" ]]; then
  echo "STOP: nowy katalog skryptów lub wyników już istnieje; nic nie nadpisano"
  exit 2
fi
if [[ "$(git hash-object "$LOCAL_SOLVER")" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: lokalnie skompilowany solver nie odpowiada upstream develop #214"
  exit 2
fi
git fetch --no-tags origin "refs/heads/$BRANCH:$REF"
if [[ "$(git show "$REF:src/solver/lis_solver_fgcrodr.c" | git hash-object --stdin)" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: obecny zdalny solver nie jest zgodny z zaakceptowanym upstream"
  exit 2
fi
# Pull ONLY public Python scripts: archive contains no symlinks or private data.
umask 077
mkdir "$SCRIPTS"
git archive "$REF" \
  test/mfamg_v8_solver_pair_compare.py \
  test/mfamg_replay_v8_repeats.py \
  test/mfamg_mode_audit.py | tar -xf - -C "$SCRIPTS"
if ! grep -q 'add_argument("--maxiter"' "$SCRIPTS/test/mfamg_v8_solver_pair_compare.py"; then
  echo "STOP: nowy runner nie ma wspólnego limitu MAXITER"
  exit 2
fi
echo "MAXITER6000_SAME_EXE_PROBE_START"
cd "$SCRIPTS"
set +e
python3 test/mfamg_v8_solver_pair_compare.py \
  --exe-v8 "$(cygpath -m "$EXE")" \
  --matrix "$(cygpath -m "$MATRIX")" \
  --rhs "$(cygpath -m "$RHS")" \
  --modes "$(cygpath -m "$MODES")" \
  --expected-modes 18 --maxiter 6000 --pairs 1 \
  --require-both-converged \
  --out-dir "$(cygpath -m "$OUT")"
RC=$?
set -e
echo "MAXITER6000_PROBE_EXIT=$RC"
if [[ -f "$OUT/assessment.json" ]]; then
  echo "=== AGGREGATE ASSESSMENT ==="
  cat "$OUT/assessment.json"
else
  echo "Brak assessment.json; lokalne logi pozostają nienaruszone w $OUT"
fi
exit "$RC"
