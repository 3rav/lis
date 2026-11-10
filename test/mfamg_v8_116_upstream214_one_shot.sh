#!/usr/bin/env bash
# PRIVATE 116-spring fixed-A acceptance screen, using the ALREADY BUILT
# accepted upstream develop #214 FGCRO-DR binary from the winfix1 checkout.
# Historical 116 control: A02/b02/modes02, 12 modes, MFAMG factor=2,
# restart=180, MAXITER=1000. Match FGMRES and FGCRODR environments.
# NO existing files overwritten, NO Git checkout/reset/clean/merge,
# NO automatic data upload and NO production CalculiX modification.
set -euo pipefail
cd "$HOME/lis-v8"

ROOT="$HOME/mfamg_replay"
BRANCH="experiment/mfamg-fgcrodr-replay-v8-20261108"
REF="refs/remotes/origin/$BRANCH"
EXE="$ROOT/build-v8-upstream214-20261109-winfix1/test/.libs/mfamgreplay.exe"
SOLVER="$ROOT/src-v8-upstream214-20261109-winfix1/src/solver/lis_solver_fgcrodr.c"
EXPECTED_BLOB="815a73ada1d18d970ec9b6d24306d23ff631fb35"
A="$ROOT/control_116/A02.mtx"
B="$ROOT/control_116/b02.mtx"
MODES="$ROOT/control_116/modes02.mtx"
SCRIPTS="$ROOT/v8_116_python_upstream214_20261110"
OUT="$ROOT/v8_116_fgmres_vs_fgcrodr_upstream214_20261110"

if [[ "$(printenv MSYSTEM || true)" != "UCRT64" ]]; then
  echo "STOP: wymagany terminal MSYS2 UCRT64"
  exit 2
fi
for cmd in git tar python3 cygpath; do
  command -v "$cmd" >/dev/null || { echo "STOP: brak $cmd"; exit 2; }
done
for file in "$EXE" "$SOLVER" "$A" "$B" "$MODES"; do
  if [[ ! -s "$file" ]]; then
    echo "STOP: brak niepustego EXE lub prywatnego pliku: $file"
    exit 2
  fi
done
for target in "$SCRIPTS" "$OUT"; do
  if [[ -e "$target" ]]; then
    echo "STOP: katalog istnieje, bez nadpisywania: $target"
    exit 2
  fi
done
if [[ "$(git hash-object "$SOLVER")" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: lokalne źródło solvera nie jest zaakceptowanym upstream #214"
  exit 2
fi

git fetch --no-tags origin "refs/heads/$BRANCH:$REF"
if [[ "$(git show "$REF:src/solver/lis_solver_fgcrodr.c" | git hash-object --stdin)" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: zdalny solver różni się od zaakceptowanego upstream #214"
  exit 2
fi

# Archive only 3 public Python scripts. The selected paths have no symlinks;
# neither matrix values nor mode columns are copied from the private inputs.
umask 077
mkdir "$SCRIPTS"
git archive "$REF" \
  test/mfamg_v8_solver_pair_compare.py \
  test/mfamg_replay_v8_repeats.py \
  test/mfamg_mode_audit.py | tar -xf - -C "$SCRIPTS"

if ! grep -q 'add_argument("--fine-agg-factor"' \
  "$SCRIPTS/test/mfamg_v8_solver_pair_compare.py"; then
  echo "STOP: runner nie ma jawnego ustawienia factor=2"
  exit 2
fi

echo "V8_116_UPSTREAM214_PAIR_START"
cd "$SCRIPTS"
set +e
python3 test/mfamg_v8_solver_pair_compare.py \
  --exe-v8 "$(cygpath -m "$EXE")" \
  --matrix "$(cygpath -m "$A")" \
  --rhs "$(cygpath -m "$B")" \
  --modes "$(cygpath -m "$MODES")" \
  --expected-modes 12 \
  --fine-agg-factor 2 \
  --maxiter 1000 \
  --pairs 3 --require-both-converged \
  --out-dir "$(cygpath -m "$OUT")"
RC=$?
set -e
echo "V8_116_BENCHMARK_EXIT=$RC"
if [[ -f "$OUT/assessment.json" ]]; then
  echo "=== 116 PRIVATE AGGREGATE ==="
  cat "$OUT/assessment.json"
else
  echo "Brak assessment.json; pozostawiono lokalne dzienniki w $OUT"
fi
exit "$RC"
