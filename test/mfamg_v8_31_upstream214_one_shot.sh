#!/usr/bin/env bash
# One guarded, PRIVATE UCRT64 replay of v8 FGMRES vs v8 FGCRODR (A04/b04, 18 Ritz).
# This script does NOT modify the user's checkout/worktree or upload contact data.
# It archives the current remote Draft branch to a NEW private source folder,
# builds into a NEW private out-of-tree folder, and writes a NEW local result folder.
# It will never delete/overwrite pre-existing output/source/build directories.
set -euo pipefail

cd "$HOME/lis-v8"
ROOT="$HOME/mfamg_replay"
BRANCH="experiment/mfamg-fgcrodr-replay-v8-20261108"
REF="refs/remotes/origin/$BRANCH"
EXPECTED_BLOB="815a73ada1d18d970ec9b6d24306d23ff631fb35"
SRC="$ROOT/src-v8-upstream214-20261109"
BUILD="$ROOT/build-v8-upstream214-20261109"
OUT="$ROOT/v8_fgmres_vs_fgcrodr_31_upstream214_20261109"
MATRIX="$ROOT/challenge_31/A04.mtx"
RHS="$ROOT/challenge_31/b04.mtx"
MODES="$ROOT/v8_mode_audit/ritz/modes_all_ritz.mtx"

MSYS_KIND="$(printenv MSYSTEM || true)"
if [[ "$MSYS_KIND" != "UCRT64" ]]; then
  echo "STOP: wymagany terminal MSYS2 UCRT64 (MSYSTEM=$MSYS_KIND)"
  exit 2
fi
for exe in git tar python3 cygpath autoreconf make; do
  command -v "$exe" >/dev/null || {
    echo "STOP: brak narzędzia $exe"
    exit 2
  }
done
for file in "$MATRIX" "$RHS" "$MODES"; do
  if [[ ! -s "$file" ]]; then
    echo "STOP: brak niepustego prywatnego pliku: $file"
    exit 2
  fi
done
for path in "$SRC" "$BUILD" "$OUT"; do
  if [[ -e "$path" ]]; then
    echo "STOP: istniejący katalog nie zostanie nadpisany: $path"
    exit 2
  fi
done
if [[ ! -d "$ROOT" || ! -w "$ROOT" ]]; then
  echo "STOP: katalog prywatny nie istnieje lub nie jest zapisywalny: $ROOT"
  exit 2
fi

# Remote-tracking ref only: never merge, reset, clean or check out locally.
git fetch --no-tags origin "refs/heads/$BRANCH:$REF"
SOURCE_SHA="$(git rev-parse "$REF")"
BLOB="$(git show "$REF:src/solver/lis_solver_fgcrodr.c" | git hash-object --stdin)"
if [[ "$BLOB" != "$EXPECTED_BLOB" ]]; then
  echo "STOP: solver różni się od zaakceptowanego upstream develop #214"
  echo "BLOB=$BLOB EXPECTED=$EXPECTED_BLOB"
  exit 2
fi
echo "SOURCE_SHA=$SOURCE_SHA"
echo "UPSTREAM_SOLVER_BLOB_OK=$BLOB"

# No mkdir -m: Windows ACL-backed MSYS2 locations may refuse chmod.
umask 077
mkdir "$SRC"
git archive "$REF" | tar -xf - -C "$SRC"
mkdir "$BUILD"
(
  cd "$SRC"
  autoreconf -fiv
)
(
  cd "$BUILD"
  "$SRC/configure" \
    --prefix=/ucrt64 \
    --enable-omp --enable-shared \
    --disable-longlong --enable-half \
    --disable-longdouble --enable-quad \
    --disable-complex --disable-mpi \
    --disable-test --with-gnu-ld
  make -j2
  make -C test mfamgreplay.exe
)
EXE="$BUILD/test/.libs/mfamgreplay.exe"
if [[ ! -s "$EXE" ]]; then
  echo "STOP: nowy mfamgreplay.exe nie powstał"
  exit 2
fi
echo "NEW_UPSTREAM214_EXE_OK"

# The user's Python may be native Windows Python, not MSYS Python.
# Use Windows-style absolute paths for every file passed to Python.
cd "$SRC"
set +e
python3 test/mfamg_v8_solver_pair_compare.py \
  --exe-v8 "$(cygpath -m "$EXE")" \
  --matrix "$(cygpath -m "$MATRIX")" \
  --rhs "$(cygpath -m "$RHS")" \
  --modes "$(cygpath -m "$MODES")" \
  --expected-modes 18 --pairs 3 --require-both-converged \
  --out-dir "$(cygpath -m "$OUT")"
RC=$?
set -e
echo "V8_31_SOLVER_BENCHMARK_EXIT=$RC"
if [[ -f "$OUT/assessment.json" ]]; then
  echo "=== PRIVATE ASSESSMENT (aggregate only) ==="
  cat "$OUT/assessment.json"
else
  echo "Brak assessment.json; dzienniki pozostają lokalne, nie usuwaj ich."
fi
exit "$RC"
