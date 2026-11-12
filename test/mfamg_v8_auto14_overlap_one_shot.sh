#!/usr/bin/env bash
# Optional read-only post-mortem of an already completed private A04 12/14/18 run.
# Uses the six offline Ritz directions as a DIAGNOSTIC ONLY, not for generating modes.
# No new solver runs, no changes to input data or local checkout, no upload.
set -euo pipefail
cd "$HOME/lis-v8"
ROOT="$HOME/mfamg_replay"
BRANCH="experiment/mfamg-fgcrodr-replay-v8-20261108"
REF="refs/remotes/origin/$BRANCH"
BASE="$ROOT/challenge_31/modes04.mtx"
GENERATED="$ROOT/v8_auto14_small_singular_upstream214_20261110.mtx"
OFFLINE="$ROOT/v8_mode_audit/ritz/modes_all_ritz.mtx"
MANIFEST="$ROOT/v8_auto14_compare_upstream214_20261110/manifest.json"
EXPECTED_BLOB="815a73ada1d18d970ec9b6d24306d23ff631fb35"

if [[ "$(printenv MSYSTEM || true)" != UCRT64 ]]; then
  echo "V8_AUTO14_OVERLAP=BLOCKED_REQUIRES_UCRT64"
  exit 2
fi
for command in git tar python3 cygpath mktemp; do
  command -v "$command" >/dev/null || { echo "V8_AUTO14_OVERLAP=BLOCKED_MISSING_TOOL"; exit 2; }
done
for file in "$BASE" "$GENERATED" "$OFFLINE" "$MANIFEST"; do
  if [[ ! -s "$file" ]]; then
    echo "V8_AUTO14_OVERLAP=BLOCKED_MISSING_PRIOR_RUN_INPUT"
    exit 2
  fi
done
# Require the exact three mode files fingerprinted by the completed benchmark.
python3 - "$(cygpath -m "$MANIFEST")" "$(cygpath -m "$BASE")" \
    "$(cygpath -m "$GENERATED")" "$(cygpath -m "$OFFLINE")" <<'PY'
import hashlib, json, pathlib, sys
try:
    manifest = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
    if manifest.get("schema") != "mfamg-v8-automatic-candidate-control-v1":
        raise ValueError("unexpected manifest schema")
    inputs = manifest["input_fingerprints"]
    for name, file in zip(("baseline12", "generated14", "offline18"), sys.argv[2:]):
        path = pathlib.Path(file)
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                digest.update(block)
        expected = inputs[name]
        if digest.hexdigest() != expected["sha256"] or path.stat().st_size != expected["bytes"]:
            raise ValueError("a mode input differs from completed replay")
except (OSError, ValueError, KeyError, TypeError, UnicodeError):
    print("V8_AUTO14_OVERLAP=BLOCKED_PRIOR_RUN_FINGERPRINT_MISMATCH")
    sys.exit(2)
print("V8_AUTO14_OVERLAP_PRIOR_RUN_MODES_UNCHANGED=1")
PY

git fetch --no-tags origin "refs/heads/$BRANCH:$REF"
if [[ "$(git show "$REF:src/solver/lis_solver_fgcrodr.c" | git hash-object --stdin)" != "$EXPECTED_BLOB" ]]; then
  echo "V8_AUTO14_OVERLAP=BLOCKED_UPSTREAM_SOURCE_DIVERGENCE"
  exit 2
fi
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT
umask 077
git archive "$REF" \
  test/mfamg_mode_audit.py \
  test/mfamg_mode_ritz_selection.py \
  test/mfamg_v8_auto14_overlap.py | tar -xf - -C "$TMP"
python3 "$TMP/test/mfamg_v8_auto14_overlap.py" \
  --baseline "$(cygpath -m "$BASE")" \
  --generated "$(cygpath -m "$GENERATED")" \
  --offline-control "$(cygpath -m "$OFFLINE")"
