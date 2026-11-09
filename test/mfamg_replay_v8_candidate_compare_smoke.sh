#!/usr/bin/env bash
# Public synthetic only: independent 24-row modes and mocked replay.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
echo mock > "$tmp/A.mtx"
echo mock > "$tmp/b.mtx"
python3 - "$tmp" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
n = 24
cols = [[float(i == j) for i in range(n)] for j in range(18)]
def write(name, vectors):
    text = "%%MatrixMarket matrix array real general\n"
    text += "%% synthetic CI only\n%d %d\n" % (n, len(vectors))
    text += "".join("%.17g\n" % x for col in vectors for x in col)
    (p / name).write_text(text, encoding="ascii")
write("baseline.mtx", cols[:12])
write("generated.mtx", cols[:14])
write("offline.mtx", cols[:18])
bad = [[float(i == 23) for i in range(n)]] + cols[1:14]
write("bad_generated.mtx", bad)
PY
cat > "$tmp/replay.py" <<'PY'
#!/usr/bin/env python3
import os, sys
modefile = os.path.basename(sys.argv[3])
if os.environ.get("LIS_FORBIDDEN_AMBIENT"):
    sys.exit(2)
isgen = modefile == "generated.mtx"
isoff = modefile == "offline.mtx"
modes = 18 if isoff else 14 if isgen else 12
print("REPLAY n=24 nnz=24 modes=%d restart=180 maxiter=3000 tol=1e-5 method=fgcrodr recycle_k=8" % modes)
if os.environ.get("MOCK_MISSING") == "1" and isgen:
    sys.exit(1)
fail = (not isgen and not isoff) or (
    isgen and os.environ.get("MOCK_FAIL_GENERATED") == "1") or (
    isoff and os.environ.get("MOCK_FAIL_OFFLINE") == "1")
print("REPLAY_RECYCLE_STATS builds=1 updates=3 uses=1 rejects=0")
print("REPLAY_RECYCLE_CHECK rc=0 orth_error=1e-12 relation_error=2e-12")
print("REPLAY_RESULT solve_rc=0 status=%d iter=%d reported_residual=8e-6 true_relative=%.9g wall_seconds=%g" %
      (4 if fail else 0, 3001 if fail else 120, 2e-5 if fail else 8e-6,
       9 if isgen else 12 if isoff else 18))
sys.exit(1 if fail else 0)
PY
chmod +x "$tmp/replay.py"
call() {
  python3 "$here/mfamg_replay_v8_candidate_compare.py" \
    --exe "$tmp/replay.py" --matrix "$tmp/A.mtx" --rhs "$tmp/b.mtx" \
    --baseline "$tmp/baseline.mtx" \
    --generated "${2:-$tmp/generated.mtx}" \
    --offline-control "$tmp/offline.mtx" \
    --out-dir "$1" --repetitions 2 --require-offline-control-converged
}
LIS_FORBIDDEN_AMBIENT=1 call "$tmp/success" > "$tmp/success.out"
grep -q '^V8_CANDIDATE_ASSESSMENT=CONTROL_CONVERGED$' "$tmp/success.out"
python3 - "$tmp/success" <<'PY'
import csv, json, pathlib, sys
p = pathlib.Path(sys.argv[1])
with (p / "candidate_runs.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
assert [r["variant"] for r in rows] == [
    "baseline12", "generated14", "offline18",
    "offline18", "generated14", "baseline12"]
assert [int(r["converged"]) for r in rows] == [0, 1, 1, 1, 1, 0]
report = json.loads((p / "assessment.json").read_text())
assert report["assessment"] == "CONTROL_CONVERGED"
assert report["generated_vs_offline_percent_median_wall_reduction"] == 25
assert report["metrics"]["baseline12"]["median_wall_if_all_converged"] is None
assert report["metrics"]["generated14"]["all_converged"]
manifest = json.loads((p / "manifest.json").read_text())
assert set(manifest["input_fingerprints"]) == {
    "exe", "matrix", "rhs", "baseline12", "generated14", "offline18"}
assert len(list(p.glob("*.log"))) == 6
PY
if call "$tmp/success" > "$tmp/overwrite.out" 2>&1; then
  echo "ERROR: output overwrite accepted" >&2
  exit 1
fi
grep -q 'output directory already exists' "$tmp/overwrite.out"
if call "$tmp/bad" "$tmp/bad_generated.mtx" > "$tmp/bad.out" 2>&1; then
  echo "ERROR: altered baseline accepted" >&2
  exit 1
fi
grep -q 'first 12 columns must equal baseline' "$tmp/bad.out"
test ! -e "$tmp/bad"
MOCK_FAIL_GENERATED=1 call "$tmp/nonconverged" > "$tmp/nonconverged.out"
python3 - "$tmp/nonconverged/assessment.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
assert report["metrics"]["generated14"]["converged"] == 0
assert report["generated_vs_offline_percent_median_wall_reduction"] is None
PY
if MOCK_FAIL_OFFLINE=1 call "$tmp/offline_fail" > "$tmp/offline_fail.out" 2>&1; then
  echo "ERROR: offline nonconvergence was accepted" >&2
  exit 1
fi
grep -q 'V8_CANDIDATE_ASSESSMENT=INCONCLUSIVE_OFFLINE_CONTROL' "$tmp/offline_fail.out"
if MOCK_MISSING=1 call "$tmp/missing" > "$tmp/missing.out" 2>&1; then
  echo "ERROR: missing replay result was accepted" >&2
  exit 1
fi
grep -q 'V8_CANDIDATE_ERROR=expected exactly one' "$tmp/missing.out"
test -s "$tmp/missing/candidate_runs.csv"
test -s "$tmp/missing/rep01_run2_generated14.log"
echo V8_CANDIDATE_COMPARISON_SMOKE_OK
