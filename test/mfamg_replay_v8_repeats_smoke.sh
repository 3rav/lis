#!/usr/bin/env bash
# Fake executable only: no contact snapshots or private matrices.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
for file in A04.mtx b04.mtx; do
  echo mock > "$tmp/$file"
done
# Real Matrix Market arrays, not dummy files: 17 must equal all 18
# after deleting the FIRST additional Ritz direction (column 13).
python3 - "$tmp" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
n = 24
vectors = [[float(i == j) for i in range(n)] for j in range(18)]
def write(path, cols):
    lines = ["%%MatrixMarket matrix array real general", f"{n} {len(cols)}"]
    lines.extend(str(x) for col in cols for x in col)
    (root / path).write_text("\n".join(lines) + "\n", encoding="ascii")
write("modes_omit_ritz_01.mtx", vectors[:12] + vectors[13:])
write("modes_all_ritz.mtx", vectors)
write("baseline12.mtx", vectors[:12])
write("bad_omit06.mtx", vectors[:-1])
PY
cat > "$tmp/replay.py" <<'PY'
#!/usr/bin/env python3
import os, sys
if "LIS_STALE_EXPERIMENT_FLAG" in os.environ:
    print("ERROR: inherited ambient LIS option was not sanitized")
    sys.exit(2)
is17 = "omit_ritz_01" in sys.argv[3]
modes = (18 if os.environ.get("MOCK_BAD_MODES") == "1" else 17) if is17 else 18
header = "REPLAY n=24 nnz=8 modes=%d restart=180 maxiter=3000 tol=1e-5 method=fgcrodr recycle_k=8" % modes
if os.environ.get("MOCK_BAD_TOL") == "1":
    header = header.replace("tol=1e-5", "tol=2e-5")
print(header)
if os.environ.get("MOCK_MUTATE") == "1" and is17:
    with open(sys.argv[3], "a", encoding="ascii") as f:
        f.write("% changed mid-run\n")
if (os.environ.get("MOCK_MISSING") == "1" or
        (os.environ.get("MOCK_MISSING_18") == "1" and not is17)):
    sys.exit(1)
fail = os.environ.get("MOCK_MAXITER") == "1" and is17
uses = 0 if os.environ.get("MOCK_NO_REUSE") == "1" and is17 else 1
print("REPLAY_RECYCLE_STATS builds=1 updates=3 uses=%d rejects=0" % uses)
print("REPLAY_RECYCLE_CHECK rc=0 orth_error=1e-12 relation_error=2e-12")
print("REPLAY_RESULT solve_rc=0 status=%d iter=%d reported_residual=8e-6 true_relative=%g wall_seconds=%g" %
      (4 if fail else 0, 2833 if is17 else 2805, 3e-5 if fail else 8e-6, 21.0 if is17 else 24.0))
sys.exit(1 if fail else 0)
PY
chmod +x "$tmp/replay.py"
call() {
  python3 "$here/mfamg_replay_v8_repeats.py" \
    --exe "$tmp/replay.py" --matrix "$tmp/A04.mtx" --rhs "$tmp/b04.mtx" \
    --modes-17 "${2:-$tmp/modes_omit_ritz_01.mtx}" --modes-18 "$tmp/modes_all_ritz.mtx" \
    --baseline-modes "$tmp/baseline12.mtx" \
    --pairs 3 --out-dir "$1" --require-converged
}
LIS_STALE_EXPERIMENT_FLAG=1 call "$tmp/ok" > "$tmp/ok.out"
grep -q 'V8_REPEAT_ASSESSMENT=MATCHED_CONVERGED' "$tmp/ok.out"
python3 - "$tmp/ok" <<'PY'
import csv, json, pathlib, sys
path = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader((path / "matched_runs.csv").open()))
assert [r["variant"] for r in rows] == [
    "omit01_17","all18","all18","omit01_17","omit01_17","all18"]
assert len(list(path.glob("*.log"))) == 6
assert all(r["converged"] == "1" for r in rows)
data = json.loads((path / "assessment.json").read_text())
assert data["omit01_vs_all18_percent_median_wall_reduction"] == 12.5
assert data["paired_median_wall_reduction_percent"] == 12.5
assert data["paired_wall_reduction_percent_by_pair"] == [12.5] * 3
manifest = json.loads((path / "manifest.json").read_text())
assert manifest["mode_preflight"] == "omit01_verified_against_rotated_all18"
assert manifest["matrix_rows_from_modes"] == 24
assert set(manifest["input_fingerprints"]) == {
    "exe", "matrix", "rhs", "modes_17", "modes_18", "baseline_modes"}
assert all(len(info["sha256"]) == 64 and info["bytes"] > 0
           for info in manifest["input_fingerprints"].values())
assert data["manifest_file"] == "manifest.json"
assert all(r["modes"] == ("17" if r["variant"] == "omit01_17" else "18") for r in rows)
PY
if call "$tmp/ok" > "$tmp/dup.out" 2>&1; then exit 1; fi
grep -q 'output directory already exists' "$tmp/dup.out"
if MOCK_MAXITER=1 call "$tmp/maxiter" > "$tmp/max.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE' "$tmp/max.out"
test -s "$tmp/maxiter/matched_runs.csv"
python3 - "$tmp/maxiter/assessment.json" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
assert data["medians"]["all18"]["median_wall_seconds"] is None
assert data["medians"]["omit01_17"]["median_wall_seconds"] is None
assert data["paired_median_wall_reduction_percent"] is None
PY
if MOCK_MISSING=1 call "$tmp/missing" > "$tmp/missing.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ERROR=expected exactly one' "$tmp/missing.out"
if MOCK_NO_REUSE=1 call "$tmp/inactive" > "$tmp/inactive.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE' "$tmp/inactive.out"
if MOCK_BAD_MODES=1 call "$tmp/bad_modes" > "$tmp/bad_modes.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ERROR=unexpected replay configuration' "$tmp/bad_modes.out"
if MOCK_MISSING_18=1 call "$tmp/partial" > "$tmp/partial.out" 2>&1; then exit 1; fi
python3 - "$tmp/partial/matched_runs.csv" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 1 and rows[0]["variant"] == "omit01_17"
PY
test -s "$tmp/partial/pair01_run2_all18.log"
if call "$tmp/bad_basis" "$tmp/bad_omit06.mtx" > "$tmp/bad_basis.out" 2>&1; then
  echo "ERROR: wrong Ritz deletion was accepted" >&2
  exit 1
fi
grep -q '17-mode input is not rotated 18-mode Ritz basis' "$tmp/bad_basis.out"
test ! -e "$tmp/bad_basis"
if MOCK_MUTATE=1 call "$tmp/mutated" > "$tmp/mutated.out" 2>&1; then
  echo "ERROR: input mutation was accepted" >&2
  exit 1
fi
grep -q 'input/EXE changed during benchmark' "$tmp/mutated.out"
test -s "$tmp/mutated/manifest.json"
test -s "$tmp/mutated/pair01_run1_omit01_17.log"
if MOCK_BAD_TOL=1 call "$tmp/wrong_tolerance" > "$tmp/wrong_tolerance.out" 2>&1; then
  echo "ERROR: mismatched replay tolerance was accepted" >&2
  exit 1
fi
grep -q 'V8_REPEAT_ERROR=unexpected replay tolerance' "$tmp/wrong_tolerance.out"
test -s "$tmp/wrong_tolerance/matched_runs.csv"
echo V8_MATCHED_REPEATS_SMOKE_OK
