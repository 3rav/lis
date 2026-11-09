#!/usr/bin/env bash
# Fake executable only: no contact snapshots or private matrices.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
for file in A04.mtx b04.mtx modes_omit_ritz_01.mtx modes_all_ritz.mtx; do
  echo mock > "$tmp/$file"
done
cat > "$tmp/replay.py" <<'PY'
#!/usr/bin/env python3
import os, sys
is17 = "omit_ritz_01" in sys.argv[3]
modes = (18 if os.environ.get("MOCK_BAD_MODES") == "1" else 17) if is17 else 18
print("REPLAY n=4 nnz=8 modes=%d restart=180 maxiter=3000 tol=1e-5 method=fgcrodr recycle_k=8" % modes)
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
    --modes-17 "$tmp/modes_omit_ritz_01.mtx" --modes-18 "$tmp/modes_all_ritz.mtx" \
    --pairs 3 --out-dir "$1" --require-converged
}
call "$tmp/ok" > "$tmp/ok.out"
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
assert all(r["modes"] == ("17" if r["variant"] == "omit01_17" else "18") for r in rows)
PY
if call "$tmp/ok" > "$tmp/dup.out" 2>&1; then exit 1; fi
grep -q 'output directory already exists' "$tmp/dup.out"
if MOCK_MAXITER=1 call "$tmp/maxiter" > "$tmp/max.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE' "$tmp/max.out"
test -s "$tmp/maxiter/matched_runs.csv"
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
echo V8_MATCHED_REPEATS_SMOKE_OK
