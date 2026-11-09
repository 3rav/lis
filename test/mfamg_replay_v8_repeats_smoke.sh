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
if os.environ.get("MOCK_MISSING") == "1":
    print("REPLAY n=4 method=fgcrodr")
    sys.exit(1)
fail = os.environ.get("MOCK_MAXITER") == "1" and is17
print("REPLAY_RECYCLE_STATS builds=1 updates=3 uses=1 rejects=0")
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
PY
if call "$tmp/ok" > "$tmp/dup.out" 2>&1; then exit 1; fi
grep -q 'output directory already exists' "$tmp/dup.out"
if MOCK_MAXITER=1 call "$tmp/maxiter" > "$tmp/max.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE' "$tmp/max.out"
test -s "$tmp/maxiter/matched_runs.csv"
if MOCK_MISSING=1 call "$tmp/missing" > "$tmp/missing.out" 2>&1; then exit 1; fi
grep -q 'V8_REPEAT_ERROR=expected exactly one' "$tmp/missing.out"
echo V8_MATCHED_REPEATS_SMOKE_OK
