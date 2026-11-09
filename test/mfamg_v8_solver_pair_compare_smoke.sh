#!/usr/bin/env bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
{
 echo '%%MatrixMarket matrix array real general'
 echo '24 2'
 for ((j=1;j<=2;j++)); do
   for ((i=1;i<=24;i++)); do
     if (( i==j )); then echo 1; else echo 0; fi
   done
 done
} > "$tmp/modes.mtx"
printf 'fake A\n' > "$tmp/A.mtx"
printf 'fake b\n' > "$tmp/b.mtx"
cat > "$tmp/mock.py" <<'PY'
#!/usr/bin/env python3
import os,sys
from pathlib import Path
method=os.environ.get("LIS_MFAMG_REPLAY_METHOD")
v7=Path(sys.argv[0]).name=="mock_v7.py"
if v7:
    assert method is None
    print("REPLAY n=24 nnz=70 modes=2 restart=180 maxiter=3000 tol=1e-5")
else:
    assert method in ("fgmres","fgcrodr")
    announced="fgmres" if os.environ.get("MOCK_BAD_METHOD")=="1" else method
    print("REPLAY n=24 nnz=70 modes=2 restart=180 maxiter=3000 tol=1e-5 method=%s recycle_k=%d" % (announced,8 if method=="fgcrodr" else 0))
assert os.environ.get("LIS_UNWANTED_FLAG") is None
assert os.environ.get("LIS_MFAMG_REPLAY_MAXITER")=="3000"
if os.environ.get("MOCK_MUTATE")=="1":
    Path(sys.argv[2]).write_text("changed\n")
if method=="fgcrodr":
    uses=0 if os.environ.get("MOCK_NO_REUSE")=="1" else 1
    print("REPLAY_RECYCLE_STATS builds=1 updates=3 uses=%d rejects=0" % uses)
    print("REPLAY_RECYCLE_CHECK rc=0 orth_error=1e-12 relation_error=2e-12")
bad=method=="fgcrodr" and os.environ.get("MOCK_MAXITER")=="1"
print("REPLAY_RESULT solve_rc=0 status=%d iter=%d reported_residual=8e-6 true_relative=%g wall_seconds=%g" % (4 if bad else 0,3000 if bad else (60 if method=="fgcrodr" else 80),4e-5 if bad else 8e-6,24 if method=="fgcrodr" else 30))
sys.exit(1 if bad else 0)
PY
chmod +x "$tmp/mock.py"
cp "$tmp/mock.py" "$tmp/mock_v7.py"
chmod +x "$tmp/mock_v7.py"
call() {
 out="$1"
 shift
 python3 "$here/mfamg_v8_solver_pair_compare.py" \
  --exe-v8 "$tmp/mock.py" --matrix "$tmp/A.mtx" --rhs "$tmp/b.mtx" \
  --modes "$tmp/modes.mtx" --expected-modes 2 --pairs 3 \
  --require-both-converged --out-dir "$out" "$@"
}
LIS_UNWANTED_FLAG=1 call "$tmp/ok" > "$tmp/ok.out"
grep -q 'V8_SOLVER_PAIR_ASSESSMENT=MATCHED_CONVERGED' "$tmp/ok.out"
python3 - "$tmp/ok" <<'PY'
import sys,csv,json
from pathlib import Path
p=Path(sys.argv[1])
rows=list(csv.DictReader((p/"paired_runs.csv").open()))
assert len(rows)==6
assert [r["variant"] for r in rows]==["fgmres_control","fgcrodr_v8","fgcrodr_v8","fgmres_control","fgmres_control","fgcrodr_v8"]
assert all(r["converged"]=="1" for r in rows)
assert len(list(p.glob("*.log")))==6
d=json.loads((p/"assessment.json").read_text())
assert d["baseline_type"]=="same_v8_binary_fgmres"
assert d["median_paired_v8_wall_reduction_percent"]==20
assert d["v8_wall_reduction_percent_from_separate_medians"]==20
assert d["metrics"]["fgcrodr_v8"]["median_iterations"]==60
PY
call "$tmp/v7" --exe-v7 "$tmp/mock_v7.py" > "$tmp/v7.out"
grep -q 'V8_SOLVER_PAIR_ASSESSMENT=MATCHED_CONVERGED' "$tmp/v7.out"
grep -q 'historical_v7_binary' "$tmp/v7/assessment.json"
if call "$tmp/duplicate" --exe-v7 "$tmp/mock.py" > "$tmp/duplicate.out" 2>&1; then exit 1; fi
grep -q 'all paths must refer to distinct files' "$tmp/duplicate.out"
if call "$tmp/ok" > "$tmp/existing.out" 2>&1; then exit 1; fi
grep -q 'output directory already exists' "$tmp/existing.out"
if MOCK_MAXITER=1 call "$tmp/nonconverged" > "$tmp/max.out" 2>&1; then exit 1; fi
grep -q 'V8_SOLVER_PAIR_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE' "$tmp/max.out"
grep -q '"paired_v8_wall_reduction_percent_by_pair": null' "$tmp/nonconverged/assessment.json"
if MOCK_NO_REUSE=1 call "$tmp/inactive" > "$tmp/inactive.out" 2>&1; then exit 1; fi
grep -q 'V8_SOLVER_PAIR_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE' "$tmp/inactive.out"
if MOCK_BAD_METHOD=1 call "$tmp/bad" > "$tmp/bad.out" 2>&1; then exit 1; fi
grep -q 'unexpected replay configuration' "$tmp/bad.out"
if MOCK_MUTATE=1 call "$tmp/mutated" > "$tmp/mutated.out" 2>&1; then exit 1; fi
grep -q 'input/EXE changed during benchmark' "$tmp/mutated.out"
echo V8_FGMRES_FGCRODR_PAIRED_SMOKE_OK
