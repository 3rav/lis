#!/usr/bin/env bash
# LOCAL sequence orchestrator mock regression: public synthetic data only.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
python3 - "$tmp" <<'PY'
from pathlib import Path
import sys
d = Path(sys.argv[1])
for name in ("A1", "A2", "A3", "b1", "b2", "b3"):
    (d / (name + ".mtx")).write_text("synthetic " + name + "\n", encoding="ascii")
(d / "steps.csv").write_text("step,matrix,rhs\none,A1.mtx,b1.mtx\n"
                            "two,A2.mtx,b2.mtx\nthree,A3.mtx,b3.mtx\n",
                            encoding="ascii")
cols = [[float(i == j) for i in range(24)] for j in range(12)]
(d / "modes.mtx").write_text(
    "%%MatrixMarket matrix array real general\n24 12\n"
    + "".join("%g\n" % v for c in cols for v in c),
    encoding="ascii")
PY
cat > "$tmp/replay.py" <<'PY'
#!/usr/bin/env python3
import os
from pathlib import Path
import sys
args = sys.argv[1:]
if os.environ.get("LIS_UNEXPECTED"):
    print("ERROR inherited unrelated LIS variable")
    sys.exit(3)
if args[0] == "--compare-csr":
    left = Path(args[1]).stem
    right = Path(args[2]).stem
    if not Path(args[2]).is_file():
        print("missing input")
        sys.exit(2)
    if os.environ.get("MOCK_BAD_PREFLIGHT") and right == "A2":
        print("unexpected binary failure, not graph mismatch")
        sys.exit(2)
    if left == right:
        print("REPLAY_NEWTON_ERROR replace-values requires a genuinely changed operator")
        print("REPLAY_CSR_COMPARE compatible=0")
        sys.exit(1)
    if right == "A3":
        if os.environ.get("MOCK_GRAPH3_SAME"):
            print("REPLAY_CSR_COMPARE compatible=1 n=24 nnz=70 changed_nnz=5")
            sys.exit(0)
        print("REPLAY_NEWTON_ERROR replace-values CSR column graph mismatch offset=4")
        print("REPLAY_CSR_COMPARE compatible=0")
        sys.exit(1)
    print("REPLAY_CSR_COMPARE compatible=1 n=24 nnz=70 changed_nnz=4")
    sys.exit(0)
first = Path(args[0]).stem
second = os.environ.get("LIS_MFAMG_REPLAY_NEXT_A")
print("REPLAY n=24 nnz=70 modes=12 restart=180 maxiter=3000 tol=1e-5 method=fgcrodr recycle_k=8")
print("REPLAY_RECYCLE_CHECK rc=0 orth_error=1e-13 relation_error=1e-13")
print("REPLAY_RESULT solve_rc=0 status=0 iter=8 reported_residual=4e-6 true_relative=4e-6 wall_seconds=2.0")
if not second:
    sys.exit(0)
next_name = Path(second).stem
replace = os.environ.get("LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES") == "1"
if replace:
    if next_name != "A2":
        sys.exit(9)
    changed = 5 if os.environ.get("MOCK_BAD_CHANGED") else 4
    print("REPLAY_NEWTON_MATRIX mode=inplace_replace_values changed_nnz=%d" % changed)
    print("REPLAY_NEWTON_POLICY cross_uses_delta=1 keeps_delta=1 transports_delta=0 rejects_delta=0")
else:
    if next_name != "A3":
        sys.exit(9)
    print("REPLAY_NEWTON_MATRIX mode=separate_object")
    print("REPLAY_NEWTON_POLICY cross_uses_delta=0 keeps_delta=0 transports_delta=0 rejects_delta=1")
print("REPLAY_NEWTON_STEP step=2 n=24 nnz=70 modes=12 method=fgcrodr recycle_k=8")
print("REPLAY_NEWTON_RECYCLE_CHECK rc=0 orth_error=1e-12 relation_error=2e-12")
fail = bool(os.environ.get("MOCK_FAIL_SECOND") and next_name == "A3")
print("REPLAY_NEWTON_RESULT step=2 solve_rc=0 status=%d iter=%d "
      "reported_residual=8e-6 true_relative=%s wall_seconds=2.4" %
      (4 if fail else 0, 3001 if fail else 9, "3e-5" if fail else "8e-6"))
if os.environ.get("MOCK_MUTATE") and next_name == "A2":
    with Path(second).open("a", encoding="ascii") as out:
        out.write("unsafe changed\n")
third = os.environ.get("LIS_MFAMG_REPLAY_THIRD_A")
if third:
    assert Path(third).stem == "A3"
    same = os.environ.get("LIS_MFAMG_REPLAY_THIRD_REPLACE_VALUES") == "1"
    mode3 = "inplace_replace_values" if same else "separate_object"
    changed3 = (6 if os.environ.get("MOCK_BAD_3_CHANGED") else 5)
    print("REPLAY_SEQUENCE_MATRIX step=3 mode=%s%s" %
          (mode3, " changed_nnz=%d" % changed3 if same else ""))
    print("REPLAY_SEQUENCE_STEP step=3 n=24 nnz=70 modes=12 method=fgcrodr recycle_k=8")
    cross = int(same or bool(os.environ.get("MOCK_STALE_THIRD")))
    reject = 0 if cross else 1
    print("REPLAY_SEQUENCE_POLICY step=3 cross_uses_delta=%d keeps_delta=%d transports_delta=0 rejects_delta=%d" %
          (cross, cross, reject))
    print("REPLAY_SEQUENCE_RECYCLE_CHECK step=3 rc=0 orth_error=1e-12 relation_error=2e-12")
    fail3 = bool(os.environ.get("MOCK_FAIL_THIRD"))
    print("REPLAY_SEQUENCE_RESULT step=3 solve_rc=0 status=%d iter=%d "
          "reported_residual=8e-6 true_relative=%s wall_seconds=3.1" %
          (4 if fail3 else 0, 3001 if fail3 else 12,
           "3e-5" if fail3 else "7e-6"))
    if fail3:
        sys.exit(1)
sys.exit(1 if fail else 0)
PY
chmod +x "$tmp/replay.py"
call() {
  python3 "$here/mfamg_v8_snapshot_sequence.py" \
    --exe "$tmp/replay.py" --modes "$tmp/modes.mtx" \
    --snapshots "$tmp/steps.csv" --out-dir "$1" --fresh-control
}

LIS_UNEXPECTED=1 call "$tmp/success" > "$tmp/success.out"
grep -q '^V8_SEQUENCE_ASSESSMENT=PAIRWISE_PASS' "$tmp/success.out"
python3 - "$tmp/success" <<'PY'
import csv
import json
from pathlib import Path
import sys
d=Path(sys.argv[1])
summary=json.loads((d/"summary.json").read_text())
assert summary["assessment"] == "PAIRWISE_PASS"
assert summary["completed_pairs"] == 2
assert summary["passed_pairs"] == 2
manifest=json.loads((d/"manifest.json").read_text())
assert manifest["step_ids"] == ["one","two","three"]
assert manifest["n_from_modes"] == 24
assert manifest["warning"].startswith("Fresh solver per adjacent PAIR")
with (d/"transitions.csv").open(newline="") as f:
    rows=list(csv.DictReader(f))
assert [r["graph"] for r in rows] == ["same_csr","changed_graph"]
assert [r["mode"] for r in rows] == ["inplace_replace_values","separate_object"]
assert [int(r["cross_uses_delta"]) for r in rows] == [1,0]
assert [int(r["rejects_delta"]) for r in rows] == [0,1]
assert all(int(r["pair_pass"]) == 1 for r in rows)
assert all(float(r["second_true_relative"]) <= 1e-5 for r in rows)
assert all(float(r["fresh_true_relative"]) <= 1e-5 for r in rows)
assert len(list(d.glob("*.log"))) == 6
PY
if call "$tmp/success" > "$tmp/overwrite.out" 2>&1; then
  echo "ERROR: overwrite accepted" >&2; exit 1
fi
grep -q 'output directory already exists' "$tmp/overwrite.out"

if MOCK_FAIL_SECOND=1 call "$tmp/maxiter" > "$tmp/maxiter.out" 2>&1; then
  echo "ERROR: MAXITER accepted as successful sequence" >&2; exit 1
fi
python3 - "$tmp/maxiter" <<'PY'
import json,sys
from pathlib import Path
d=Path(sys.argv[1]); r=json.loads((d/"summary.json").read_text())
assert r["assessment"]=="INCONCLUSIVE_OR_FAILED"
assert r["completed_pairs"]==2 and r["passed_pairs"]==1
PY

if MOCK_BAD_PREFLIGHT=1 call "$tmp/preflight_error" > "$tmp/preflight_error.out" 2>&1; then
  echo "ERROR: bad preflight masked as graph mismatch" >&2; exit 1
fi
python3 - "$tmp/preflight_error" <<'PY'
import csv,json,sys
from pathlib import Path
d=Path(sys.argv[1])
with (d/"transitions.csv").open(newline="") as f:r=list(csv.DictReader(f))
assert len(r)==1
assert "expected exactly one REPLAY_CSR_COMPARE" in r[0]["error"]
assert not (d/"pair01_replay.log").exists()
PY

if python3 "$here/mfamg_v8_snapshot_sequence.py" \
  --exe "$tmp/replay.py" --modes "$tmp/modes.mtx" \
  --snapshots "$tmp/steps.csv" --out-dir "$tmp/invalid_tol" \
  --tol 1e-4 > "$tmp/invalid.out" 2>&1; then
  echo "ERROR: incompatible tolerance accepted" >&2; exit 1
fi
test ! -e "$tmp/invalid_tol"

# In each triplet the same mocked process provides three sequential
# step diagnostics (the separate baseline runner only keeps two).
triplet() {
  python3 "$here/mfamg_v8_triplet_windows.py" \
    --exe "$tmp/replay.py" --modes "$tmp/modes.mtx" \
    --snapshots "$tmp/steps.csv" --out-dir "$1" --fresh-third
}
triplet "$tmp/triplet" > "$tmp/triplet.out"
grep -q '^V8_TRIPLET_ASSESSMENT=PERSISTENT_TRIPLETS_PASS' "$tmp/triplet.out"
python3 - "$tmp/triplet" <<'PY'
import csv,json,sys
from pathlib import Path
d=Path(sys.argv[1]); r=json.loads((d/"summary.json").read_text())
assert r["assessment"]=="PERSISTENT_TRIPLETS_PASS"
assert r["completed_windows"]==1 and r["passed_windows"]==1
with (d/"triplets.csv").open(newline="") as f: rows=list(csv.DictReader(f))
assert len(rows)==1
x=rows[0]
assert x["graph12"]=="same_csr" and x["graph23"]=="changed_graph"
assert x["mode2"]=="inplace_replace_values"
assert x["mode3"]=="separate_object"
assert int(x["third_cross_uses_delta"])==0
assert int(x["third_rejects_delta"])==1
assert float(x["third_true_relative"]) <= 1e-5
assert x["pass"]=="1"
assert len(list(d.glob("*.log")))==4
PY
MOCK_GRAPH3_SAME=1 triplet "$tmp/triplet_keep" > "$tmp/triplet_keep.out"
python3 - "$tmp/triplet_keep" <<'PY'
import csv,sys
from pathlib import Path
with (Path(sys.argv[1])/"triplets.csv").open(newline="") as f:
    rows=list(csv.DictReader(f))
assert rows[0]["mode3"]=="inplace_replace_values"
assert rows[0]["graph23"]=="same_csr"
assert rows[0]["changed23"]=="5"
assert rows[0]["third_cross_uses_delta"]=="1"
assert rows[0]["pass"]=="1"
PY
if MOCK_FAIL_THIRD=1 triplet "$tmp/triplet_maxiter" > "$tmp/triplet_maxiter.out" 2>&1; then
  echo "ERROR: third-step MAXITER accepted" >&2; exit 1
fi
if MOCK_STALE_THIRD=1 triplet "$tmp/triplet_stale" > "$tmp/triplet_stale.out" 2>&1; then
  echo "ERROR: stale third-step cross-graph recycling accepted" >&2; exit 1
fi
grep -q '^V8_TRIPLET_ASSESSMENT=INCONCLUSIVE_OR_FAILED' "$tmp/triplet_stale.out"
if MOCK_GRAPH3_SAME=1 MOCK_BAD_3_CHANGED=1 \
  triplet "$tmp/triplet_bad_nnz" > "$tmp/triplet_bad_nnz.out" 2>&1; then
  echo "ERROR: third-step changed CSR count mismatch accepted" >&2; exit 1
fi
grep -q '^V8_TRIPLET_ASSESSMENT=INCONCLUSIVE_OR_FAILED' "$tmp/triplet_bad_nnz.out"

if MOCK_MUTATE=1 call "$tmp/mutated" > "$tmp/mutated.out" 2>&1; then
  echo "ERROR: modified input accepted" >&2; exit 1
fi
python3 - "$tmp/mutated" <<'PY'
import csv,sys
from pathlib import Path
p=Path(sys.argv[1])
with (p/"transitions.csv").open(newline="") as f:r=list(csv.DictReader(f))
assert len(r)==1
assert "input/EXE changed during pair:" in r[0]["error"]
assert (p/"pair01_replay.log").is_file()
PY
if MOCK_BAD_CHANGED=1 call "$tmp/bad_changed" > "$tmp/bad_changed.out" 2>&1; then
  echo "ERROR: mismatch of CSR preflight and replay changed_nnz accepted" >&2; exit 1
fi
python3 - "$tmp/bad_changed" <<'PY'
import csv,sys
from pathlib import Path
with (Path(sys.argv[1])/"transitions.csv").open(newline="") as f:
    rows=list(csv.DictReader(f))
assert len(rows)==1
assert "changed_nnz disagrees with CSR preflight" in rows[0]["error"]
PY
echo MFAMG_V8_SEQUENCE_LOCAL_ORCHESTRATOR_SMOKE_OK
