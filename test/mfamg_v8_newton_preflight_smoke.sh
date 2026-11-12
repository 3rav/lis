#!/usr/bin/env bash
# Synthetic only: no real contact data, no writes by preflight itself.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
python3 - "$tmp" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
n=24
for kind in ("A1","A2","A3"):
    with (p/(kind+".mtx")).open("w") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write("24 24 24\n")
        for i in range(1,n+1):
            f.write("%d %d %d\n" % (i,i,4+(i%2)))
for kind in ("b1","b2","b3"):
    (p/(kind+".mtx")).write_text(
        "%%MatrixMarket matrix array real general\n24 1\n"+"1\n"*24)
cols=[[float(i==j) for i in range(n)] for j in range(12)]
(p/"modes.mtx").write_text(
    "%%MatrixMarket matrix array real general\n24 12\n"+
    "".join("%g\n"%v for c in cols for v in c))
(p/"steps.csv").write_text(
    "step,matrix,rhs\none,A1.mtx,b1.mtx\ntwo,A2.mtx,b2.mtx\nthree,A3.mtx,b3.mtx\n")
(p/"onlytwo.csv").write_text(
    "step,matrix,rhs\none,A1.mtx,b1.mtx\ntwo,A2.mtx,b2.mtx\n")
(p/"duplicate.csv").write_text(
    "step,matrix,rhs\none,A1.mtx,b1.mtx\ntwo,A1.mtx,b2.mtx\nthree,A3.mtx,b3.mtx\n")
PY
cat > "$tmp/mock.py" <<'PY'
#!/usr/bin/env python3
import os
from pathlib import Path
import sys
if sys.argv[1]!="--compare-csr":
    sys.exit(8)
if os.environ.get("MOCK_MUTATE"):
    with Path(sys.argv[3]).open("a") as f:
        f.write("bad mutation\n")
if os.environ.get("MOCK_MALFORMED"):
    print("bad log with no preflight result")
    sys.exit(2)
left,right=(Path(x).stem for x in sys.argv[2:])
if right=="A2":
    if os.environ.get("MOCK_NO_CHANGE"):
        print("REPLAY_NEWTON_ERROR replace-values requires a genuinely changed operator")
        print("REPLAY_CSR_COMPARE compatible=0")
        sys.exit(1)
    print("REPLAY_CSR_COMPARE compatible=1 n=24 nnz=24 changed_nnz=3")
    sys.exit(0)
if right=="A3":
    if os.environ.get("MOCK_NO_CHANGE"):
        print("REPLAY_NEWTON_ERROR replace-values requires a genuinely changed operator")
        print("REPLAY_CSR_COMPARE compatible=0")
        sys.exit(1)
    print("REPLAY_NEWTON_ERROR replace-values CSR column graph mismatch offset=4")
    print("REPLAY_CSR_COMPARE compatible=0")
    sys.exit(1)
print("ERROR unexpected matrix names")
sys.exit(8)
PY
chmod +x "$tmp/mock.py"
run() {
    local csv="$1"
    shift
    python3 "$here/mfamg_v8_newton_preflight.py" \
      --exe "$tmp/mock.py" --modes "$tmp/modes.mtx" \
      --snapshots "$csv" "$@"
}
LIS_UNEXPECTED=1 run "$tmp/steps.csv" > "$tmp/ok.log"
grep -qx 'V8_NEWTON_PREFLIGHT=USER_DECLARED_SEQUENCE_READY' "$tmp/ok.log"
grep -qx 'V8_NEWTON_PREFLIGHT_GRAPHS=same_csr,changed_graph' "$tmp/ok.log"
grep -qx 'V8_NEWTON_CHRONOLOGY=CSV_ORDER_DECLARED_NOT_INDEPENDENTLY_PROVEN' "$tmp/ok.log"
grep -qx 'V8_NEWTON_DATA_PRIVACY=LOCAL_READ_ONLY' "$tmp/ok.log"
run "$tmp/steps.csv" --max-windows 0 > "$tmp/all.log"
grep -qx 'V8_NEWTON_PREFLIGHT_WINDOWS=1' "$tmp/all.log"
test "$(find "$tmp" -maxdepth 1 -type d | wc -l)" -eq 1
if run "$tmp/onlytwo.csv" > "$tmp/two.log" 2>&1; then exit 1; fi
grep -qx 'V8_NEWTON_PREFLIGHT=BLOCKED' "$tmp/two.log"
if run "$tmp/duplicate.csv" > "$tmp/dup.log" 2>&1; then exit 1; fi
grep -qx 'V8_NEWTON_PREFLIGHT=BLOCKED' "$tmp/dup.log"
if MOCK_NO_CHANGE=1 run "$tmp/steps.csv" > "$tmp/nochange.log" 2>&1; then exit 1; fi
grep -qx 'V8_NEWTON_PREFLIGHT=BLOCKED' "$tmp/nochange.log"
if MOCK_MALFORMED=1 run "$tmp/steps.csv" > "$tmp/malformed.log" 2>&1; then exit 1; fi
grep -qx 'V8_NEWTON_PREFLIGHT=BLOCKED' "$tmp/malformed.log"
if MOCK_MUTATE=1 run "$tmp/steps.csv" > "$tmp/mutation.log" 2>&1; then exit 1; fi
grep -qx 'V8_NEWTON_PREFLIGHT=BLOCKED' "$tmp/mutation.log"
echo MFAMG_V8_NEWTON_PREFLIGHT_SMOKE_OK
