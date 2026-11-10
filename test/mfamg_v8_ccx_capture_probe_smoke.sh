#!/usr/bin/env bash
# Synthetic fixture only; no CalculiX build or private file access.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$tmp/ccx/src"
cat > "$tmp/ccx/src/nonlingeo.c" <<'C'
void nonlingeo(void) {
    int iinc = 1, iit = 2;
    mafillsmmain();
    preiter();
    spooles();
    matrixstorage();
}
C
cat > "$tmp/synthetic.inp" <<'INP'
** synthetic public-only example
*HEADING
*STEP
*STATIC
*CONTACT PAIR
*INCLUDE, INPUT=public_synthetic_mesh.inp
*END STEP
INP
python3 "$here/mfamg_v8_ccx_capture_probe.py" \
  --source-root "$tmp/ccx" --input-file "$tmp/synthetic.inp" > "$tmp/good.out"
grep -qx 'CCX_CAPTURE_PROBE=SOURCE_AND_MODEL_READABLE' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_SOURCE_ASSEMBLY=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_SOURCE_DIRECT_SOLVER=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_SOURCE_ITERATIVE_SOLVER=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_INP_STEP=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_INP_STATIC=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_INP_CONTACT=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_INP_INCLUDE=1' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_EXPORT_ENABLED=0' "$tmp/good.out"
grep -qx 'CCX_CAPTURE_PROVENANCE=NOT_YET_ESTABLISHED' "$tmp/good.out"
if grep -F "$tmp" "$tmp/good.out"; then
  echo "ERROR: exposed local path" >&2; exit 1
fi
if python3 "$here/mfamg_v8_ccx_capture_probe.py" \
   --source-root "$tmp/missing" --input-file "$tmp/synthetic.inp" \
   > "$tmp/missing.out" 2> "$tmp/missing.err"; then
   echo "ERROR: missing source accepted" >&2; exit 1
fi
grep -qx 'CCX_CAPTURE_PROBE=BLOCKED' "$tmp/missing.err"
if grep -F "$tmp" "$tmp/missing.err"; then
  echo "ERROR: failed probe exposed path" >&2; exit 1
fi
echo CCX_CAPTURE_PROBE_SMOKE_OK
