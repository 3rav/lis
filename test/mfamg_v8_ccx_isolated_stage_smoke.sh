#!/usr/bin/env bash
# Tests source staging on SYNTHETIC local fixture only.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
mkdir -p "$tmp/install/left/src" "$tmp/install/right/src" "$tmp/install/ccx" "$tmp/install/MKL_PARDISO"
echo 'not a real executable' > "$tmp/install/ccx/ccx.exe"
echo 'BUILDER=synthetic-only' > "$tmp/install/left/Makefile"
echo 'BUILDER=synthetic-only' > "$tmp/install/right/Makefile"
python3 - "$tmp" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1])
for sub in ('left','right'):
    folder=root/'install'/sub/'src'
    for name in ['nonlingeo.c','mafillsmmain.c','CalculiX.h','ccx.c']:
        (folder/name).write_text(f'/* synthetic {name} */\n',encoding='ascii')
    for i in range(52):
        (folder/f'synth_{i:03d}.c').write_text('/* synthetic only */\n',encoding='ascii')
PY
run_stage() {
  python3 "$here/mfamg_v8_ccx_isolated_stage.py" \
    --source-dir "$tmp/install/left/src" \
    --peer-dir "$tmp/install/right/src" \
    --install-dir "$tmp/install" \
    --installed-exe "$tmp/install/ccx/ccx.exe" \
    --mkl-dir "$tmp/install/MKL_PARDISO" \
    --out-dir "$tmp/stage"
}
run_stage >"$tmp/ok.stdout"
grep -qx 'CCX_STAGE=ISOLATED_SOURCE_COPY_READY' "$tmp/ok.stdout"
grep -qx 'CCX_STAGE_SOURCE_CANDIDATES_MATCH=1' "$tmp/ok.stdout"
grep -qx 'CCX_STAGE_SOURCE_FILES=56' "$tmp/ok.stdout"
grep -qx 'CCX_STAGE_ROOT_MAKEFILES=1' "$tmp/ok.stdout"
grep -qx 'CCX_STAGE_BUILD=NOT_ATTEMPTED' "$tmp/ok.stdout"
grep -qx 'CCX_STAGE_BINARY_CORRESPONDENCE=NOT_PROVEN' "$tmp/ok.stdout"
test -s "$tmp/stage/src/nonlingeo.c"
test -s "$tmp/stage/Makefile"
test ! -e "$tmp/stage/ccx.exe"
if run_stage >"$tmp/blocked.stdout" 2>"$tmp/blocked.stderr"; then
  echo 'ERROR: existing output overwritten' >&2; exit 1
fi
grep -qx 'CCX_STAGE=BLOCKED' "$tmp/blocked.stderr"
test -s "$tmp/stage/src/nonlingeo.c"
# Original sources must be identical, mismatch must reject before writes.
echo '/* modified */' >>"$tmp/install/right/src/nonlingeo.c"
if python3 "$here/mfamg_v8_ccx_isolated_stage.py" \
    --source-dir "$tmp/install/left/src" \
    --peer-dir "$tmp/install/right/src" \
    --install-dir "$tmp/install" \
    --installed-exe "$tmp/install/ccx/ccx.exe" \
    --mkl-dir "$tmp/install/MKL_PARDISO" \
    --out-dir "$tmp/other_stage" >"$tmp/mismatch.stdout" 2>"$tmp/mismatch.stderr"; then
  echo 'ERROR: different source peers accepted' >&2; exit 1
fi
test ! -e "$tmp/other_stage"
grep -qx 'CCX_STAGE=BLOCKED' "$tmp/mismatch.stderr"
if grep -F "$tmp" "$tmp/ok.stdout" "$tmp/blocked.stderr" "$tmp/mismatch.stderr"; then
  echo 'ERROR: printed private path' >&2; exit 1
fi
echo CCX_ISOLATED_STAGE_SMOKE_OK
