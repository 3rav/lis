#!/usr/bin/env bash
# Independent 8-row offline basis / mocked runner test.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
python3 - "$tmp" <<'PY'
import pathlib, sys
root = pathlib.Path(sys.argv[1])
n = 8
def write(path, cols):
    with (root / path).open('w', encoding='ascii') as f:
        f.write('%%MatrixMarket matrix array real general\n')
        f.write('%d %d\n' % (n, len(cols)))
        for j in cols:
            for i in range(n):
                f.write('%d\n' % (i == j))
write('base.mtx', [0, 1])
write('full.mtx', list(range(n)))
write('negative.mtx', list(range(n-1)))
PY
mkdir -p "$tmp/v8_mode_audit/loo" "$tmp/challenge_31"
python3 "$here/mfamg_mode_ritz_leave_one_out.py" \
  --baseline "$tmp/base.mtx" --rotated "$tmp/full.mtx" \
  --negative-control "$tmp/negative.mtx" \
  --output-dir "$tmp/v8_mode_audit/loo" > "$tmp/generated.txt"
grep -q '^RITZ_LOO_OMIT_06_MATCHES_PREVIOUS_NEGATIVE_CONTROL=YES$' "$tmp/generated.txt"
grep -q '^RITZ_LOO_REPLAY_NEW_CASES=01,02,03,04,05$' "$tmp/generated.txt"
PYTHONPATH="$here" python3 - "$tmp" <<'PY'
import pathlib, sys
from mfamg_mode_audit import read_modes
root = pathlib.Path(sys.argv[1])
_, full = read_modes(root/'full.mtx')
_, negative = read_modes(root/'negative.mtx')
for i in range(1, 7):
    f = root/'v8_mode_audit'/'loo'/('modes_omit_ritz_%02d.mtx' % i)
    n, v = read_modes(f)
    assert n == 8 and len(v) == 7
    assert v == full[:i+1] + full[i+2:]
_, last = read_modes(root/'v8_mode_audit'/'loo'/'modes_omit_ritz_06.mtx')
assert last == negative
PY
if python3 "$here/mfamg_mode_ritz_leave_one_out.py" \
  --baseline "$tmp/base.mtx" --rotated "$tmp/full.mtx" \
  --negative-control "$tmp/negative.mtx" --output-dir "$tmp/v8_mode_audit/loo" \
  > "$tmp/retry.out" 2> "$tmp/retry.err"; then
  echo "ERROR: overwrite permitted" >&2; exit 1
fi
grep -q '^RITZ_LOO_ERROR=output exists:' "$tmp/retry.err"
python3 - "$tmp" <<'PY'
import pathlib, sys
root = pathlib.Path(sys.argv[1])
lines = (root/'full.mtx').read_text().splitlines()
lines[2] = '0'
(root/'altered.mtx').write_text('\n'.join(lines)+'\n')
PY
mkdir "$tmp/invalid"
if python3 "$here/mfamg_mode_ritz_leave_one_out.py" \
  --baseline "$tmp/base.mtx" --rotated "$tmp/altered.mtx" \
  --negative-control "$tmp/negative.mtx" --output-dir "$tmp/invalid" \
  > "$tmp/bad.out" 2> "$tmp/bad.err"; then
  echo "ERROR: corrupt baseline accepted" >&2; exit 1
fi
grep -q '^RITZ_LOO_ERROR=rotated full basis changed' "$tmp/bad.err"
test "$(find "$tmp/invalid" -type f | wc -l)" -eq 0
echo '% synthetic A' > "$tmp/challenge_31/A04.mtx"
echo '% synthetic b' > "$tmp/challenge_31/b04.mtx"
cat > "$tmp/mock-replay" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
test "${LIS_MFAMG_FINE_AGG_FACTOR:-}" = "4"
test "${LIS_MFAMG_REPLAY_MAXITER:-}" = "3000"
test "${LIS_MFAMG_REPLAY_METHOD:-}" = "fgcrodr"
test "${LIS_MFAMG_REPLAY_RECYCLE_DIM:-}" = "8"
test "${LIS_FGCRODR_INTRA_SOLVE_RECYCLE:-}" = "1"
test "${LIS_MFAMG_REPLAY_VALIDATE_RECYCLE:-}" = "1"
test -s "$1" && test -s "$2" && test -s "$3"
echo "REPLAY n=8 nnz=8 modes=7 restart=180 maxiter=3000 tol=1e-5 method=fgcrodr recycle_k=8"
echo "REPLAY_RECYCLE_CHECK rc=0 orth_error=0 relation_error=0"
echo "REPLAY_RESULT solve_rc=0 status=4 iter=3001 reported_residual=1e-4 true_relative=1e-4 wall_seconds=0"
exit 1
MOCK
chmod +x "$tmp/mock-replay"
MFAMG_REPLAY_DATA="$tmp" MFAMG_REPLAY_EXE="$tmp/mock-replay" \
  bash "$here/mfamg_replay_ritz_leave_one_out.sh" > "$tmp/run.txt"
grep -q '^RITZ_LOO_REPLAY_COMPLETE$' "$tmp/run.txt"
test "$(find "$tmp/v8_mode_audit/loo/logs" -name '*.log' | wc -l)" -eq 5
if MFAMG_REPLAY_DATA="$tmp" MFAMG_REPLAY_EXE="$tmp/mock-replay" \
  bash "$here/mfamg_replay_ritz_leave_one_out.sh" \
  > "$tmp/rerun.out" 2> "$tmp/rerun.err"; then
  echo "ERROR: replay log overwrite permitted" >&2; exit 1
fi
grep -q '^RITZ_LOO_REPLAY_ERROR=log_already_exists' "$tmp/rerun.err"
echo RITZ_LOO_SMOKE_OK
