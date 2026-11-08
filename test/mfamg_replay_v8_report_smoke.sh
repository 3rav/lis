#!/usr/bin/env bash
# Pure log-parser regression; no actual contact data and no solver required.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
work="$(mktemp -d)"
trap 'rm -rf -- "$work"' EXIT
cat >"$work/a.log" <<'LOG'
REPLAY_SWEEP_META case=31_12m_f2_fgmres snapshot=31 factor=2
REPLAY n=60265 nnz=12345 modes=12 restart=180 maxiter=1000 tol=1e-5 method=fgmres recycle_k=0
REPLAY_RESULT solve_rc=4 status=4 iter=1000 reported_residual=1.832690841870e-03 true_relative=1.832690841870e-03 wall_seconds=12.345
REPLAY_SWEEP_EXIT process_exit=1
LOG
cat >"$work/b.log" <<'LOG'
REPLAY_SWEEP_META case=116_12m_f2_fgcrodr_k8 snapshot=116 factor=2
LIS_FGCRODR_INTRA_RECYCLE iter=180 k=8 relative=9.000000e-03
LIS_FGCRODR_INTRA_RECYCLE iter=360 k=8 relative=8.000000e-04
REPLAY n=90000 nnz=24680 modes=12 restart=180 maxiter=1000 tol=1e-5 method=fgcrodr recycle_k=8
REPLAY_RECYCLE dim=8 action=2 quality=1.000e+00 min_quality=1.000e+00
REPLAY_RECYCLE_STATS builds=1 updates=2 uses=1 rejects=0
REPLAY_RECYCLE_CHECK rc=0 orth_error=4.440892098501e-16 relation_error=1.882138143748e-13
REPLAY_RESULT solve_rc=0 status=0 iter=592 reported_residual=9.000000000000e-06 true_relative=9.000000000000e-06 wall_seconds=22.500
REPLAY_SWEEP_EXIT process_exit=0
LOG
awk -f "$here/mfamg_replay_v8_report.awk" "$work/a.log" "$work/b.log" >"$work/report.csv"
test "$(wc -l < "$work/report.csv")" -eq 3
grep -Fq '31_12m_f2_fgmres,31,60265,12,2,fgmres,0,180,1000,1,4,4,1000,1.832690841870e-03,1.832690841870e-03,12.345' "$work/report.csv"
grep -Fq '116_12m_f2_fgcrodr_k8,116,90000,12,2,fgcrodr,8,180,1000,0,0,0,592,9.000000000000e-06,9.000000000000e-06,22.500,8,1,2,1,0,2,0,4.440892098501e-16,1.882138143748e-13,b.log' "$work/report.csv"
# A missing result must not silently look like a success.
echo 'REPLAY n=240 method=fgmres' >"$work/malformed.log"
if awk -f "$here/mfamg_replay_v8_report.awk" "$work/malformed.log" >"$work/bad.csv" 2>/dev/null; then
  echo "ERROR: malformed input accepted" >&2
  exit 1
fi
echo 'REPLAY_CSV_REPORT_SMOKE_OK'
