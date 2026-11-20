#!/usr/bin/env python3
"""Compare matched full-double-MGS and DGKS-style adaptive GCRODR trials.

Mathematical correctness is validated for every individual solve by
audit_dynamic_recycle.py first; this script reports performance and actual
reorthogonalization decision counts without gating timer differences.
"""
import csv
import math
import re
import statistics
import sys
from pathlib import Path

PRECONS = ("jacobi", "ilutp")
SIDES = (28, 40)
RESTARTS = (12, 24)
MODES = ("full", "dgks", "full_repeat")
METRIC = re.compile(
    r"GCRODR_ADAPTIVE_MGS_COUNTS c_reorth=(\d+) c_skipped=(\d+) "
    r"v_reorth=(\d+) v_skipped=(\d+) iter=(\d+) ret=(-?\d+)"
)


def read_cases(root, precon, side, restart, mode):
    path = root / mode / f"{precon}-n{side}-r{restart}.csv"
    values = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (row["method"], int(row["k"]), int(row["stride"]),
                   int(row["intra"]), int(row["trial"]), int(row["step"]))
            if key in values:
                raise ValueError(f"duplicate {key} in {path}")
            values[key] = (float(row["cpu_seconds"]), int(row["iterations"]),
                           float(row["true_relres"]), int(row["api_error"]),
                           int(row["solver_status"]))
    if len(values) != 11*4*6:
        raise ValueError(f"{path}: found {len(values)} rows; expected 264")
    if any(not math.isfinite(v[0]) or v[0] < 0 or
           v[2] > 1e-8 or v[3] or v[4] for v in values.values()):
        raise ValueError(f"{path}: bad residual, API or timing")
    return values


def median_sequence(rows, method, k, stride, intra):
    cpu, iters = [], []
    for trial in (1, 2, 3):
        recs = [rows[(method, k, stride, intra, trial, step)]
                for step in range(6)]
        cpu.append(sum(v[0] for v in recs))
        iters.append(sum(v[1] for v in recs))
    return statistics.median(cpu), statistics.median(iters)


def main():
    if len(sys.argv) != 2:
        print("usage: audit_gcrodr_adaptive_mgs.py RESULTS", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    stats = [0, 0, 0, 0, 0]
    case_speedups = []
    case_itdelta = []
    print("precon side restart k stride full_ms dgks_ms speedup "
          "baseline_repeat_drift iter_delta")
    for precon in PRECONS:
        for side in SIDES:
            for restart in RESTARTS:
                data = {mode:read_cases(root, precon, side, restart, mode)
                        for mode in MODES}
                for k in (4, 8):
                    for stride in (1, 4):
                        full, fi = median_sequence(data["full"], "gcrodr",
                                                   k, stride, 0)
                        adaptive, ai = median_sequence(data["dgks"], "gcrodr",
                                                       k, stride, 0)
                        rep, _ = median_sequence(data["full_repeat"],
                                                 "gcrodr", k, stride, 0)
                        speedup=full/adaptive
                        drift=rep/full
                        case_speedups.append(speedup)
                        case_itdelta.append(ai-fi)
                        print(f"{precon} {side} {restart} {k} {stride} "
                              f"{full*1000:.4f} {adaptive*1000:.4f} "
                              f"{speedup:.4f} {drift:.4f} {ai-fi:+.0f}")
                log = root / "dgks" / f"{precon}-n{side}-r{restart}.stderr"
                for line in log.read_text(encoding="utf-8").splitlines():
                    match=METRIC.search(line)
                    if match:
                        ns=list(map(int,match.groups()))
                        for i in range(5):
                            stats[i]+=ns[i]
                        if ns[5] != 0:
                            raise ValueError(f"{log}: non-success MGS solve {ns}")
    if stats[1]+stats[0]==0 or stats[3]+stats[2]==0:
        raise ValueError("adaptive GCRO-DR never executed C/V basis projections")
    if stats[1]==0 or stats[3]==0:
        raise ValueError("adaptive MGS did not skip any second pass")
    print(f"\nAdaptive C: second-pass executed={stats[0]}, skipped={stats[1]}")
    print(f"Adaptive V: second-pass executed={stats[2]}, skipped={stats[3]}")
    print(f"Median GCRO-DR speedup (32 combinations): "
          f"{statistics.median(case_speedups):.3f}x")
    print(f"Faster cases: {sum(x>1 for x in case_speedups)}/{len(case_speedups)}")
    print(f"Median change in iterations: "
          f"{statistics.median(case_itdelta):+.1f}")
    print("NOTE: synthetic shared-runner CPU; no production claim or time gate.")
    return 0


if __name__=="__main__":
    try:
        sys.exit(main())
    except (OSError,ValueError,ZeroDivisionError,KeyError) as exc:
        print(f"FAIL: {exc}",file=sys.stderr)
        sys.exit(1)
