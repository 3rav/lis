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
MODES = ("full", "dgks", "one_c", "full_repeat")
METRIC = re.compile(
    r"GCRODR_ADAPTIVE_MGS_COUNTS mode=(\d+) c_reorth=(\d+) c_skipped=(\d+) "
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
    stats = {mode:[0,0,0,0,0] for mode in ("dgks","one_c")}
    case_speedups = {"dgks":[], "one_c":[]}
    case_itdelta = {"dgks":[], "one_c":[]}
    print("mode precon side restart k stride full_ms candidate_ms speedup "
          "baseline_repeat_drift iter_delta")
    for precon in PRECONS:
        for side in SIDES:
            for restart in RESTARTS:
                data = {mode:read_cases(root, precon, side, restart, mode)
                        for mode in MODES}
                for mode in ("dgks","one_c"):
                    for k in (4, 8):
                        for stride in (1, 4):
                            full, fi = median_sequence(data["full"], "gcrodr",
                                                       k, stride, 0)
                            candidate, ai = median_sequence(data[mode], "gcrodr",
                                                            k, stride, 0)
                            rep, _ = median_sequence(data["full_repeat"],
                                                     "gcrodr", k, stride, 0)
                            speedup=full/candidate
                            drift=rep/full
                            case_speedups[mode].append(speedup)
                            case_itdelta[mode].append(ai-fi)
                            print(f"{mode} {precon} {side} {restart} {k} {stride} "
                                  f"{full*1000:.4f} {candidate*1000:.4f} "
                                  f"{speedup:.4f} {drift:.4f} {ai-fi:+.0f}")
                    log = root / mode / f"{precon}-n{side}-r{restart}.stderr"
                    for line in log.read_text(encoding="utf-8").splitlines():
                        match=METRIC.search(line)
                        if match:
                            ns=list(map(int,match.groups()))
                            if ns[0] != (1 if mode=="dgks" else 2):
                                raise ValueError(f"{log}: wrong instrument mode {ns}")
                            for i in range(5):
                                stats[mode][i]+=ns[i+1]
                            if ns[6] != 0:
                                raise ValueError(f"{log}: unsuccessful solve {ns}")
    for mode in ("dgks","one_c"):
        st=stats[mode]
        if st[1]+st[0]==0:
            raise ValueError(f"{mode}: adaptive GCRO-DR did not perform C projections")
        if st[1]==0:
            raise ValueError(f"{mode}: did not skip any C second pass")
        if mode=="dgks" and st[3]+st[2]==0:
            raise ValueError("DGKS mode never executed V projections")
        if mode=="one_c" and (st[0]!=0 or st[2]!=0 or st[3]!=0):
            raise ValueError(f"{mode}: unexpected extra-pass counter {st}")
        print(f"\n{mode}: C second executed={st[0]}, skipped={st[1]}, "
              f"V second executed={st[2]}, skipped={st[3]}")
        print(f"{mode}: median GCRO-DR speedup (32 combinations): "
              f"{statistics.median(case_speedups[mode]):.3f}x")
        print(f"{mode}: faster cases "
              f"{sum(x>1 for x in case_speedups[mode])}/"
              f"{len(case_speedups[mode])}")
        print(f"{mode}: median change in iterations: "
              f"{statistics.median(case_itdelta[mode]):+.1f}")
    print("NOTE: synthetic shared-runner CPU; no production claim or time gate.")
    return 0


if __name__=="__main__":
    try:
        sys.exit(main())
    except (OSError,ValueError,ZeroDivisionError,KeyError) as exc:
        print(f"FAIL: {exc}",file=sys.stderr)
        sys.exit(1)
