#!/usr/bin/env python3
"""Validate and summarize a diagnostic Ritz-update-cadence comparison.

Benchmarks are synthetic Jacobi-preconditioned CPU timings on shared CI
runners. Never interpret these medians as production solver speedups.
"""
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

METHODS = {
    ("fgmres", 0, 0),
    ("gcrodr", 0, 0),
    ("gcrodr", 4, 0),
    ("gcrodr", 8, 0),
    ("fgcrodr", 0, 0),
    ("fgcrodr", 4, 0),
    ("fgcrodr", 8, 0),
    ("fgcrodr", 4, 1),
}
RESTARTS = (12, 24)
MODES = ("1", "4", "once", "1_repeat")
SIZES = (28, 40)
SOLVES = range(4)


def load(folder: Path):
    records = {}
    for side in SIZES:
        for mode in MODES:
            path = folder / f"n{side}-{mode}.csv"
            if not path.is_file():
                raise ValueError(f"missing {path}")
            with path.open(newline="", encoding="utf-8") as file:
                rows = list(csv.DictReader(file))
            if len(rows) != len(METHODS)*len(RESTARTS)*len(SOLVES):
                raise ValueError(f"{path}: expected 64 rows, got {len(rows)}")
            for row in rows:
                method = row["method"]
                k = int(row["k"])
                intra = int(row["intra"])
                restart = int(row["restart"])
                solve = int(row["solve_index"])
                size = int(row["side"])
                cpu = float(row["cpu_seconds"])
                wall = float(row["wall_seconds"])
                residual = float(row["true_relres"])
                iters = int(row["iterations"])
                api = int(row["api_error"])
                status = int(row["solver_status"])
                key = (size, mode, method, k, intra, restart, solve)
                if size != side or (method, k, intra) not in METHODS:
                    raise ValueError(f"{path}: unexpected method or size: {row}")
                if restart not in RESTARTS or solve not in SOLVES or key in records:
                    raise ValueError(f"{path}: missing/duplicated key {key}")
                if (not all(map(math.isfinite, (cpu, wall, residual)))
                    or cpu < 0 or wall < 0 or iters <= 0
                    or api != 0 or status != 0 or residual > 1e-8):
                    raise ValueError(f"{path}: failure or invalid values: {row}")
                records[key] = (cpu, wall, iters, residual)
    return records


def med(rec, side, mode, method, k, intra, restart, field):
    return statistics.median(rec[(side, mode, method, k, intra, restart, solve)][field]
                             for solve in (1, 2, 3))


def main():
    if len(sys.argv) != 2:
        print("usage: audit_recycle_cadence.py ARTIFACT_DIR", file=sys.stderr)
        return 2
    records = load(Path(sys.argv[1]))
    print(f"PASS: {len(records)} converged numerical records; "
          f"all true residuals <= 1e-8; complete baseline/4/once/repeat coverage")
    for side in SIZES:
        for restart in RESTARTS:
            baseline = med(records, side, "1", "fgmres", 0, 0, restart, 0)
            print(f"\nside={side}, n={side*side}, restart={restart}, "
                  f"FGMRES median={baseline*1e3:.3f} ms")
            print("method              k intra | base_ms  every4_ms  once_ms | "
                  "every4/base once/base | base_iter every4_iter once_iter")
            for method, k, intra in sorted(METHODS):
                times = [med(records, side, mode, method, k, intra, restart, 0)*1e3
                         for mode in ("1", "4", "once")]
                iterations = [med(records, side, mode, method, k, intra, restart, 2)
                              for mode in ("1", "4", "once")]
                print(f"{method:<18} {k:>2} {intra:>5} | "
                      f"{times[0]:7.3f} {times[1]:10.3f} {times[2]:8.3f} | "
                      f"{times[1]/times[0]:11.3f} {times[2]/times[0]:9.3f} | "
                      f"{iterations[0]:9g} {iterations[1]:11g} {iterations[2]:9g}")
            drift = med(records, side, "1_repeat", "fgmres", 0, 0, restart, 0)/baseline
            print(f"FGMRES baseline-repeat drift factor: {drift:.3f}")
    # Do not gate experimental performance or iteration changes; only
    # mathematical convergence and complete comparability.
    print("\nNOTE: synthetic shared-runner CPU times; no production claims.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (KeyError, ValueError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
