#!/usr/bin/env python3
"""Audit a paired evolving-operator GCRO-DR / FGCRO-DR / FGMRES benchmark.

All CSV rows must converge with an independently computed true residual.
Measurements are synthetic and shared-runner CPU; speedups never gate CI.
"""
import csv
import math
import statistics
import sys
from pathlib import Path

METHODS = [
    ("fgmres", 0, 1, 0),
    ("gcrodr", 4, 1, 0),
    ("gcrodr", 4, 4, 0),
    ("gcrodr", 8, 1, 0),
    ("gcrodr", 8, 4, 0),
    ("fgcrodr", 0, 1, 0),
    ("fgcrodr", 4, 1, 0),
    ("fgcrodr", 4, 4, 0),
    ("fgcrodr", 8, 1, 0),
    ("fgcrodr", 8, 4, 0),
    ("fgcrodr", 4, 4, 1),
]
FIELDS = {"precon", "side", "restart", "trial", "step", "nstep",
          "method", "k", "stride", "intra", "wall_seconds",
          "cpu_seconds", "iterations", "api_error", "solver_status",
          "true_relres", "recycle_builds", "recycle_updates",
          "recycle_cross_uses", "precon_apply_wall"}
NSTEP = 6


def audited_csv(path: Path, side: int, restart: int, precon: str,
                trials: int):
    got = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or ()) != FIELDS:
            raise ValueError(f"{path}: bad header: {reader.fieldnames}")
        for row in reader:
            method = row["method"]
            k, stride, intra = (int(row[key]) for key in ("k", "stride", "intra"))
            trial, step = int(row["trial"]), int(row["step"])
            if (row["precon"] != precon or int(row["side"]) != side
                or int(row["restart"]) != restart
                or int(row["nstep"]) != NSTEP
                or (method, k, stride, intra) not in METHODS
                or trial not in range(trials) or step not in range(NSTEP)):
                raise ValueError(f"{path}: unexpected case: {row}")
            key = (method, k, stride, intra, trial, step)
            if key in got:
                raise ValueError(f"{path}: duplicate {key}")
            cpu, wall, resid, ptime = (
                float(row[x]) for x in
                ("cpu_seconds", "wall_seconds", "true_relres",
                 "precon_apply_wall")
            )
            iterations = int(row["iterations"])
            err, status = int(row["api_error"]), int(row["solver_status"])
            updates = int(row["recycle_updates"])
            cross = int(row["recycle_cross_uses"])
            if (not all(map(math.isfinite, (cpu, wall, resid, ptime)))
                or min(cpu, wall, resid, ptime) < 0
                or err != 0 or status != 0 or not (0 < iterations <= 500)
                or resid > 1e-8 or updates < 0 or cross < 0):
                raise ValueError(f"{path}: numerical or timing failure: {row}")
            got[key] = (cpu, wall, iterations, resid, updates, cross)
    expected = len(METHODS) * trials * NSTEP
    if len(got) != expected:
        raise ValueError(f"{path}: expected {expected}, got {len(got)}")
    return got


def round_values(rows, method, k, stride, intra, trial):
    return [rows[(method, k, stride, intra, trial, step)]
            for step in range(NSTEP)]


def median_time(rows, method, k, stride, intra, trials):
    return statistics.median(sum(v[0] for v in
        round_values(rows, method, k, stride, intra, trial))
        for trial in range(1, trials))


def median_iters(rows, method, k, stride, intra, trials):
    return statistics.median(sum(v[2] for v in
        round_values(rows, method, k, stride, intra, trial))
        for trial in range(1, trials))


def median_updates(rows, method, k, stride, intra, trials):
    return statistics.median(
        round_values(rows, method, k, stride, intra, trial)[-1][4]
        for trial in range(1, trials)
    )


def main():
    if len(sys.argv) != 3:
        print("usage: audit_dynamic_recycle.py RESULTS_DIR TRIALS",
              file=sys.stderr)
        return 2
    root, trials = Path(sys.argv[1]), int(sys.argv[2])
    if trials < 2 or trials > 8:
        raise ValueError("trials out of range")
    count, maxresid = 0, 0.0
    for precon in ("jacobi", "ilutp"):
        for side in (28, 40):
            for restart in (12, 24):
                path = root / f"{precon}-n{side}-r{restart}.csv"
                data = audited_csv(path, side, restart, precon, trials)
                count += len(data)
                maxresid = max(maxresid, *(v[3] for v in data.values()))
                base = median_time(data, "fgmres", 0, 1, 0, trials)
                print(f"\nprecon={precon}, n={side*side}, restart={restart}, "
                      f"FGMRES sequence CPU={base*1000:.3f} ms")
                print("solver       k stride intra | seq_cpu_ms x_fgmres "
                      "seq_iters updates")
                for name, k, stride, intra in METHODS:
                    m = median_time(data, name, k, stride, intra, trials)
                    it = median_iters(data, name, k, stride, intra, trials)
                    upd = median_updates(data, name, k, stride, intra, trials)
                    print(f"{name:<12} {k:>2} {stride:>6} {intra:>5} | "
                          f"{m*1e3:>10.3f} {m/base:>8.3f} "
                          f"{it:>9g} {upd:>7g}")
                for name in ("gcrodr", "fgcrodr"):
                    for k in (4, 8):
                        t1 = median_time(data, name, k, 1, 0, trials)
                        t4 = median_time(data, name, k, 4, 0, trials)
                        u1 = median_updates(data, name, k, 1, 0, trials)
                        u4 = median_updates(data, name, k, 4, 0, trials)
                        print(f"cadence_effect {name} k={k}: "
                              f"seq_cpu_speedup={t1/t4:.3f}x "
                              f"ritz_updates={u1:g}->{u4:g}")
    print(f"\nPASS: {count} converged A-changing solves, "
          f"max true relative residual={maxresid:.3g}")
    print("NOTE: synthetic contact-like value updates, shared Actions "
          "runner CPU, not actual FEA/contact snapshots.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError, ZeroDivisionError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
