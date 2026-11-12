#!/usr/bin/env python3
"""Local paired FGMRES-control vs opt-in FGCRO-DR replay on one fixed A/b.

Default: SAME v8 executable for both methods (isolates solver algorithm).
Optional --exe-v7: historical v7 FGMRES executable, distinct binary.
Neither comparison alone proves whole-CalculiX/Newton speedup.
No matrices, modes, or raw logs are uploaded.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

from mfamg_mode_audit import read_modes
from mfamg_replay_v8_repeats import (
    ENV, assert_inputs_unchanged, fingerprint_inputs, integer, kv, parse, real,
)

VARIANTS = ("fgmres_control", "fgcrodr_v8")
FIELDS = ("pair", "sequence", "variant", "n", "nnz", "modes",
          "process_exit", "solve_rc", "status", "iter", "true_relative",
          "wall_seconds", "converged", "builds", "updates", "uses",
          "rejects", "orth_error", "relation_error", "log_file")


def baseline_parse(output, exit_code, tol, logfile, expected_modes, historical_v7, maxiter=3000):
    """A v7 replay does NOT print method/recycle fields; v8 FGMRES does."""
    groups = {tag: [] for tag in ("REPLAY", "REPLAY_RESULT",
                                   "REPLAY_RECYCLE_STATS", "REPLAY_RECYCLE_CHECK")}
    for line in output.splitlines():
        for tag in groups:
            if line.startswith(tag + " "):
                groups[tag].append(kv(line))
    if (len(groups["REPLAY"]) != 1 or len(groups["REPLAY_RESULT"]) != 1 or
            groups["REPLAY_RECYCLE_STATS"] or groups["REPLAY_RECYCLE_CHECK"]):
        raise ValueError("FGMRES control requires one header/result and no recycle diagnostic")
    h, rr = groups["REPLAY"][0], groups["REPLAY_RESULT"][0]
    if historical_v7:
        if "method" in h or "recycle_k" in h:
            raise ValueError("historical v7 baseline has unexpected v8 header fields")
    elif h.get("method") != "fgmres" or integer(h, "recycle_k") != 0:
        raise ValueError("same-executable FGMRES control did not select fgmres/k0")
    n, nnz, modes = (integer(h, key) for key in ("n", "nnz", "modes"))
    declared_tol = real(h, "tol")
    if (n <= 0 or nnz <= 0 or modes != expected_modes or
            integer(h, "restart") != 180 or integer(h, "maxiter") != maxiter or
            not math.isclose(declared_tol, tol, rel_tol=1e-12, abs_tol=0.0)):
        raise ValueError("unexpected FGMRES control configuration")
    rc, status, iters = (integer(rr, key) for key in ("solve_rc", "status", "iter"))
    relative, wall = real(rr, "true_relative"), real(rr, "wall_seconds")
    if (exit_code not in (0, 1) or
            (exit_code == 0 and (rc != 0 or status != 0)) or
            iters < 0 or relative < 0 or wall < 0):
        raise ValueError("inconsistent FGMRES control result/exit status")
    converged = int(exit_code == 0 and rc == 0 and status == 0 and
                    iters > 0 and relative <= tol and wall > 0)
    return dict(n=n, nnz=nnz, modes=modes, process_exit=exit_code,
                solve_rc=rc, status=status, iter=iters, true_relative=relative,
                wall_seconds=wall, converged=converged,
                builds=None, updates=None, uses=None, rejects=None,
                orth_error=None, relation_error=None, log_file=logfile)


def assess(rows, pairs, historical_v7):
    groups = {name: [r for r in rows if r["variant"] == name] for name in VARIANTS}
    complete = (len(rows) == pairs * 2 and all(
        len(groups[name]) == pairs and
        all(sum(r["pair"] == idx for r in groups[name]) == 1
            for idx in range(1, pairs + 1))
        for name in VARIANTS))
    accepted = complete and all(r["converged"] for r in rows)
    metrics = {}
    for name, group in groups.items():
        metrics[name] = dict(
            runs=len(group),
            converged=sum(r["converged"] for r in group),
            median_iterations=statistics.median(r["iter"] for r in group)
            if accepted else None,
            median_wall_seconds=statistics.median(r["wall_seconds"] for r in group)
            if accepted else None,
            median_true_relative=statistics.median(r["true_relative"] for r in group)
            if group else None,
        )
    paired = None
    median_difference = None
    if accepted:
        paired = []
        for idx in range(1, pairs + 1):
            baseline = next(r["wall_seconds"] for r in groups["fgmres_control"]
                            if r["pair"] == idx)
            v8 = next(r["wall_seconds"] for r in groups["fgcrodr_v8"]
                      if r["pair"] == idx)
            paired.append(100.0 * (1.0 - v8 / baseline))
        median_difference = 100.0 * (
            1.0 - metrics["fgcrodr_v8"]["median_wall_seconds"] /
            metrics["fgmres_control"]["median_wall_seconds"])
    return dict(
        assessment="MATCHED_CONVERGED" if accepted else "INCONCLUSIVE_NONCONVERGENCE",
        baseline_type="historical_v7_binary" if historical_v7 else "same_v8_binary_fgmres",
        metrics=metrics,
        paired_v8_wall_reduction_percent_by_pair=paired,
        median_paired_v8_wall_reduction_percent=statistics.median(paired)
        if paired is not None else None,
        v8_wall_reduction_percent_from_separate_medians=median_difference,
        caution=("Fixed-A replay only; three or more exploratory pairs do not establish "
                 "statistical significance, changing-operator safety, or CalculiX speedup. "
                 "Same-v8-binary FGMRES control is NOT a historical v7 comparison.")
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("exe-v8", "matrix", "rhs", "modes", "out-dir"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--exe-v7", type=Path, help="Optional true v7 FGMRES binary; default uses v8 FGMRES")
    p.add_argument("--expected-modes", type=int, default=18)
    p.add_argument("--pairs", type=int, default=3)
    p.add_argument("--maxiter", type=int, default=3000, help="Common FGMRES/FGCRODR MAXITER (100..6000)")
    p.add_argument("--fine-agg-factor", type=int, choices=(2, 4), default=4,
                   help="Identical MFAMG aggregation factor (2 for historical 116-spring, 4 for 31-spring)")
    p.add_argument("--tol", type=float, default=1e-5)
    p.add_argument("--require-both-converged", action="store_true")
    a = p.parse_args()
    if not (1 <= a.pairs <= 30 and 100 <= a.maxiter <= 6000 and 1 <= a.expected_modes <= 256 and
            math.isfinite(a.tol) and 0 < a.tol < 1):
        p.error("pairs 1..30, maxiter 100..6000, expected modes 1..256, tol finite in (0,1)")
    try:
        paths = dict(exe_v8=a.exe_v8, matrix=a.matrix, rhs=a.rhs, modes=a.modes)
        if a.exe_v7 is not None:
            paths["exe_v7"] = a.exe_v7
        for label, path in paths.items():
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError("missing/empty input " + label + ": " + str(path))
        if len({path.resolve() for path in paths.values()}) != len(paths):
            raise ValueError("all paths must refer to distinct files")
        n, modes = read_modes(a.modes)
        if n <= 0 or len(modes) != a.expected_modes:
            raise ValueError("unexpected mode dimensions")
        hashes = fingerprint_inputs(paths)
        if a.out_dir.exists():
            raise ValueError("output directory already exists: " + str(a.out_dir))
        a.out_dir.mkdir(parents=True)
        # Keep POSIX permission restriction where available. Windows uses ACLs.
        if os.name != "nt":
            os.chmod(a.out_dir, 0o700)
        fixed = {k: v for k, v in ENV.items()
                 if k not in ("LIS_MFAMG_REPLAY_METHOD",
                              "LIS_MFAMG_REPLAY_RECYCLE_DIM",
                              "LIS_FGCRODR_INTRA_SOLVE_RECYCLE",
                              "LIS_FGCRODR_INTRA_TRACE",
                              "LIS_MFAMG_REPLAY_VALIDATE_RECYCLE")}
        fixed["LIS_MFAMG_REPLAY_MAXITER"] = str(a.maxiter)
        fixed["LIS_MFAMG_FINE_AGG_FACTOR"] = str(a.fine_agg_factor)
        cfg = {k: v for k, v in os.environ.items() if not k.startswith("LIS_")}
        manifest = dict(
            schema="mfamg-fgmres-vs-fgcrodr-paired-v1",
            inputs_sha256_and_size=hashes,
            baseline_type="historical_v7_binary" if a.exe_v7 else "same_v8_binary_fgmres",
            pairs=a.pairs, expected_modes=a.expected_modes, maxiter=a.maxiter,
            fine_agg_factor=a.fine_agg_factor, mode_rows=n,
            tolerance=a.tol, common_LIS_environment=fixed,
            experiment="single_fixed_A_b; independent fresh processes; no Newton",
            note="Private matrix/modes/logs remain local; do not upload",
        )
        (a.out_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        rows = []
        with (a.out_dir / "paired_runs.csv").open("x", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=FIELDS)
            writer.writeheader()
            output.flush()
            for pair in range(1, a.pairs + 1):
                order = VARIANTS if pair % 2 else tuple(reversed(VARIANTS))
                for seq, variant in enumerate(order, 1):
                    env = dict(cfg)
                    env.update(fixed)
                    if variant == "fgcrodr_v8":
                        env.update({k: ENV[k] for k in (
                            "LIS_MFAMG_REPLAY_METHOD", "LIS_MFAMG_REPLAY_RECYCLE_DIM",
                            "LIS_FGCRODR_INTRA_SOLVE_RECYCLE", "LIS_FGCRODR_INTRA_TRACE",
                            "LIS_MFAMG_REPLAY_VALIDATE_RECYCLE")})
                        exe = a.exe_v8
                    else:
                        exe = a.exe_v7 or a.exe_v8
                        if a.exe_v7 is None:
                            env["LIS_MFAMG_REPLAY_METHOD"] = "fgmres"
                    logfile = "pair%02d_run%d_%s.log" % (pair, seq, variant)
                    result = subprocess.run(
                        [str(exe.resolve()), str(a.matrix.resolve()),
                         str(a.rhs.resolve()), str(a.modes.resolve())],
                        cwd=str(a.matrix.resolve().parent), env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        text=True, errors="replace", check=False)
                    (a.out_dir / logfile).write_text(result.stdout, encoding="utf-8")
                    assert_inputs_unchanged(paths, hashes)
                    if variant == "fgcrodr_v8":
                        row = parse(result.stdout, result.returncode, a.tol, logfile,
                                    a.expected_modes, maxiter=a.maxiter)
                    else:
                        row = baseline_parse(result.stdout, result.returncode, a.tol,
                                             logfile, a.expected_modes, a.exe_v7 is not None,
                                             maxiter=a.maxiter)
                    if row["n"] != n or (rows and (row["n"], row["nnz"]) !=
                                         (rows[0]["n"], rows[0]["nnz"])):
                        raise ValueError("inconsistent imported matrix dimensions/nnz")
                    row.update(pair=pair, sequence=seq, variant=variant)
                    rows.append(row)
                    writer.writerow(row)
                    output.flush()
                    print("V8_SOLVER_PAIR pair=%d case=%s status=%d iter=%d wall=%.3f converged=%d" %
                          (pair, variant, row["status"], row["iter"],
                           row["wall_seconds"], row["converged"]), flush=True)
        assert_inputs_unchanged(paths, hashes)
        report = assess(rows, a.pairs, a.exe_v7 is not None)
        report.update(pairs=a.pairs, manifest_file="manifest.json")
        (a.out_dir / "assessment.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("V8_SOLVER_PAIR_ASSESSMENT=" + report["assessment"])
        print("V8_SOLVER_PAIR_SUMMARY=" + str(a.out_dir / "assessment.json"))
        return int(a.require_both_converged and report["assessment"] != "MATCHED_CONVERGED")
    except (OSError, ValueError, OverflowError) as exc:
        print("V8_SOLVER_PAIR_ERROR=" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
