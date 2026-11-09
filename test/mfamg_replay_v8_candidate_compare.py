#!/usr/bin/env python3
"""Local 12/14/18 MFAMG candidate replay on a FIXED operator, not Newton steps."""
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
from mfamg_mode_ritz_selection import orthonormalize
from mfamg_replay_v8_repeats import ENV, parse, fingerprint_inputs, assert_inputs_unchanged

VARIANTS = ("baseline12", "generated14", "offline18")
FIELDS = ("repetition", "sequence", "variant", "n", "nnz", "modes",
          "process_exit", "solve_rc", "status", "iter", "true_relative",
          "wall_seconds", "builds", "updates", "uses", "rejects",
          "orth_error", "relation_error", "converged", "log_file")


def preflight(baseline, generated, offline):
    n0, base = read_modes(baseline)
    ng, auto = read_modes(generated)
    no, control = read_modes(offline)
    if n0 != ng or n0 != no or len(base) != 12 or len(auto) != 14 or len(control) != 18:
        raise ValueError("expected common row count and 12/14/18 mode columns")
    if auto[:12] != base or control[:12] != base:
        raise ValueError("generated and offline first 12 columns must equal baseline exactly")
    original_q = orthonormalize(base, tolerance=1e-10)
    orthonormalize(auto[12:], prior=original_q, tolerance=1e-10)
    orthonormalize(control[12:], prior=original_q, tolerance=1e-10)
    return n0


def assess(rows, repeats):
    groups = {name: [r for r in rows if r["variant"] == name] for name in VARIANTS}
    complete = len(rows) == 3 * repeats and all(len(group) == repeats for group in groups.values())
    metrics = {}
    for name, group in groups.items():
        all_ok = complete and all(r["converged"] for r in group)
        metrics[name] = {
            "runs": len(group),
            "converged": sum(r["converged"] for r in group),
            "all_converged": all_ok,
            "median_true_relative": statistics.median(
                r["true_relative"] for r in group) if group else None,
            "median_iterations_if_all_converged": statistics.median(
                r["iter"] for r in group) if all_ok else None,
            "median_wall_if_all_converged": statistics.median(
                r["wall_seconds"] for r in group) if all_ok else None,
        }
    comparison = None
    g, o = metrics["generated14"], metrics["offline18"]
    if g["all_converged"] and o["all_converged"]:
        comparison = 100 * (1 - g["median_wall_if_all_converged"] /
                            o["median_wall_if_all_converged"])
    status = ("CONTROL_CONVERGED" if metrics["offline18"]["all_converged"]
              else "INCONCLUSIVE_OFFLINE_CONTROL")
    return {
        "assessment": status if complete else "INCONCLUSIVE_INCOMPLETE",
        "metrics": metrics,
        "generated_vs_offline_percent_median_wall_reduction": comparison,
        "caution": "Descriptive fixed-A replay only; unequal modes, no proof of speedup or Newton safety.",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("exe", "matrix", "rhs", "baseline", "generated", "offline-control", "out-dir"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--repetitions", type=int, default=1)
    p.add_argument("--tol", type=float, default=1e-5)
    p.add_argument("--require-offline-control-converged", action="store_true")
    args = p.parse_args()
    if not (1 <= args.repetitions <= 10 and math.isfinite(args.tol)
            and 0 < args.tol < 1):
        p.error("repetitions must be 1..10 and tolerance in (0,1)")
    try:
        paths = {
            "exe": args.exe, "matrix": args.matrix, "rhs": args.rhs,
            "baseline12": args.baseline, "generated14": args.generated,
            "offline18": args.offline_control,
        }
        for key, path in paths.items():
            if not path.is_file() or not path.stat().st_size:
                raise ValueError("missing/empty %s: %s" % (key, path))
        if len({path.resolve() for path in paths.values()}) != len(paths):
            raise ValueError("input paths must all be distinct")
        n = preflight(args.baseline, args.generated, args.offline_control)
        fingerprints = fingerprint_inputs(paths)
        if args.out_dir.exists():
            raise ValueError("output directory already exists: " + str(args.out_dir))
        args.out_dir.mkdir(parents=True)
        os.chmod(args.out_dir, 0o700)
        manifest = {
            "schema": "mfamg-v8-automatic-candidate-control-v1",
            "input_fingerprints": fingerprints,
            "fixed_solver_environment": ENV,
            "preflight": "baseline_preserved_and_extra_modes_independent",
            "n": n,
            "repetitions": args.repetitions,
            "note": "Local only. Fresh-process fixed-A solves, NOT a Newton sequence",
        }
        (args.out_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        env = {k: v for k, v in os.environ.items() if not k.startswith("LIS_")}
        env.update(ENV)
        modes = dict(zip(VARIANTS, (args.baseline, args.generated, args.offline_control)))
        rows = []
        with (args.out_dir / "candidate_runs.csv").open(
                "x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            stream.flush()
            for repetition in range(1, args.repetitions + 1):
                order = VARIANTS if repetition % 2 else tuple(reversed(VARIANTS))
                for sequence, variant in enumerate(order, 1):
                    path = modes[variant]
                    logfile = "rep%02d_run%d_%s.log" % (repetition, sequence, variant)
                    result = subprocess.run(
                        [str(args.exe.resolve()), str(args.matrix.resolve()),
                         str(args.rhs.resolve()), str(path.resolve())],
                        cwd=str(args.matrix.resolve().parent), env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        text=True, errors="replace", check=False)
                    (args.out_dir / logfile).write_text(result.stdout, encoding="utf-8")
                    assert_inputs_unchanged(paths, fingerprints)
                    row = parse(result.stdout, result.returncode, args.tol,
                                logfile, {"baseline12": 12, "generated14": 14,
                                          "offline18": 18}[variant])
                    if row["n"] != n:
                        raise ValueError("replay rows do not match mode rows")
                    if rows and (row["n"], row["nnz"]) != (rows[0]["n"], rows[0]["nnz"]):
                        raise ValueError("replay matrix dimensions/nnz differ between variants")
                    row.update(repetition=repetition, sequence=sequence, variant=variant)
                    rows.append(row)
                    writer.writerow(row)
                    stream.flush()
                    print("V8_CANDIDATE rep=%d variant=%s iter=%d residual=%.9g status=%d converged=%d" %
                          (repetition, variant, row["iter"], row["true_relative"],
                           row["status"], row["converged"]), flush=True)
        assert_inputs_unchanged(paths, fingerprints)
        result = assess(rows, args.repetitions)
        result["manifest_file"] = "manifest.json"
        result["repetitions"] = args.repetitions
        (args.out_dir / "assessment.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("V8_CANDIDATE_ASSESSMENT=" + result["assessment"])
        print("V8_CANDIDATE_SUMMARY=" + str(args.out_dir / "assessment.json"))
        if args.require_offline_control_converged and result["assessment"] != "CONTROL_CONVERGED":
            return 1
        return 0
    except (OSError, ValueError, OverflowError) as exc:
        print("V8_CANDIDATE_ERROR=" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
