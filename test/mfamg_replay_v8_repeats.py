#!/usr/bin/env python3
"""Matched private-contact 17/18 offline Ritz replay: AB/BA order.

Outputs local logs and statistics; never uploads matrix or mode data.
This is a reproducibility tool, NOT online mode discovery.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import platform
from pathlib import Path
import statistics
import subprocess
import sys

from mfamg_mode_audit import read_modes

ENV = {
    "LIS_MFAMG_EXPERIMENT_FINE_AGG": "1",
    "LIS_MFAMG_EXPERIMENT_DEEP_MERGE": "1",
    "LIS_MFAMG_FINE_AGG_FACTOR": "4",
    "LIS_MFAMG_STAR_SCHUR": "0",
    "LIS_MFAMG_VCYCLE_SMOOTH_STEPS": "3",
    "LIS_MFAMG_VCYCLE_DIAG": "0",
    "LIS_MFAMG_VCYCLE_TRACE": "0",
    "LIS_MFAMG_DIAG": "1",
    "LIS_FGMRES_TRACE": "1",
    "LIS_FGMRES_PROJECTED_RELATIVE": "1",
    "LIS_MFAMG_REPLAY_METHOD": "fgcrodr",
    "LIS_MFAMG_REPLAY_RECYCLE_DIM": "8",
    "LIS_FGCRODR_INTRA_SOLVE_RECYCLE": "1",
    "LIS_FGCRODR_INTRA_TRACE": "1",
    "LIS_MFAMG_REPLAY_VALIDATE_RECYCLE": "1",
    "LIS_MFAMG_REPLAY_RESTART": "180",
    "LIS_MFAMG_REPLAY_MAXITER": "3000",
}
FIELDS = ("pair", "sequence", "variant", "n", "nnz", "modes", "process_exit", "solve_rc", "status",
          "iter", "true_relative", "wall_seconds", "builds", "updates",
          "uses", "rejects", "orth_error", "relation_error",
          "converged", "log_file")



def validate_omit01_basis(path17, path18, baseline=None):
    """Fail closed: the ONLY variable is deletion of Ritz column 01.

    The first 12 columns must match; the remaining five 17-mode columns
    must equal all-18 Ritz columns 02..06, in exactly the same order.
    """
    n17, m17 = read_modes(path17)
    n18, m18 = read_modes(path18)
    if n17 != n18 or len(m17) != 17 or len(m18) != 18:
        raise ValueError("expected 17 and 18 Matrix Market columns on identical row counts")
    if m17 != m18[:12] + m18[13:]:
        raise ValueError("17-mode input is not rotated 18-mode Ritz basis with direction 01 omitted")
    if baseline is not None:
        n0, m0 = read_modes(baseline)
        if n0 != n18 or len(m0) != 12 or m0 != m18[:12]:
            raise ValueError("12 baseline columns do not match the Ritz input")
    return n18


def fingerprint(path):
    """Local-only content digest, not an upload or an input-data export."""
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return {"sha256": h.hexdigest(), "bytes": path.stat().st_size}


def fingerprint_inputs(paths):
    return {label: fingerprint(path) for label, path in paths.items()}


def assert_inputs_unchanged(paths, expected):
    if fingerprint_inputs(paths) != expected:
        raise ValueError("input/EXE changed during benchmark; outputs are not comparable")


def kv(line):
    result = {}
    for word in line.split()[1:]:
        if "=" in word:
            key, val = word.split("=", 1)
            if key in result:
                raise ValueError("duplicate field " + key)
            result[key] = val
    return result


def integer(data, key):
    value = data.get(key, "")
    if not value or not value.strip("+-").isdigit():
        raise ValueError("invalid/missing integer field " + key)
    return int(value)


def real(data, key):
    try:
        value = float(data[key])
    except (KeyError, ValueError, OverflowError) as exc:
        raise ValueError("invalid/missing floating field " + key) from exc
    if not math.isfinite(value):
        raise ValueError("nonfinite field " + key)
    return value


def parse(text, exit_code, tol, log, expected_modes):
    groups = {name: [] for name in ("REPLAY", "REPLAY_RESULT", "REPLAY_RECYCLE_STATS",
                                    "REPLAY_RECYCLE_CHECK")}
    for line in text.splitlines():
        for name in groups:
            if line.startswith(name + " "):
                groups[name].append(kv(line))
    if any(len(group) != 1 for group in groups.values()):
        raise ValueError("expected exactly one REPLAY header, REPLAY_RESULT, REPLAY_RECYCLE_STATS and REPLAY_RECYCLE_CHECK")
    hdr = groups["REPLAY"][0]
    n, nnz, modes = (integer(hdr, k) for k in ("n", "nnz", "modes"))
    if (n <= 0 or nnz <= 0 or modes != expected_modes or
            hdr.get("method") != "fgcrodr" or integer(hdr, "recycle_k") != 8 or
            integer(hdr, "restart") != 180 or integer(hdr, "maxiter") != 3000):
        raise ValueError("unexpected replay configuration (matrix/modes/solver)")
    rr = groups["REPLAY_RESULT"][0]
    st = groups["REPLAY_RECYCLE_STATS"][0]
    va = groups["REPLAY_RECYCLE_CHECK"][0]
    status, rc, iterations = (integer(rr, k) for k in ("status", "solve_rc", "iter"))
    residual, wall = real(rr, "true_relative"), real(rr, "wall_seconds")
    orth, relation = real(va, "orth_error"), real(va, "relation_error")
    if exit_code not in (0, 1) or (exit_code == 0 and status != 0):
        raise ValueError("inconsistent/unexpected process exit or status")
    builds, updates, uses, rejects = (
        integer(st, k) for k in ("builds", "updates", "uses", "rejects"))
    if min(builds, updates, uses, rejects) < 0:
        raise ValueError("negative recycling counter")
    # A zero-dimensional/unused space can trivially satisfy algebra checks:
    # it must never be counted as positive evidence of intra-solve recycling.
    converged = (exit_code == 0 and status == 0 and rc == 0 and iterations > 0
                 and 0 <= residual <= tol and wall > 0 and integer(va, "rc") == 0
                 and 0 <= orth <= 1e-8 and 0 <= relation <= 1e-8
                 and builds >= 1 and updates >= 1 and uses >= 1)
    return dict(n=n, nnz=nnz, modes=modes,
                process_exit=exit_code, solve_rc=rc, status=status, iter=iterations,
                true_relative=residual, wall_seconds=wall,
                builds=builds, updates=updates, uses=uses, rejects=rejects,
                orth_error=orth, relation_error=relation,
                converged=int(converged), log_file=log)


def assess(rows, expected_pairs):
    groups = {name: [r for r in rows if r["variant"] == name]
              for name in ("omit01_17", "all18")}
    complete = (len(rows) == 2 * expected_pairs
                and all(sum(r["pair"] == p and r["sequence"] in (1, 2)
                            and r["variant"] == v for r in rows) == 1
                        for p in range(1, expected_pairs + 1)
                        for v in ("omit01_17", "all18")))
    accepted = complete and all(row["converged"] for row in rows)
    metrics = {}
    for name, group in groups.items():
        metrics[name] = {
            "runs": len(group),
            "converged": sum(r["converged"] for r in group),
            "median_wall_seconds": statistics.median(r["wall_seconds"] for r in group) if accepted else None,
            "median_iterations": statistics.median(r["iter"] for r in group) if accepted else None,
            "wall_seconds_range": [min(r["wall_seconds"] for r in group),
                                   max(r["wall_seconds"] for r in group)] if group else None,
        }
    gain = None
    if accepted and metrics["all18"]["median_wall_seconds"] > 0:
        gain = 100.0 * (1.0 - metrics["omit01_17"]["median_wall_seconds"] /
                        metrics["all18"]["median_wall_seconds"])
    paired_reductions = None
    if accepted:
        paired_reductions = []
        for p in range(1, expected_pairs + 1):
            t17 = next(r["wall_seconds"] for r in rows
                       if r["pair"] == p and r["variant"] == "omit01_17")
            t18 = next(r["wall_seconds"] for r in rows
                       if r["pair"] == p and r["variant"] == "all18")
            paired_reductions.append(100.0 * (1.0 - t17 / t18))
    return {
        "paired_wall_reduction_percent_by_pair": paired_reductions,
        "paired_median_wall_reduction_percent": (
            statistics.median(paired_reductions) if paired_reductions is not None else None),
        "assessment": "MATCHED_CONVERGED" if accepted else "INCONCLUSIVE_NONCONVERGENCE",
        "medians": metrics,
        "omit01_vs_all18_percent_median_wall_reduction": gain,
        "caution": "Matched pairs are exploratory; no statistical speedup or Newton-sequence claim.",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for option in ("exe", "matrix", "rhs", "modes-17", "modes-18", "out-dir"):
        p.add_argument("--" + option, type=Path, required=True)
    p.add_argument("--baseline-modes", type=Path, help="Optional original 12-mode verification")
    p.add_argument("--pairs", type=int, default=3)
    p.add_argument("--tol", type=float, default=1e-5)
    p.add_argument("--require-converged", action="store_true")
    a = p.parse_args()
    if not (1 <= a.pairs <= 50 and math.isfinite(a.tol) and 0 < a.tol < 1):
        p.error("pairs must be 1..50, tol must be finite in (0,1)")
    try:
        for path in (a.exe, a.matrix, a.rhs, a.modes_17, a.modes_18):
            if not path.is_file() or not path.stat().st_size:
                raise ValueError("missing/empty input: " + str(path))
        if a.baseline_modes is not None and (
                not a.baseline_modes.is_file() or not a.baseline_modes.stat().st_size):
            raise ValueError("missing/empty baseline modes: " + str(a.baseline_modes))
        inputs = {"exe": a.exe, "matrix": a.matrix, "rhs": a.rhs,
                  "modes_17": a.modes_17, "modes_18": a.modes_18}
        if a.baseline_modes is not None:
            inputs["baseline_modes"] = a.baseline_modes
        if len({path.resolve() for path in inputs.values()}) != len(inputs):
            raise ValueError("all input and executable paths must be distinct")
        expected_n = validate_omit01_basis(a.modes_17, a.modes_18,
                                           a.baseline_modes)
        input_hashes = fingerprint_inputs(inputs)
        if a.out_dir.exists():
            raise ValueError("output directory already exists: " + str(a.out_dir))
        a.out_dir.mkdir(parents=True)
        os.chmod(a.out_dir, 0o700)
        manifest = {
            "schema": "mfamg-v8-matched-inputs-v1",
            "input_fingerprints": input_hashes,
            "mode_preflight": "omit01_verified_against_rotated_all18",
            "matrix_rows_from_modes": expected_n,
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "pairs": a.pairs,
            "tolerance": a.tol,
            "fixed_solver_environment": ENV,
            "note": "Local only; never upload private files or raw contact logs",
        }
        (a.out_dir / "manifest.json").write_text(
            json.dumps(manifest, sort_keys=True, indent=2) + "
", encoding="utf-8")
        # Do not inherit unrelated experimental LIS_* toggles from the
        # caller's shell; they can silently alter the compared solver.
        env = {key: value for key, value in os.environ.items()
               if not key.startswith("LIS_")}
        env.update(ENV)
        rows = []
        # Journal completed runs immediately: a late failing replay must not
        # discard evidence from earlier expensive, private local solves.
        with (a.out_dir / "matched_runs.csv").open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            stream.flush()
            for pair in range(1, a.pairs + 1):
                order = ("omit01_17", "all18") if pair % 2 else ("all18", "omit01_17")
                for sequence, variant in enumerate(order, 1):
                    modes = a.modes_17 if variant == "omit01_17" else a.modes_18
                    logfile = "pair%02d_run%d_%s.log" % (pair, sequence, variant)
                    result = subprocess.run(
                        [str(a.exe.resolve()), str(a.matrix.resolve()),
                         str(a.rhs.resolve()), str(modes.resolve())],
                        cwd=str(a.matrix.resolve().parent), env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        text=True, errors="replace", check=False)
                    (a.out_dir / logfile).write_text(result.stdout, encoding="utf-8")
                    assert_inputs_unchanged(inputs, input_hashes)
                    row = parse(result.stdout, result.returncode, a.tol, logfile,
                                17 if variant == "omit01_17" else 18)
                    if row["n"] != expected_n:
                        raise ValueError("replay row count disagrees with preflight mode files")
                    if rows and (row["n"], row["nnz"]) != (rows[0]["n"], rows[0]["nnz"]):
                        raise ValueError("matrix shape/nnz changed between paired runs")
                    row.update(pair=pair, sequence=sequence, variant=variant)
                    rows.append(row)
                    writer.writerow(row)
                    stream.flush()
                    print("V8_REPEAT pair=%d case=%s status=%d iter=%d wall=%.3f converged=%d" %
                          (pair, variant, row["status"], row["iter"],
                           row["wall_seconds"], row["converged"]), flush=True)
        assert_inputs_unchanged(inputs, input_hashes)
        report = assess(rows, a.pairs)
        report["manifest_file"] = "manifest.json"
        report["pairs"] = a.pairs
        report["fixed_solver_environment"] = ENV
        (a.out_dir / "assessment.json").write_text(
            json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        print("V8_REPEAT_ASSESSMENT=" + report["assessment"])
        print("V8_REPEAT_SUMMARY=" + str(a.out_dir / "assessment.json"))
        return 1 if a.require_converged and report["assessment"] != "MATCHED_CONVERGED" else 0
    except (OSError, ValueError) as exc:
        print("V8_REPEAT_ERROR=" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
