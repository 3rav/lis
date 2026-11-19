#!/usr/bin/env python3
"""Deterministic public ILUTP top-k A/B benchmark. Python stdlib only."""
import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
from datetime import datetime, timezone

FIELDS = (
    "variant", "ilutp", "equil", "tau", "status", "iterations",
    "nnz_L", "nnz_U", "moved_columns", "setup_cpu_s", "solve_cpu_s",
    "physical_relres", "api_error", "solver_status",
)
VARIANTS = ("legacy_ilut", "ilutp_0p5", "ilutp_1p0")
BASE_SHA = "3d027903ab079158d99bfc91f8ab6867c212dc20"
PATCH_SHA = "78e694a82fe15575b22df11977b7277617e9b039"

def values_for_row(i, n, far, seed):
    """Fixed nonsymmetric sparse pattern with strict row diagonal dominance."""
    columns = set(range(max(0, i - 3), min(n, i + 4)))
    columns.discard(i)
    for k in range(far):
        j = (i * 1103515245 + k * 2654435761 + seed * 97) % n
        while j in columns or j == i:
            j = (j + 1) % n
        columns.add(j)
    off = []
    for j in sorted(columns):
        strength = 0.025 + ((i * 43 + j * 71 + seed) % 53) / 1000.0
        v = -strength if (j < i or (i + j) % 3 == 0) else strength
        off.append((j, v))
    diagonal = 1.0 + math.fsum(abs(v) for _, v in off)
    return sorted(off + [(i, diagonal)])

def generate_inputs(n, far, seed, target):
    if n < 128 or not 0 <= far <= 32 or n <= 2 * far + 7:
        raise ValueError("require n>=128, 0<=far<=32, n>2*far+7")
    target.mkdir(parents=True, exist_ok=True)
    a_path, b_path = target / "synthetic_A.mtx", target / "synthetic_b.mtx"
    rows = [values_for_row(i, n, far, seed) for i in range(n)]
    nnz = sum(map(len, rows))
    with a_path.open("w", encoding="ascii", newline="\n") as a, \
         b_path.open("w", encoding="ascii", newline="\n") as b:
        a.write("%%MatrixMarket matrix coordinate real general\n"
                "% deterministic synthetic sparse system, NOT CalculiX data\n"
                f"{n} {n} {nnz}\n")
        b.write("%%MatrixMarket matrix array real general\n"
                "% manufactured b=A*x, x[j]=1+((17*j)%97)/100\n"
                f"{n} 1\n")
        for i, row in enumerate(rows):
            for j, v in row:
                a.write(f"{i+1} {j+1} {v:.17g}\n")
            rhs = math.fsum(v * (1 + ((17 * j) % 97) / 100.0)
                            for j, v in row)
            b.write(f"{rhs:.17g}\n")
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    return a_path, b_path, {"n": n, "nnz": nnz, "far_per_row": far,
                            "seed": seed, "matrix_sha256": sha(a_path),
                            "rhs_sha256": sha(b_path)}

def invoke(executable, matrix, rhs, rate, run, order, label, out, config):
    cmd = [str(executable), str(matrix), str(rhs),
           "--rate", str(rate), "--maxiter", str(config["maxiter"]),
           "--restart", str(config["restart"]),
           "--rtol", str(config["rtol"]), "--drop", str(config["drop"]),
           "--pivot-tol", str(config["pivot_tol"]),
           "--rcm", "0", "--row-scale", "0", "--equil", "0"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", LC_ALL="C")
    env["LD_LIBRARY_PATH"] = (str(executable.parent.parent / "src" / ".libs")
                              + os.pathsep + env.get("LD_LIBRARY_PATH", ""))
    p = subprocess.run(cmd, text=True, capture_output=True, env=env, check=False)
    tag = f"rate{rate}-run{run}-{order}-{label}"
    (out / f"{tag}.stderr.log").write_text(p.stderr, encoding="utf-8")
    (out / f"{tag}.stdout.csv").write_text(p.stdout, encoding="utf-8")
    if p.returncode != 0:
        raise RuntimeError(f"benchmark {tag} exited {p.returncode}; see logs")
    # LIS prints its matrix-size banner before the benchmark CSV header.
    lines = p.stdout.splitlines(keepends=True)
    header = ",".join(FIELDS)
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == header)
    except StopIteration:
        raise RuntimeError(f"benchmark {tag} emitted no expected CSV header; see logs")
    rd = csv.DictReader(io.StringIO("".join(lines[start:])))
    if rd.fieldnames != list(FIELDS):
        raise RuntimeError(f"unexpected CSV schema for {tag}: {rd.fieldnames}")
    result = list(rd)
    if [r["variant"] for r in result] != list(VARIANTS):
        raise RuntimeError(f"unexpected variants in {tag}")
    return [{"build": label, "rate": rate, "run": run, "order": order,
             **r} for r in result]

def check_pair(baseline, candidate, rate, run):
    errors = []
    a = {r["variant"]: r for r in baseline}
    b = {r["variant"]: r for r in candidate}
    for variant in VARIANTS:
        left, right = a[variant], b[variant]
        for field in ("status", "iterations", "nnz_L", "nnz_U",
                      "moved_columns", "api_error", "solver_status"):
            if left[field] != right[field]:
                errors.append(f"rate={rate} run={run} {variant} {field}: "
                              f"{left[field]} != {right[field]}")
        for label, record in (("baseline", left), ("candidate", right)):
            if (record["status"] != "converged" or
                    record["api_error"] != "0"):
                errors.append(f"rate={rate} run={run} {variant} "
                              f"{label} not converged: {record['status']}")
            try:
                resid = float(record["physical_relres"])
                if not math.isfinite(resid) or resid > 1e-7 or resid < 0:
                    errors.append(f"rate={rate} run={run} {variant} "
                                  f"{label} invalid residual {resid}")
            except ValueError:
                errors.append(f"rate={rate} run={run} {variant} "
                              f"{label} invalid residual text")
        try:
            ra, rb = float(left["physical_relres"]), float(right["physical_relres"])
            if not math.isclose(ra, rb, rel_tol=1e-8, abs_tol=5e-12):
                errors.append(f"rate={rate} run={run} {variant} "
                              f"residual mismatch {ra} vs {rb}")
        except ValueError:
            pass
    return errors

def summarize(rows):
    summary = []
    for rate in sorted({int(r["rate"]) for r in rows}):
        for variant in VARIANTS:
            vals = {}
            for build in ("baseline", "candidate"):
                chosen = [r for r in rows if r["rate"] == rate and
                          r["variant"] == variant and r["build"] == build]
                vals[build] = {
                    field: statistics.median(float(r[field]) for r in chosen)
                    for field in ("setup_cpu_s", "solve_cpu_s")}
            a, b = vals["baseline"], vals["candidate"]
            # Distributions use matched A/B runs, not ratios of independently
            # chosen medians. The legacy ILUT row is a negative timing control.
            by_build = {
                name: {int(r["run"]): r for r in rows
                       if r["rate"] == rate and r["variant"] == variant
                       and r["build"] == name}
                for name in ("baseline", "candidate")
            }
            runs = sorted(set(by_build["baseline"]) & set(by_build["candidate"]))
            def paired_changes(field):
                result = []
                for run in runs:
                    baseline = float(by_build["baseline"][run][field])
                    candidate = float(by_build["candidate"][run][field])
                    if baseline <= 0:
                        raise ValueError(f"nonpositive {field}: rate={rate} run={run}")
                    result.append((candidate / baseline - 1) * 100)
                return result
            setup_pairs = paired_changes("setup_cpu_s")
            solve_pairs = paired_changes("solve_cpu_s")
            summary.append({"rate": rate, "variant": variant,
                            "paired_runs": len(runs),
                            "baseline_setup_median_s": a["setup_cpu_s"],
                            "candidate_setup_median_s": b["setup_cpu_s"],
                            "setup_delta_percent": ((b["setup_cpu_s"] / a["setup_cpu_s"] - 1) * 100
                                                    if a["setup_cpu_s"] > 0 else ""),
                            "paired_setup_delta_median_percent": statistics.median(setup_pairs),
                            "paired_setup_delta_min_percent": min(setup_pairs),
                            "paired_setup_delta_max_percent": max(setup_pairs),
                            "baseline_solve_median_s": a["solve_cpu_s"],
                            "candidate_solve_median_s": b["solve_cpu_s"],
                            "solve_delta_percent": ((b["solve_cpu_s"] / a["solve_cpu_s"] - 1) * 100
                                                    if a["solve_cpu_s"] > 0 else ""),
                            "paired_solve_delta_median_percent": statistics.median(solve_pairs),
                            "paired_solve_delta_min_percent": min(solve_pairs),
                            "paired_solve_delta_max_percent": max(solve_pairs)})
    return summary

def write_csv(filename, records):
    if not records:
        return
    with filename.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline-bin", type=Path, required=True)
    p.add_argument("--candidate-bin", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("ilutp-topk-results"))
    p.add_argument("--n", type=int, default=4000)
    p.add_argument("--far", type=int, default=8)
    p.add_argument("--seed", type=int, default=20261119)
    p.add_argument("--rates", default="1,2,4,8")
    p.add_argument("--repeats", type=int, default=7)
    p.add_argument("--warmups", type=int, default=1)
    p.add_argument("--baseline-ref", default=BASE_SHA)
    p.add_argument("--candidate-ref", default=PATCH_SHA)
    args = p.parse_args()
    if args.repeats < 1 or args.warmups < 0:
        p.error("repeats must be positive and warmups nonnegative")
    try:
        rates = [int(x) for x in args.rates.split(",")]
        if not rates or len(rates) != len(set(rates)) or any(x < 1 or x > 8 for x in rates):
            raise ValueError
    except ValueError:
        p.error("rates must be unique integers in [1,8], comma-separated")
    binaries = {}
    for label, path in (("baseline", args.baseline_bin),
                        ("candidate", args.candidate_bin)):
        binary = path.resolve()
        if not binary.is_file() or not os.access(binary, os.X_OK):
            p.error(f"{label} binary not executable: {binary}")
        binaries[label] = binary
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    a, b, data_manifest = generate_inputs(args.n, args.far, args.seed, out / "inputs")
    config = {"maxiter": 220, "restart": 50, "rtol": 1e-8,
              "drop": 0, "pivot_tol": 1e-9}
    manifest = {"utc_started": datetime.now(timezone.utc).isoformat(),
                "host": platform.platform(), "python": sys.version,
                "cpu": platform.processor(), "baseline_ref": args.baseline_ref,
                "candidate_ref": args.candidate_ref,
                "binaries": {k: str(v) for k, v in binaries.items()},
                "input": data_manifest, "solver_options": config,
                "rates": rates, "repeats": args.repeats, "warmups": args.warmups,
                "measurement": "per-process CPU time (clock), not wall time",
                "policy": "numeric parity gates only; speedup never gates CI"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows, errors = [], []
    try:
        for rate in rates:
            for rep in range(-args.warmups, args.repeats):
                labels = ("baseline", "candidate") if rep % 2 == 0 else ("candidate", "baseline")
                pair = {}
                for order, label in enumerate(labels):
                    chunk = invoke(binaries[label], a, b, rate, rep, order,
                                   label, out, config)
                    pair[label] = chunk
                    if rep >= 0:
                        rows.extend(chunk)
                if rep >= 0:
                    errors.extend(check_pair(pair["baseline"], pair["candidate"], rate, rep))
    except Exception as exc:
        errors.append(f"runner failure: {exc}")
    write_csv(out / "raw.csv", rows)
    if rows and all(sum(r["build"] == build and r["rate"] == rate and
                        r["variant"] == variant for r in rows) == args.repeats
                    for rate in rates for variant in VARIANTS
                    for build in ("baseline", "candidate")):
        summary = summarize(rows)
        write_csv(out / "summary.csv", summary)
        print("rate,variant,baseline_setup_s,candidate_setup_s,delta_percent")
        for r in summary:
            print(f"{r['rate']},{r['variant']},{r['baseline_setup_median_s']:.6g},"
                  f"{r['candidate_setup_median_s']:.6g},"
                  f"{r['setup_delta_percent']:.3f}")
    (out / "validation.txt").write_text(
        ("PASS: identical status, iterations, factor counts and close true residuals\n"
         if not errors else "FAIL:\n" + "\n".join(errors) + "\n"))
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Numerical parity PASS. Results: {out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
