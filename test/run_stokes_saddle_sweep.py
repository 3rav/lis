#!/usr/bin/env python3
"""Compare ILUT, Equil-only, ILUTP-only, and ILUTP+Equil on Stokes saddles.

Uses native ./test/benchilutp, preserves the generated original A and b,
checks the *physical* residual, and compares actual stored L+U fill.
Do not label a result a synergy merely because it beats unscaled ILUT.
"""
import argparse
import csv
import itertools
import json
import math
import subprocess
import tempfile
from pathlib import Path

from gen_stokes_saddle import make

BENCHMARK_MODES = (
    ("plain", 0, 0),
    ("row_only", 1, 0),
    ("equil", 0, 1),
)


def converged(row, limit):
    try:
        return (row["status"] == "converged"
                and int(row["api_error"]) == 0
                and int(row["solver_status"]) == 0
                and float(row["physical_relres"]) <= limit)
    except (KeyError, ValueError):
        return False


def run_bench(bench, matrix, rhs, output, args, rate, drop, row_scale, equil):
    cmd = [str(bench), str(matrix), str(rhs),
           "--maxiter", str(args.maxiter), "--restart", str(args.restart),
           "--rtol", str(args.rtol), "--drop", str(drop), "--rate", str(rate),
           "--pivot-tol", str(args.pivot_tol),
           "--rcm", "0", "--row-scale", str(row_scale),
           "--equil", str(equil), "--csv", str(output)]
    job = subprocess.run(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True,
                         timeout=args.timeout)
    if job.returncode != 0 or not output.is_file():
        raise RuntimeError(
            f"benchilutp exited {job.returncode}: {job.stderr[-1000:]}")
    with output.open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    if len(rows) != 3:
        raise RuntimeError(f"expected 3 variants; got {len(rows)}")
    for row in rows:
        if row.get("equil") != str(equil):
            raise RuntimeError(f"requested Equil={equil} but row={row}")
    return rows


def evaluate_variant(rows, limit):
    by_key = {(row["scaling"], row["variant"]): row for row in rows}
    standalone_pivot = by_key.get(("plain", "ilutp_0p5"))
    standalone_equil = by_key.get(("equil", "legacy_ilut"))
    combined = by_key.get(("equil", "ilutp_0p5"))
    if not (standalone_pivot and standalone_equil and combined):
        return None

    def factor_nnz(row):
        return int(row["nnz_L"]) + int(row["nnz_U"])

    p_ok = converged(standalone_pivot, limit)
    e_ok = converged(standalone_equil, limit)
    c_ok = converged(combined, limit)
    answer = {
        "pivot_ok": p_ok, "equil_ok": e_ok, "combined_ok": c_ok,
        "pivot_iter": int(standalone_pivot["iterations"]),
        "equil_iter": int(standalone_equil["iterations"]),
        "combined_iter": int(combined["iterations"]),
        "pivot_fill": factor_nnz(standalone_pivot),
        "equil_fill": factor_nnz(standalone_equil),
        "combined_fill": factor_nnz(combined),
        "combined_res": float(combined["physical_relres"]),
    }
    answer["strict_synergy"] = c_ok and not p_ok and not e_ok
    answer["both_faster"] = (
        c_ok and p_ok and e_ok
        and answer["combined_iter"] < min(answer["pivot_iter"],
                                          answer["equil_iter"]))
    answer["bounded_fill"] = (
        answer["combined_fill"] <=
        1.20 * max(answer["pivot_fill"], answer["equil_fill"]))
    return answer


def run(args):
    bench = args.bench.resolve()
    if not bench.exists():
        raise FileNotFoundError(bench)
    args.output.mkdir(parents=True, exist_ok=True)
    samples = [
        (side, contrast, units, ordering, gamma, rate, drop)
        for side, contrast, units, ordering, gamma, rate, drop
        in itertools.product(
            args.sides, args.contrasts, args.units, args.orderings,
            args.gammas, args.rates, args.drops)
    ]
    if args.quick:
        samples = [sample for sample in samples
                   if sample[0] == args.sides[0]
                   and sample[1] == args.contrasts[0]
                   and sample[4] == args.gammas[0]
                   and sample[5] == args.rates[0]
                   and sample[6] == args.drops[0]]
    collected = []
    reports = []
    with tempfile.TemporaryDirectory(prefix="lis-stokes-saddle-") as temp:
        work = Path(temp)
        last_matrix = None
        for idx, (side, contrast, units, ordering, gamma, rate, drop) in enumerate(samples):
            matrix, rhs = work / "A.mtx", work / "b.mtx"
            descriptor = (side, contrast, units, ordering, gamma)
            if descriptor != last_matrix:
                meta = make(side, contrast, units, ordering, 0.0, gamma, 4.5,
                            matrix, rhs)
                last_matrix = descriptor
            else:
                meta = {"n": 3 * side * side, "nnz": None}
            rows = []
            for scaling, row_scale, equil in BENCHMARK_MODES:
                filename = work / f"{scaling}.csv"
                try:
                    result = run_bench(bench, matrix, rhs, filename, args,
                                       rate, drop, row_scale, equil)
                except (RuntimeError, subprocess.TimeoutExpired) as exc:
                    result = []
                    print(f"SADDLE_ERROR idx={idx} scaling={scaling}: {exc}",
                          flush=True)
                for row in result:
                    row.update({
                        "scaling": scaling,
                        "side": side, "n": 3 * side * side,
                        "contrast": contrast, "units": units,
                        "ordering": ordering, "gamma": gamma,
                        "rate": rate, "drop": drop, "case_id": idx,
                    })
                rows.extend(result)
            collected.extend(rows)
            result = evaluate_variant(rows, args.accept_residual)
            if result is None:
                continue
            result.update({
                "case_id": idx, "side": side, "n": 3 * side * side,
                "contrast": contrast, "units": units, "ordering": ordering,
                "gamma": gamma, "rate": rate, "drop": drop})
            reports.append(result)
            label = ("STRICT" if result["strict_synergy"] else
                     "FASTER" if result["both_faster"] else "no")
            if label != "no":
                print("SADDLE_CANDIDATE " + json.dumps(result, sort_keys=True),
                      flush=True)
                if result["strict_synergy"] or (
                        result["both_faster"] and result["bounded_fill"]):
                    fixture = args.output / f"positive_{idx}"
                    fixture.mkdir(parents=True, exist_ok=True)
                    (fixture / "A.mtx").write_bytes(matrix.read_bytes())
                    (fixture / "b.mtx").write_bytes(rhs.read_bytes())
            if idx % 20 == 0:
                print(f"SADDLE_PROGRESS {idx + 1}/{len(samples)}", flush=True)

    def write_csv(path, data):
        if not data:
            path.write_text("", encoding="utf-8")
            return
        keys = list(dict.fromkeys(key for row in data for key in row.keys()))
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(data)
    write_csv(args.output / "all_native_variants.csv", collected)
    write_csv(args.output / "all_saddle_comparisons.csv", reports)
    summary = {
        "tested_settings": len(samples),
        "completed_comparisons": len(reports),
        "strict_synergy": sum(x["strict_synergy"] for x in reports),
        "faster_than_both": sum(x["both_faster"] for x in reports),
        "faster_with_bounded_fill": sum(
            x["both_faster"] and x["bounded_fill"] for x in reports),
        "combined_converged": sum(x["combined_ok"] for x in reports),
        "one_or_more_controls_failed": sum(
            x["combined_ok"] and not(x["pivot_ok"] and x["equil_ok"])
            for x in reports),
        "top_candidates": sorted(
            (r for r in reports if r["both_faster"] or r["strict_synergy"]),
            key=lambda r: (
                -int(r["strict_synergy"]),
                -int(r["both_faster"] and r["bounded_fill"]),
                r["combined_iter"],
                r["combined_fill"]),
        )[:20]
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print("SADDLE_SUMMARY " + json.dumps(
        {k: v for k, v in summary.items() if k != "top_candidates"},
        sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", type=Path, default=Path("./test/benchilutp"))
    parser.add_argument("--output", type=Path, default=Path("stokes-results"))
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--sides", type=int, nargs="+", default=[4, 6])
    parser.add_argument("--contrasts", type=float, nargs="+", default=[30, 300])
    parser.add_argument("--units", nargs="+", default=[
        "pressure_checker", "pressure_uniform", "pressure_weak", "mixed_checker"])
    parser.add_argument("--orderings", nargs="+", default=[
        "pressure_first", "interleaved"])
    parser.add_argument("--gammas", type=float, nargs="+", default=[0.3, 1])
    parser.add_argument("--rates", type=float, nargs="+", default=[1, 3])
    parser.add_argument("--drops", type=float, nargs="+", default=[0.005, 0.05])
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--accept-residual", type=float, default=1e-5)
    parser.add_argument("--maxiter", type=int, default=140)
    parser.add_argument("--restart", type=int, default=50)
    parser.add_argument("--pivot-tol", type=float, default=0.001)
    parser.add_argument("--timeout", type=float, default=90)
    run(parser.parse_args())
