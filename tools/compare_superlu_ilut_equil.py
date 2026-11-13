#!/usr/bin/env python3
"""Controlled SciPy/SuperLU ILU ablation: Equil x column ordering.

Only Equil and permc_spec vary across the four cases per matrix.  This is
an experiment, not a substitute for LIS's own ILUT implementation.
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmread
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse.linalg import LinearOperator, gmres, spilu


FIELDS = (
    "matrix", "rhs_source", "n", "nnz", "rcm", "row_max_scale",
    "equil", "permc_spec", "status", "info", "iterations", "setup_s",
    "solve_s", "rel_residual", "L_nnz", "U_nnz", "row_perm_moved",
    "error",
)


def read_matrix(path):
    a = mmread(path)
    if not sparse.issparse(a):
        a = sparse.coo_matrix(a)
    a = a.tocsc().astype(np.result_type(a.dtype, np.float64))
    a.sum_duplicates()
    a.eliminate_zeros()
    if a.shape[0] != a.shape[1]:
        raise ValueError(f"{path}: square matrix required, got {a.shape}")
    if a.shape[0] == 0 or not np.all(np.isfinite(a.data)):
        raise ValueError(f"{path}: empty or non-finite matrix")
    return a


def read_rhs(path, n):
    if path.lower().endswith(".npy"):
        b = np.load(path, allow_pickle=False)
    else:
        b = mmread(path)
        if sparse.issparse(b):
            b = b.toarray()
    b = np.asarray(b).reshape(-1)
    if b.size != n or not np.all(np.isfinite(b)):
        raise ValueError(f"{path}: RHS must have {n} finite entries")
    return b


def prepare_system(a, b, rcm, row_max_scale):
    n = a.shape[0]
    perm = (reverse_cuthill_mckee(a.tocsr(), symmetric_mode=False)
            if rcm else np.arange(n))
    work = a[perm, :][:, perm].tocsc()
    rhs = b[perm].copy()
    if row_max_scale:
        maxabs = np.abs(work).max(axis=1).toarray().reshape(-1)
        invscale = np.ones(n, dtype=np.float64)
        nonzero = maxabs > 0
        invscale[nonzero] = 1.0 / maxabs[nonzero]
        if not np.all(np.isfinite(invscale)):
            raise ValueError("row max-norm scaling overflowed")
        work = (sparse.diags(invscale, format="csc") @ work).tocsc()
        rhs *= invscale
    return work, rhs, perm


def run_case(a, b, work, rhs, perm, args, equil, ordering, rhs_source):
    result = dict(
        matrix=str(args.current_matrix), rhs_source=rhs_source,
        n=a.shape[0], nnz=a.nnz, rcm=args.rcm,
        row_max_scale=args.row_max_scale, equil=equil, permc_spec=ordering,
        status="failed", info="", iterations=0, setup_s="", solve_s="",
        rel_residual="", L_nnz="", U_nnz="", row_perm_moved="", error="",
    )
    try:
        settings = {"Equil": bool(equil)}
        kwargs = dict(
            drop_tol=args.drop_tol, fill_factor=args.fill_factor,
            permc_spec=ordering, options=settings,
        )
        if args.diag_pivot_thresh is not None:
            kwargs["diag_pivot_thresh"] = args.diag_pivot_thresh
        t0 = time.perf_counter()
        ilu = spilu(work, **kwargs)
        result["setup_s"] = time.perf_counter() - t0
        result["L_nnz"] = ilu.L.nnz
        result["U_nnz"] = ilu.U.nnz
        result["row_perm_moved"] = int(np.count_nonzero(
            ilu.perm_r != np.arange(a.shape[0])
        ))
        precon = LinearOperator(work.shape, matvec=ilu.solve, dtype=work.dtype)
        iteration_count = [0]

        def record_iteration(_):
            iteration_count[0] += 1

        t0 = time.perf_counter()
        x_perm, info = gmres(
            work, rhs, M=precon, atol=args.atol, rtol=args.rtol,
            restart=args.restart, maxiter=args.maxiter,
            callback=record_iteration, callback_type="pr_norm",
        )
        result["solve_s"] = time.perf_counter() - t0
        result["iterations"] = iteration_count[0]
        result["info"] = int(info)
        x = np.empty_like(x_perm)
        x[perm] = x_perm
        rel_res = np.linalg.norm(b - a @ x) / np.linalg.norm(b)
        result["rel_residual"] = float(rel_res)
        result["status"] = "converged" if info == 0 and np.isfinite(rel_res) and rel_res <= args.rtol * args.residual_allowance else "not_converged"
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrices", nargs="+", type=Path, help="Matrix Market A files")
    parser.add_argument("--rhs", action="append", type=Path, default=None,
                        help="explicit RHS path; repeat once per matrix")
    parser.add_argument("--manufactured-rhs", action="store_true",
                        help="explicitly use A @ random(seed=42) instead of physical RHS")
    parser.add_argument("--rcm", action="store_true", help="same symmetric RCM on rows and columns for all four cases")
    parser.add_argument("--row-max-scale", action="store_true", help="same external max-norm row scaling in all cases")
    parser.add_argument("--drop-tol", type=float, default=1e-4)
    parser.add_argument("--fill-factor", type=float, default=10.0)
    parser.add_argument("--diag-pivot-thresh", type=float, default=None,
                        help="fix SuperLU pivoting threshold, if needed")
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--restart", type=int, default=40)
    parser.add_argument("--maxiter", type=int, default=200, help="maximum GMRES restart cycles")
    parser.add_argument("--residual-allowance", type=float, default=10.,
                        help="physical relative residual must be <= rtol * this factor")
    parser.add_argument("--csv", type=Path, default=Path("superlu_equil_ablation.csv"))
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args(argv)
    if args.rhs is not None and len(args.rhs) != len(args.matrices):
        parser.error("repeat --rhs once per matrix (or omit for auto-detection)")
    if args.rhs is not None and args.manufactured_rhs:
        parser.error("--rhs and --manufactured-rhs are mutually exclusive")
    if args.drop_tol < 0 or args.fill_factor < 1 or args.rtol <= 0 or args.atol < 0 or args.restart < 1 or args.maxiter < 1 or args.residual_allowance <= 0:
        parser.error("invalid drop tolerance, fill, solver tolerance or iteration limit")
    if args.diag_pivot_thresh is not None and not 0 <= args.diag_pivot_thresh <= 1:
        parser.error("--diag-pivot-thresh must be within [0, 1]")
    return args


def main(argv=None):
    args = parse_args(argv)
    results = []
    for i, matrix_path in enumerate(args.matrices):
        args.current_matrix = matrix_path
        a = read_matrix(matrix_path)
        if args.manufactured_rhs:
            rng = np.random.default_rng(42)
            b = a @ rng.standard_normal(a.shape[0])
            source = "manufactured(seed=42)"
        else:
            if args.rhs is not None:
                rhs_path = args.rhs[i]
            elif matrix_path.name.endswith("_A.mtx"):
                rhs_path = matrix_path.with_name(matrix_path.name[:-6] + "_b.mtx")
            elif matrix_path.name == "A.mtx":
                rhs_path = matrix_path.with_name("b.mtx")
            else:
                raise ValueError(f"{matrix_path}: no predictable RHS; pass --rhs or --manufactured-rhs")
            if not rhs_path.is_file():
                raise ValueError(f"{rhs_path}: missing RHS; pass --rhs or --manufactured-rhs")
            b = read_rhs(str(rhs_path), a.shape[0]).astype(a.dtype, copy=False)
            source = str(rhs_path)
        if not np.linalg.norm(b):
            raise ValueError(f"{matrix_path}: zero RHS; provide --rhs")
        work, rhs, perm = prepare_system(a, b, args.rcm, args.row_max_scale)
        for ordering in ("NATURAL", "COLAMD"):
            for equil in (False, True):
                row = run_case(a, b, work, rhs, perm, args, equil, ordering, source)
                results.append(row)
                print(f"{matrix_path.name} {ordering:7} Equil={int(equil)} "
                      f"{row['status']:13} its={row['iterations']} "
                      f"relres={row['rel_residual']} {row['error']}")
    with args.csv.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(results)
    if args.json is not None:
        args.json.write_text(json.dumps({
            "method": "SciPy spilu + left-preconditioned GMRES; two-factor Equil/permc_spec ablation",
            "fixed": {"rcm": args.rcm, "row_max_scale": args.row_max_scale,
                      "drop_tol": args.drop_tol, "fill_factor": args.fill_factor,
                      "diag_pivot_thresh": args.diag_pivot_thresh,
                      "rtol": args.rtol, "atol": args.atol,
                      "restart": args.restart, "maxiter_cycles": args.maxiter},
            "results": results,
        }, indent=2), encoding="utf-8")
    print(f"Wrote {len(results)} cases to {args.csv}")
    if any(row["status"] != "converged" for row in results):
        print("NOTE: some variants did not converge; inspect per-case status and physical residual", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
