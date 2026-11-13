# SuperLU ILU equilibration / ordering ablation (experimental)

This is an **offline control experiment**, not a change to LIS's ILUT.
It tests whether SuperLU's **Equil** option adds measurable value when all
other settings, especially pivoting, remain fixed. It also measures the
independent effect of `COLAMD` versus `NATURAL` column ordering.

## Inputs

Use **actual assembled systems** exported with
[`lis_snapshot_export_mm()`](lis_snapshot_export.md), e.g.
`snapshots/step03_iter07_A.mtx` and
`snapshots/step03_iter07_b.mtx`. Matrix/RHS pairs with names
`A.mtx`/`b.mtx` in a common directory are also recognized. Alternatively,
supply explicit `--rhs PATH` arguments (once per matrix, in matching order).

There are no five CCX contact matrices in this repository. Their successful
5/5 SuperLU outcome cannot be decomposed into individual causes until those
same snapshots (including RHS) are provided and evaluated.

## Run

Requirements: Python 3, NumPy and SciPy with `spilu` and the modern
`gmres(rtol=..., callback_type=...)` API.

```sh
python3 tools/compare_superlu_ilut_equil.py \
  snapshots/*_A.mtx \
  --rcm --row-max-scale \
  --drop-tol 1e-4 --fill-factor 10 \
  --rtol 1e-8 --restart 40 --maxiter 200 \
  --csv superlu_equil_ablation.csv \
  --json superlu_equil_ablation.json
```

Match `--drop-tol`, `--fill-factor`, `--restart`,
`--maxiter`, `--rtol`, and any `--diag-pivot-thresh` to the **original**
successful run. Omit `--diag-pivot-thresh` to keep SciPy's SuperLU default
unchanged across all four cases. `--rcm` applies a symmetric permutation
to the rows and columns; `--row-max-scale` applies a fixed external
row-max-norm scaling. Both are held constant **between variants**.

Per matrix, the script tests exactly these four combinations:

| `permc_spec` | `options['Equil']` |
| --- | --- |
| NATURAL | False |
| NATURAL | True |
| COLAMD | False |
| COLAMD | True |

All variants use the **same** `spilu` drop tolerance, fill factor,
SuperLU row-pivoting policy and left-preconditioned GMRES configuration.
The original RHS is reordered/scaled with the matrix. The solution is
unpermuted and validated against the **original**, unscaled `A` and `b`.
The CSV records setup and solve time, GMRES iterations, full-system relative
residual `||b - A*x||_2 / ||b||_2`, factor NNZ, SciPy's
`perm_r` nonidentity count, success/failure and any exception. Row
permutation count is **not** a pivoting-on/off experiment.

If a real RHS is unavailable, `--manufactured-rhs` explicitly substitutes
`b = A @ random_vector(seed=42)`. Such a run is a **smoke test only**, not
a substitute for reproducing contact-system convergence.

## Interpret the results

Compare Equil True/False **within each ordering**; compare NATURAL/COLAMD
**within each Equil setting**. Treat only runs with acceptable **original
physical residuals** as converged. Look for consistently rescued cases or
substantial iteration reductions without regressions, not just lower
factorization time. Four cases per matrix means **20 rows for five snapshots**.

The tool is not a pivoting ablation: pivoting stays enabled in all cases.
If Equil adds little, the previously reported 5/5 result may depend on
ordering, pivoting or interactions between them. If Equil does help, it
justifies a **separate**, default-off `-ilut_equil 1` implementation
proposal for LIS with tests for the regular and adjoint apply paths,
FP16/quad, OpenMP/MPI and frozen-pattern preconditioner refresh.
SciPy/SuperLU ILU and LIS ILUT are **different algorithms**, so even a
clear SuperLU gain is not a proof that the LIS port will match it.

The LIS `row_scale` variable currently used in
`src/precon/lis_precon_ilut.c` controls *pivot regularization*; it does
not implement two-sided equilibration. The top-level LIS `-scale` option
is likewise a different solver scaling mechanism.

## Smoke tests

```sh
python3 -m unittest discover -s tools -p 'test_compare_superlu_ilut_equil.py' -v
python3 tools/compare_superlu_ilut_equil.py \
  test/testmat.mtx --manufactured-rhs --rcm --row-max-scale \
  --csv /tmp/superlu_equil_smoke.csv
```

The smoke test is intended to validate the experiment's mechanics, **not**
to demonstrate a numerical advantage for Equil.
