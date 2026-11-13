# SuperLU ILU max-norm equilibration / ordering ablation

This is an **offline control experiment**, not a change to LIS ILUT or the
CCX adapter. It evaluates the effect of explicit max-norm two-sided scaling
`R*A*C` versus the existing fixed row scaling, and `NATURAL` versus
`COLAMD`, with all other SuperLU/GMRES settings held constant.

## Critical correction: SciPy `spilu(options={"Equil": ...})`

**Simply toggling the `Equil` option passed to SciPy `spilu` is NOT an
equilibration experiment.** SciPy 1.17's
[`_superluobject.c`](https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/linalg/_dsolve/_superluobject.c)
passes the option to the low-level
[`dgsitrf`](https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/linalg/_dsolve/SuperLU/SRC/dgsitrf.c)
factorization path, which does not carry out the driver equilibration
normally performed around LU/ILU. The low-level option alone does not
apply any `R` or `C` factors. Measurements made by changing only this
option are **not evidence** for or against equilibration. The original
draft of this experiment used that ineffective toggle and was corrected.

**The current script explicitly constructs** the max-norm scaling:
`R_ii = 1 / max_j |A_ij|` (for nonzero rows), followed by
`C_jj = 1 / max_i |(RA)_ij|` (for nonzero columns).
If `--row-max-scale` is already specified, `R` is common to all four
cases and the equilibration-on variant adds only `C`. If it is not
specified, the equilibration-on variant computes both factors.
This is a **manual max-norm analogue**, not an exact reproduction of
SuperLU's conditional `dgsequ`/`dlaqgs` driver or optional matching
and pivoting strategies.

To compare preconditioners fairly, GMRES's linear operator and right-hand
side **remain the same across all four cases**. The preconditioner applies
`M^{-1}=C*ILU(R*A*C)^{-1}*R` (with `R` already included in
the solver operator when `--row-max-scale` is selected).
In particular, the right-hand side is **not additionally transformed by C**
when running GMRES. The returned physical solution is unpermuted and its
true original `||b-A*x||_2/||b||_2` residual is checked.

## Inputs and execution

Use actual assembled systems exported with
[`lis_snapshot_export_mm()`](lis_snapshot_export.md), e.g.
`snapshots/step03_iter07_A.mtx` and
`snapshots/step03_iter07_b.mtx`. Alternately provide one explicit
`--rhs PATH` per matrix, in the same order.

No raw CCX snapshots are committed to the public LIS repository.

Requirements: Python 3, NumPy and SciPy with `spilu` and modern
`gmres(rtol=..., callback_type=...)`.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/compare_superlu_ilut_equil.py \
  snapshots/*_A.mtx \
  --rcm --row-max-scale \
  --drop-tol 1e-4 --fill-factor 10 \
  --rtol 1e-8 --restart 40 --maxiter 40 \
  --csv superlu_equil_ablation.csv \
  --json superlu_equil_ablation.json
```

`--rcm` applies an identical symmetric permutation for every variant.
`--row-max-scale` applies identical external max-norm row scaling in
each variant, isolating the incremental value of `C`.
`spilu` uses the same `drop_tol`, `fill_factor`, row-pivoting policy
(`diag_pivot_thresh`, default unless specified), drop rules and other
settings in all variants. The inactive `spilu` Equil option is set to
False in every case to prevent misleading interpretation.

`--manufactured-rhs` explicitly substitutes
`b = A @ random_vector(seed=42)` if actual RHS data are unavailable.
Such a run is a smoke test, not evidence on physical contact systems.

## Measured contact snapshots (2026-11-13)

An uploaded ten-step series (`drop001`) was evaluated locally on the
original A/b Matrix Market pairs (not committed to GitHub). Measurements:
SciPy **1.17.0**, NumPy **2.3.5**, 64-bit real, one BLAS/OpenMP thread,
`drop_tol=1e-4`, `fill_factor=10`, default SuperLU pivot threshold,
symmetric RCM, maximum-entry row normalization, restarted GMRES(40),
`rtol=1e-8`, `atol=0`. The five odd steps share an identical A
(SHA-256 verified), but their RHS differ by proportional scaling.
There are **six distinct A matrices** and **ten RHS**. Factors were
reused across identical A in the local test; each RHS was solved.

| Newton step | NATURAL R-only → R+C | COLAMD R-only → R+C | Columns with nontrivial C |
| ---: | ---: | ---: | ---: |
| 1 | 6 → 6 | 6 → 6 | 2 |
| 2 | 19 → 19 | 13 → 13 | 409 |
| 3 | 6 → 6 | 6 → 6 | 2 |
| 4 | 16 → 16 | 17 → 18 | 7 |
| 5 | 6 → 6 | 6 → 6 | 2 |
| 6 | 18 → 18 | 19 → 19 | 2 |
| 7 | 6 → 6 | 6 → 6 | 2 |
| 8 | 18 → 18 | 22 → 22 | 2 |
| 9 | 6 → 6 | 6 → 6 | 2 |
| 10 | 18 → 18 | 23 → 23 | 2 |

All **40/40** preconditioned GMRES runs converged under the
physical-residual acceptance threshold of `10*rtol`. Max observed
physical relative residual was approximately `1.25e-8`.
The additional max-norm column scale did not reduce iterations in
any case; the only change was one extra iteration with COLAMD at
step 4. The column correction is nearly trivial: 2/60,265 columns
at most steps, 7/60,265 in step 4, and 409/60,265 in step 2.

With row-only scaling, COLAMD reduced `nnz(L)+nnz(U)` by
roughly **5.7%–13.0%** versus NATURAL across the six distinct matrices.
However, COLAMD improved iteration counts **only in step 2**
(19 → 13), worsened them in steps 4, 6, 8, 10, and had no effect
in five repeated odd steps. Column ordering is therefore not a
universal convergence improvement here.

**Decision:** these 10 contact systems do **not** justify
introducing a new LIS default-off `-ilut_equil 1` path merely
to obtain additional max-norm column scaling. They also do not
isolate the benefits of SuperLU row pivoting, which was enabled
under the same policy throughout. LIS ILUT and SuperLU ILU are
different algorithms: a transfer must be separately validated
before performance claims can be made.

## Smoke tests

```sh
python3 -m unittest discover -s tools \
  -p 'test_compare_superlu_ilut_equil.py' -v

python3 tools/compare_superlu_ilut_equil.py \
  test/testmat.mtx --manufactured-rhs --rcm --row-max-scale \
  --csv /tmp/superlu_equil_smoke.csv
```

The small synthetic test explicitly includes a row with an off-diagonal
larger than its diagonal, making the column scaling nontrivial.
Its assertions prevent accidentally reverting to the no-op
`spilu(options={"Equil": ...})` comparison.
