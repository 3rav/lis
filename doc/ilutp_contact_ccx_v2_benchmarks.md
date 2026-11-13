# Native LIS ILUTP validation on CCX contact Newton matrices

**Status: mixed outcome; keep ILUTP experimental and opt-in.** These results are from an actual serial LIS executable built by GitHub Actions, not SciPy or SuperLU.

## Inputs and methodology

- The private CCX 2.24 snapshot archive `lis_snapshots_v2.7z` contains real physical Newton matrices and RHS vectors. Only statistics are published here, **not the underlying snapshots**.
- Snapshot `03` (256 contact springs): 60,265 x 60,265 and 3,639,839 stored entries.
- Snapshot `05` (187 contact springs): 60,265 x 60,265 and 3,634,871 stored entries.
- Binary: `test/benchilutp`, compiled from this PR in [GitHub Actions run 43486091537](https://github.com/3rav/lis/actions/runs/43486091537).
- Fixed settings: serial CSR, FGMRES, one thread, restart 120, maxiter 240, tol 1e-5, no LIS global scaling, RCM, row-max scaling, ILUT drop 0.005, rate 12, pivot regularization 1e-9.
- Actual residual: `||A*x-b||_2 / ||b||_2` on the original physical A and b. Reported times are single-thread CPU seconds. Factor setup is measured with `lis_solver_sequence_prepare_precon`; Krylov time is measured with `lis_solve_kernel` using exactly the borrowed factors, **without automatic sequence fallback/rebuilding**. Older measurements using `lis_solve_sequence` could include an unreported additional factorization after breakdown and are superseded by this table.
- Moved columns count the nonidentity entries of the **final** column permutation (not the number of swaps).

## 187 contact springs — ILUTP helps

| Variant | LIS status | Iter | Physical relres | Setup CPU s | Solve CPU s | L nnz | U nnz | Moved cols |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Legacy ILUT | MAXITER | 241 | 0.910575622983685 | 4.8061 | 10.0260 | 5,697,332 | 5,498,058 | 0 |
| **ILUTP τ=0.5** | **SUCCESS** | **76** | **9.5305e-6** | **4.8797** | **2.8044** | 5,415,236 | 5,087,917 | 40 |
| **ILUTP τ=1.0** | **SUCCESS** | **76** | **9.4012e-6** | **4.9938** | **2.8920** | 5,419,048 | 5,107,171 | 212 |

With tau=0.5, the preconditioner changed only 40 final permutation positions, yet it reached the target residual with approximately the same factor setup cost and lower fill than legacy ILUT. LIS may count maxiter+1 at the iteration limit.

## 256 contact springs — ILUTP does not fix breakdown

| Variant | LIS status | Iter | Physical relres | Setup CPU s | Solve CPU s | L nnz | U nnz | Moved cols |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Legacy ILUT | BREAKDOWN | 1 | 1 | 15.1286 | 0.0783 | 7,287,173 | 10,132,181 | 0 |
| ILUTP τ=0.5 | BREAKDOWN | 2 | 1 | 34.3652 | 0.2537 | 13,273,971 | 10,355,890 | 20,547 |
| ILUTP τ=1.0 | BREAKDOWN | 2 | 1 | 38.0714 | 0.1320 | 13,410,461 | 10,584,738 | 32,897 |

The true residual remained approximately 1 in all three cases. ILUTP caused more fill and more than doubled factor setup cost, while still breaking down. This is a *negative result* and must be visible alongside the positive 187-spring result.

Both input matrices have nonzero diagonal entries. Minimum `abs(A[i,i]) / max_j abs(A[i,j])`: 0.0815 (256 springs), 0.5858 (187 springs). The pathological pivots therefore arise during elimination rather than from simply zero diagonal entries of the original matrices.


## What causes the 256-spring breakdown?

A separate *diagnostic-only* build of LIS (not part of this PR) instrumented
the first four real-double FGMRES Arnoldi steps and scanned all stored ILUT
factor entries. Its source is on the
[diagnostic branch](https://github.com/3rav/lis/tree/diagnostic/fgmres-contact-breakdown-20261113).
For the exact 256-spring input:

- **All stored L/U/D entries are finite**, for both tested ILUTP thresholds.
- At tau=0.5, max absolute inverse diagonal = **1e9**, max |L| = **683.2**,
  max |U| = **43.97**.
- At tau=1.0, max absolute inverse diagonal = **1e9**, max |L| = **2125.8**,
  max |U| = **39.68**.
- The triangular preconditioner solve already produces a vector with
  max absolute component about **9.9e229** for tau=0.5; the unscaled
  norm computation overflows to infinity. The resulting first
  Givens rotation is nonfinite. This is an *amplification during triangular
  solves*, not nonfinite factor storage.
- An experimental FGMRES mode that normalizes the finite preconditioned
  vector before Arnoldi removes the immediate overflow, but does
  **not** fix convergence: after 120 iterations, the original
  physical relative residual is **1.009** (tau=0.5) or **0.980**
  (tau=1.0). Consecutive normalized preconditioned vectors are
  nearly collinear. **This normalization is not proposed for merging.**
- Stronger diagonal pivot regularization at 1e-3 and 1e-2 still causes
  immediate BREAKDOWN. At 1e-3 it reduces max |D| from 1e9 to 1e3
  but does not prevent triangular amplification.
- Reducing ILUT fill rate to 2 lets the legacy ILUT run for 240
  iterations without an immediate BREAKDOWN, but it is not converged
  (true relative residual ~1.059). ILUTP still BREAKDOWNs.
  Thus lower fill merely trades breakdown for stagnation in this state.

**Engineering conclusion:** the hard case needs a *numerically different
preconditioner*, e.g. stable row pivoting, matching-based ordering / scaling,
or a contact-aware block preconditioner; merely changing the column pivot
threshold, fill budget, or scaling FGMRES output does not suffice.

## Reproduce

Build this PR with serial real-double LIS and test targets. Keep the original CCX Matrix Market A and b pair and run:

    OMP_NUM_THREADS=1 ./test/benchilutp \
      /path/to/pivot1e9_test2_00000000000000000005_newton_A.mtx \
      /path/to/pivot1e9_test2_00000000000000000005_newton_b.mtx \
      --maxiter 240 --restart 120 --rtol 1e-5 \
      --drop 0.005 --rate 12 --pivot-tol 1e-9 \
      --rcm 1 --row-scale 1 --csv ccx_187.csv

Change snapshot 05 to 03 for the 256-spring experiment. The CLI reports setup/solve CPU time, original true residual, L/U fill and the final permutation moved count. The CLI was separately smoke-tested against a deterministic 100x100 Matrix Market input in CI.

## Recommendation

**Keep upstream PR Draft and ILUTP opt-in**; the 187-spring success is strong evidence for value on some contact systems, but the harder 256-spring case is still a breakdown with significant extra fill. Before promotion: diagnose early FGMRES/preconditioner breakdown, measure pivot growth and memory, and benchmark against matching-based ordering (such as MC64) and SuperLU row pivoting.

These are direct linear-solver snapshots, not a nonlinear CCX increment acceptance test. No claim is made that the overall CCX nonlinear analysis now converges.
