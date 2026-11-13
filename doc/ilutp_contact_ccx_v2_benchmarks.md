# Native LIS ILUTP validation on CCX contact Newton matrices

**Status: mixed outcome for the *limited-fill* ILUTP mode; retain this negative-control benchmark.** These results are from an actual serial LIS executable built by GitHub Actions, not SciPy or SuperLU. A later independent complete-fill mode has additional validation; see the note at the end and [15-step native results](ilutp-native-contact-v2-15steps.md).

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

**Bounded conclusion:** limited-fill ILUTP with *pivot regularization
1e-9* and the stated drop/fill settings failed. The observed
`max(abs(D))=1e9` is consistent with a pivot being limited by the
1e-9 regularization floor; this observation is not proof of the cause.
The effects of removing diagonal perturbation versus removing dropping
have **not been isolated**. Do not infer from these runs alone that a
different factorization family is required.

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

**Keep upstream PR Draft and ILUTP opt-in.** The 187-spring
limited-fill success is useful, but the separately labelled 256-spring
case remains a negative control **for the measured settings**.
Before asserting a root cause, run the two-case unregularized-pivot
experiment below. End-to-end CCX behavior and a baseline comparison
against the application's production solver are still unmeasured.

These are direct linear-solver snapshots, not a nonlinear CCX increment acceptance test. No claim is made that the overall CCX nonlinear analysis now converges.


## Unresolved 256-spring attribution: next required experiment

**Not yet performed on the original 256-spring snapshot 03**. The
available archive named `lis_snapshots(1).7z` instead contains ten
`drop001` matrices (roughly 3.46–3.48M stored entries), not the original
`lis_snapshots_v2.7z` with snapshot 03 (3,639,839 entries). Do not
substitute a same-numbered or similarly sized matrix.

After building this PR, use the actual snapshot 03 *A and b* to compare,
holding RCM, row scaling, solver tolerance and other settings fixed:

```sh
A=/path/to/snapshot03_256_A.mtx
B=/path/to/snapshot03_256_b.mtx
OMP_NUM_THREADS=1 ./test/benchilutp "$A" "$B" \
  --maxiter 240 --restart 120 --rtol 1e-5 \
  --drop 0.005 --rate 12 --pivot-tol 0 \
  --rcm 1 --row-scale 1 --csv limited-pivot0.csv

OMP_NUM_THREADS=1 ./test/benchilutp "$A" "$B" \
  --maxiter 240 --restart 120 --rtol 1e-5 \
  --drop 0.005 --rate 12 --pivot-tol 0 \
  --rcm 1 --row-scale 1 --full 1 --only-full 1 \
  --full-tau 0.1 --full-limit 60000000 --csv full-pivot0.csv
```

The limited mode emits both `tau=0.5` and `tau=1.0`; the full
mode's threshold here is `0.1`, matching the existing full-fill
contact run. **The threshold is a remaining difference** between the
two; use `--full-tau 0.5` as an optional, matching-threshold control
before attributing success specifically to fill.
Evaluate original-system residuals, setup/solve CPU, factor NNZ,
and peak RSS; a successful solver status alone is insufficient.

If limited-fill with `--pivot-tol 0` succeeds, a full-fill rescue
cannot be attributed to fill alone. If only full fill succeeds, that
supports an important role for fill/dropping under these settings but
does not isolate the threshold or guarantee end-to-end CCX convergence.

## Separate 15-step complete-fill Newton series (private)

A separately uploaded private archive `lis_snapshots_v2(1).7z`
with prefix `pivot1e9_test2` contained **15 A/b pairs but only 12
distinct matrix A hashes** (steps 1, 4, 6 and 12 share A).
Its spring count is unknown; **this is not evidence that the separately
labelled 256-spring snapshot 03 has been solved**.

With native column-pivoted full fill, RCM=1, row-max scaling,
threshold=0.1 and a cap of 60 million stored L+U off-diagonal
entries, all 15 RHS were solved independently with one FGMRES
iteration per case. Original-system residuals were in
[4.631e-15, 2.837e-12]. Each run refactorized independently,
including repeated matrices; hashes of the source inputs were
checked before and after execution. This is **15 solves,
12 distinct matrices**, not 15 independent numerical problems.

| Step | Factor CPU s | L+U nnz | Peak RSS KiB | True relative residual |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 9.59 | 23533512 | 431836 | 4.631e-15 |
| 2 | 38.30 | 39107753 | 638232 | 2.028e-14 |
| 3 | 65.04 | 51108311 | 799024 | 1.664e-12 |
| 4 | 9.11 | 23533512 | 431868 | 4.631e-15 |
| 5 | 28.61 | 33338232 | 561388 | 3.816e-14 |
| 6 | 9.79 | 23533512 | 431704 | 4.631e-15 |
| 7 | 12.99 | 26710622 | 473212 | 4.943e-14 |
| 8 | 12.13 | 26021990 | 464356 | 5.722e-14 |
| 9 | 10.55 | 24611888 | 444312 | 2.837e-12 |
| 10 | 10.32 | 24756200 | 446168 | 5.145e-14 |
| 11 | 11.31 | 25737599 | 459092 | 5.952e-13 |
| 12 | 9.50 | 23533512 | 431676 | 4.631e-15 |
| 13 | 10.88 | 24986252 | 449248 | 8.641e-14 |
| 14 | 10.29 | 24786986 | 446804 | 4.747e-13 |
| 15 | 10.10 | 24440936 | 442352 | 1.513e-13 |

The maximum L+U fill was **51,108,311 entries** (step 3),
peak RSS **799,024 KiB**, and longest factor setup **65.04 CPU s**.
Step 3 needed a 150-second rather than the initial 45-second timeout.
The 60-million fill cap leaves approximately **17% headroom relative
to the 51.1-million observed fill**, not a guaranteed memory bound.

Earlier full-fill attempts with RCM=0 and tau=1 were **resource-limited**:
they exceeded the configured fill budget or timed out; they do not
establish numerical breakdown. The subsequent successful setup changed
RCM and pivot threshold together, so neither effect can be assigned
causally from those measurements.

Reproduction on *that* private series (once available) is through
`tools/run_native_ilutp_contact.py`, with `--full-limit 60000000
--full-tau 0.1 --rcm 1 --row-scale 1 --restart 40 --rtol 1e-6
--accept-residual 1e-5 --timeout 150`.

**Important limitations:** this is only an offline linear-system
test on one private FEA model, not a CCX nonlinear simulation,
not a test of the original 256-spring failure, and not a comparison
to the production CCX solver. A production sparse direct solver
with AMD/METIS may be faster; no measurement here establishes that.
The limit counts stored factor coefficients, not total RSS.
