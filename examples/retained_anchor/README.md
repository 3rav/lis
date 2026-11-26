# Retained-anchor FGMRES example

This source-tree example demonstrates a **lagged preconditioner** for a sequence
of related sparse systems. It is not an installed API, a production contact
solver, or a recommendation to use complete-fill ILUTP for all problems.

The driver borrows a matrix P and reuses its ILUTP factors to solve A x = b.
A projected residual ratio at iteration 20 determines whether to continue;
the physical residual is checked independently of that estimate. If the
attempt fails numerically, the driver tries one fresh factorization of A.
It commits x only when the physical RHS-relative residual is at most 1e-8
(absolute residual for a zero RHS). Hard API errors are not retried.

## Requirements and ownership

- Native real-double, unscaled assembled CSR; same global DOF ordering.
- One active thread, outside an OpenMP parallel region; no MPI, complex or QUAD.
- Finite matrix entries and RHS. No aliasing among b, x and matrix storage.
- The caller owns P, A, b and x. Keep the current anchor alive and immutable
  until destroy (or successful rolling promotion). Matrix generation checks
  detect published changes, not dangling pointers or unannounced writes.
- Use synchronously, without concurrent input mutation. Recreate after a
  hard solve error. This example uses internal interfaces and must be built
  against the same source revision, with the same configuration.

## Build and run the example regression

From the repository root, with a matching static build:

    ./configure --disable-mpi --enable-omp --disable-fortran --disable-f90 --disable-saamg --disable-test --disable-shared --enable-static
    make -j2
    cc -O2 -fopenmp -DHAVE_CONFIG_H -Iinclude -Iexamples/retained_anchor examples/retained_anchor/lis_anchor_pilot.c test/testanchorpolicy.c src/.libs/liblis.a -lm -o test-anchor-policy
    OMP_NUM_THREADS=1 ./test-anchor-policy

For an already initialized LIS application with suitable P, A, b and x:

    lis_anchor_pilot *ctx = NULL;
    lis_anchor_pilot_report report;
    LIS_INT rc = lis_anchor_pilot_create_single_pass(P, &ctx);
    if (!rc) {
        rc = lis_anchor_pilot_solve(ctx, A, b, x, &report);
        /* rc == LIS_SUCCESS and report.accepted indicate a committed result. */
    }
    lis_anchor_pilot_destroy(ctx);

## Policy choices

- `create`: retain the original factor, restart after the 20-step pilot.
- `create_restart_matched`: adjust the continuation tolerance to match the
  original physical RHS-relative target.
- `create_single_pass`: retain the active Krylov basis on a successful
  checkpoint; no forced restart at iteration 20.
- `create_low_memory`: release the original factor before fresh fallback;
  rebuild lazily on the next solve.
- `create_rolling`: retain a successful fresh fallback factor as the next
  anchor. After **every** return query `lis_anchor_pilot_current_anchor(ctx)`
  before releasing either matrix. Failed calls do not change anchor identity.

The fixed reference profile is FGMRES restart 60, pilot 20, continuation 80,
fresh budget 100, ILUTP drop 0.005, rate 8, pivot tolerance 1e-9, threshold
0.5, RCM, row scaling and complete fill with limit 60000000. All policies
use zero initial guesses. These settings are deliberately explicit and
are not general tuning recommendations. Complete fill may be expensive,
and retaining factors can increase peak memory; low-memory mode is not
a strict RSS bound.

For related lagged/reused preconditioner interfaces, see PETSc
[KSPSetReusePreconditioner](https://petsc.org/release/manualpages/KSP/KSPSetReusePreconditioner/)
and [SNESSetLagPreconditioner](https://petsc.org/release/manualpages/SNES/SNESSetLagPreconditioner/).
