# FGCRO-DR: experimental within-solve harmonic Ritz recycling

**Status:** opt-in, technical review in progress ([PR #50](https://github.com/3rav/lis/pull/50)).
The feature is independent of MFAMG and does not require offline
near-nullspace vectors.

## Behavior

With the opt-in switch, FGCRO-DR activates the most recent harmonic
Ritz recycle subspace **after an Arnoldi restart within the same solve**.
It builds/updates the subspace from right-preconditioned search
directions, checks operator identity and structure/value generation,
then projects the physical residual using the saved relation A*Z=C
before proceeding with the current recycle dimension and new B stride.
True residual must still pass convergence criteria.

This is **not** an AMG mode generator and does not automatically
construct the offline six near-nullspace modes in MFAMG experiments.

## Opt-in configuration

Use normal LIS options, e.g.:

    -i fgcrodr -p jacobi -fgcrodr_k 1 -restart 4 -maxiter 300 -tol 1e-5

Process-wide environment variables (the feature switches are disabled
by default):

| Variable | Value | Effect |
| --- | --- | --- |
| LIS_FGCRODR_INTRA_SOLVE_RECYCLE | 1 | Enable recycling on later Arnoldi restarts in one solve |
| LIS_FGCRODR_INTRA_TRACE | 1 | Emit per-restart diagnostics when the above flag is also enabled |
| LIS_FGMRES_PROJECTED_RELATIVE | 1 | Use relative projected-residual semantics; true residual still certifies success |

The standalone test enables projected-relative normalization **both**
with intra-solve recycling ON and OFF. Therefore it does **not** prove
anything about the isolated benefit of the normalization switch.
Compare workloads only with equivalent normalization settings.

Recycled rank must remain smaller than restart size. Without the new
intra-solve switch, historical cross-solve FGCRO-DR remains the default.
For recycle dimension zero, the existing FGMRES path is used. Nonzero
recycle dimension is not implemented for the quad-precision FGCRO-DR
core.

## Independent public regressions

A new standalone test, test/testfgcrodrintra.c, runs automatically in
make check. It checks:

1. Restarted Arnoldi, harmonic Ritz build/update and actual reuse,
   using a deterministic 120-row nonsymmetric CSR system.
2. OFF and ON paths with identical operator, rhs, restart and tolerance.
3. True relative residual < 1e-5; no rejected directions; bounded
   orthogonality and A*Z-C errors below 1e-8.
4. A synthetic sequence A1 -> A2 with the *same* persistent solver but
   a distinct matrix and shifted diagonal, including the new operator
   generation in the saved recycling state.
5. The above synthetic operator sequence under two MPI ranks, with
   process-local row assembly.

PR #49 independently adds a MAXITER regression (test/testfgcrodrmaxiter.c)
for the active-recycle workspace cleanup path. Additional CI runs all
three targeted regression executables with AddressSanitizer,
LeakSanitizer and UndefinedBehaviorSanitizer enabled.

These checks cover specific source paths. They cannot guarantee that
every memory or numerical problem has been excluded.

## Non-goals and further validation

- No MFAMG v7/v8 source code, CalculiX integration or private snapshot
  file is included.
- The same-matrix value-update regression changes the diagonal,
  explicitly marks the value generation, reuses one solver and checks
  that the updated AZ=C relation and true residual are valid. The
  distinct-object A1 -> A2 case exercises a separate invalidation path.
- The manufactured A1 -> A2 sequence is **not** a CalculiX
  changing-contact Newton sequence.
- The feature is a technical-review candidate, not a demonstrated
  general performance improvement; real application benchmarking,
  robustness review and further integration testing remain necessary.
- Full contact replay and offline mode research stay in
  [draft PR #47](https://github.com/3rav/lis/pull/47).

## B-workspace allocation bound

When intra-solve recycling is enabled, workspace B reserves capacity
for the largest projected recycle rank. The requested rank can exceed
restart-1, but no cycle can retain more than restart-1 directions.
The allocator therefore caps prospective rows to restart-1, while
still preserving enough space for an already active cross-solve rank.

The standalone test also uses a deliberately oversized
`-fgcrodr_k 1000000` with restart=2. It checks that the solver
still converges, maintains actual rank 1, uses recycling, and passes
the true residual / A*Z=C checks, without reserving a million rows.
This is a resource-safety regression, not a performance benchmark.
