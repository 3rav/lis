# GCRO-DR: opt-in harmonic Ritz update cadence

GCRO-DR updates its harmonic Ritz recycle subspace at each eligible Arnoldi
restart. On matrices that need many restarted Krylov cycles, this can be a
substantial part of the total solve CPU time.

## Usage

The default is **unchanged**: update at every eligible restart.

```text
-i gcrodr -p jacobi -restart 12 -gcrodr_keep 4
```

To update less often, use the **per-solver** option:

```text
-i gcrodr -p jacobi -restart 12 -gcrodr_keep 4 -gcrodr_ritz_stride 4
```

`-gcrodr_ritz_stride 4` performs the first eligible Ritz BUILD/UPDATE,
then retains the current recycle subspace for the next three eligible
cycles. The update count resets at the start of each new solve; iterations
and restart counts may change relative to the default. Only positive
integer values in `1..INT_MAX` are accepted.

Skipped updates preserve the existing `U / ZU / C` relation and do not
discard, reinitialize or replace previously generated recycle vectors.
Cross-solve assessment and preconditioner refresh behavior remain active.
This is different from disabling recycling or from the flexible
FGCRO-DR `-fgcrodr_intra 1` mode.

## Storage and MPI safety

The cadence lives in a solver-owned private configuration allocation that
already exists in LIS for solver-specific options. No members are added
to the public `LIS_SOLVER_STRUCT`; no existing option indices, binary
interfaces or default settings change.

The generic solver option-profile comparison/copy path includes this
setting, so replacing the active profile with a different cadence
invalidates the old solve-dependent sequence state as appropriate.

For MPI builds, the solver checks that participating ranks agree on the
cadence, recycle target dimension, restart dimension and maximum iteration
limit **before** entering Krylov collectives. Differing values are
rejected with an error rather than risking a deadlock.

## Verification

`test/testgcrodrstride.c` exercises:

- the default cadence against explicit `stride=1` with identical
  iteration counts, Ritz updates, and cross-solve recycle counts;
- `stride=4` with fewer harmonic update operations and true residual
  below `1e-8` for multiple right-hand sides;
- invalid option values and configuration comparison/replacement;
- a deliberately inconsistent per-rank cadence under two-process MPI.

The existing GCRO-DR and FGCRO-DR regressions are also run in focused
Release, AddressSanitizer/UndefinedBehaviorSanitizer, and MPI CI.

## Limitations

This is an opt-in performance control; it does **not** guarantee faster
solves or fewer Krylov iterations. Synthetic measurements in the
diagnostic PR #106 showed improvements with Jacobi, while a strong ILUTP
preconditioner often left little amortizable Ritz work. Before deploying
in contact/nonlinear finite-element workloads, benchmark representative
matrix sequences on a stable host and verify physical true residuals and
end-to-end wall time. It is not intended to replace FGMRES universally.
