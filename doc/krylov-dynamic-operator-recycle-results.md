# Krylov recycling on numerically changing sparse operators

**Date:** 2026-11-20. **Status:** Diagnostic, not an upstream optimization.

## Purpose and method

This benchmark compares FGMRES, GCRO-DR and FGCRO-DR on a single
nonsymmetric CSR matrix whose **values are changed between six successive
solves** while the matrix object and the solver object persist. The
experiment explicitly notifies LIS of every values-only change. Selected
diagonal entries receive and release extra stiffness to mimic active-set
transitions, while the sparsity pattern stays fixed.

These are synthetic operators, **not physical contact matrices exported
from CalculiX or another FEA package**.

For each solver, four full sequences are run, discarding the first
warmup and taking the median of three measured sequences. The test
alternates method order across trials. Inputs: side lengths 28 and 40
(784 and 1600 unknowns), Arnoldi restarts 12 and 24, k=4 or 8,
Jacobi and ILUTP (fill rate 2, pivot threshold 0.5), and one OpenMP
worker. Each entire sequence solve is timed in CPU seconds, including
preconditioner setup/rebuild and cross-solve recycle assessment.
RHS construction and matrix mutation are outside timing.

An independent physical relative residual is evaluated for every
solve and required to be at most 1e-8; no performance thresholds
are used as CI gates.

[Full Release/ASan+UBSan workflow and artifacts](https://github.com/3rav/lis/actions/runs/44689959494)

**Result: 2112/2112 converged Release solves.** Maximum independently
computed relative residual: 9.95e-9. Sanitizer smoke also passed.

## Median CPU for six changing-operator solves

n=1600, restart=12; milliseconds, lower is better.

| Method | Jacobi stride 1 | Jacobi stride 4 | ILUTP stride 1 | ILUTP stride 4 |
|---|---:|---:|---:|---:|
| FGMRES | 7.04 | 7.04 | 14.86 | 14.86 |
| GCRO-DR k=4 | 23.96 | 16.82 | 20.56 | 19.47 |
| GCRO-DR k=8 | 69.21 | 33.00 | 29.42 | 24.31 |
| FGCRO-DR k=4 | 27.27 | 17.02 | 19.82 | 19.97 |
| FGCRO-DR k=8 | 48.97 | 26.25 | 25.47 | 25.54 |

## Aggregate comparison

Sixteen positive-k configurations per preconditioner:
two sizes x two restarts x two recycling solvers x k=4/8.
The following compares stride=4 with stride=1.

| Statistic | Jacobi | ILUTP |
|---|---:|---:|
| Median relative speedup from stride 4 | **1.50x** | **1.00x** |
| Configurations improved by stride 4 | 16/16 | 8/16 |
| Median optimized time relative to FGMRES | 2.82x slower | 1.53x slower |

FGCRO-DR with k=0 delegates to FGMRES and approximately matches its
CPU time. In **all** 32 positive-k configurations, FGMRES remains
faster even after the reduction in Ritz update frequency.

## Interpretation

- The cadence option reduces a genuine cost when many restarts occur
  under a cheap preconditioner. The effect usually disappears when
  ILUTP produces convergence within only a few restart cycles.
- Recycling is not automatically worthwhile. Assessing and
  transporting the retained subspace and double MGS in GCRO-DR still
  contribute overhead.
- The default cadence should remain 1. A nondefault cadence may
  change iteration counts and should be selected with representative
  sequence measurements.
- The fixed sparsity pattern means this does not yet validate active
  contact DOF remapping or changing matrix dimensions.
- End-to-end FEA claims require **actual A/b snapshot sequences**,
  recorded preprocessing and consistent convergence tolerances.

This is a diagnostic branch that combines the opt-in strides from
fork PRs #107 and #108 for comparative measurement only. No
production default or public ABI is modified.

The CA-MFAMG integration workflow's generated-file comparison fails
also on an untouched upstream develop branch:
[baseline run](https://github.com/3rav/lis/actions/runs/44689050379).
Its unrelated failure must not be presented as passing native
CA-MFAMG validation.
