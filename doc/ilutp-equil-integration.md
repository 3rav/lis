# Two-sided ILUT Equil stacked on upstream ILUTP PR #229

**Base:** [anishida/lis#229](https://github.com/anishida/lis/pull/229) at
`acc2e820f0149ad82465e3b2bc074c5a83324d94` (ready for review).
This branch **does not replace** the upstream PR and preserves its
native, bounded complete-fill (`-ilutp_full`) implementation, its
tests, USER/HYBRID/FP16 protections, and its benchmarks.

## Options

```text
-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0.5
-ilut_equil 1 -ilut_rcm 1 -ilut_drop 0.005
-ilut_rate 12 -ilu_pivot_tol 1e-9 -omp_num_threads 1 -scale none
```

- `-ilut_equil 0|1` (default 0) computes row maximum-absolute-value
  scaling followed by column maximum-absolute-value scaling **after** the
  row pass. The original `-ilut_row_scale 1` is a row-only baseline;
  enabling Equil implies the row pass exactly once.
- `-ilut_rcm 0|1` optionally symmetrically reorders `A` first.
- `-ilutp 0|1` and `-ilutp_threshold [0,1]` are inherited from #229.
- `-ilutp_full 1 -ilutp_full_limit N` remains available with Equil:
  its sparse full-fill pivoting factorizes the already equilibrated
  matrix; `N` still bounds stored L+U entries. Full-fill itself
  is a feature of #229, not introduced in this integration.
- All new settings default off; standard ILUT and ILUTP keep their
  original behavior unless explicitly requested.

## Factorization and application

With symmetric permutation P, row diagonal R and column diagonal C,
and ILUTP's dynamic **column** permutation Q:

```text
B = P A P^T
S = R B C
S Q ≈ L U

forward:    P^T C Q (L U)^(-1) R P
adjoint:    P^T R^H (L U)^(-H) Q^T C^H P
```

Column scales are indexed **before** Q: for factor-space position i,
use C[Q[i]]. Both forward and transpose paths perform the mapping.
The original FGMRES Krylov operator is not modified.

## Safety and scope

- Only serial, real-double CSR factorization; OpenMP thread count must be
  one in the ILUTP apply path.
- No new MPI, complex, FP16, or distributed factorization.
- Frozen-pattern numeric refresh is not safe after Equil, RCM or ILUTP
  transformations; require a regular rebuild.
- Changing Equil invalidates the cached sequence preconditioner.
  Solver option replacement preserves ILUTP full-fill and Equil flags.
- Column scaling uses O(n) extra memory and O(nnz) preprocessing.
- Dynamic **column** pivoting is not SuperLU partial row pivoting.

## Regressions

`test/testilutp.c` retains all upstream #229 native, full-fill,
USER/HYBRID, zero/tiny-tau and sequence tests. It also covers all
8 RCM x row-only x Equil combinations in exact forward and transpose
solves, both incomplete and bounded full-fill modes. It checks pivot
swap telemetry and Equil toggling in a sequence solve.

`test/testilutequil.c` independently checks two-sided scaling with
and without RCM, column maxima, forward/transpose solves, invalid
options and profile preservation.

`test/benchilutpequil.c` compares five profiles on one manufactured
nonsymmetric convection/diffusion stencil: ILUT, Equil, ILUTP, the
combination, and row-only control. It reports solver status, actual
L+U fill, pivot swaps, and the true residual on original A.

## Existing small-matrix results (prior integration, not this exact head)

The previous integration branch passed 12/12 CI jobs.
A synthetic 144-unknown convection-diffusion case with contrast 30,
rate=5 and drop=0.005 gave:
- ILUT: 3 iterations, 703 L+U entries, 0 swaps;
- Equil: 3 iterations, 1282 entries, 0 swaps;
- ILUTP: 11 iterations, 907 entries, 71 swaps;
- ILUTP+Equil: 3 iterations, 1282 entries, 0 swaps.

**These are not results on the new #229 head**. They show interaction
but no independent convergence benefit over Equil alone, and no
comparable-fill advantage. Larger public-sparse/contact benchmarks,
true residual and full wall-time comparisons remain necessary.

Keep any downstream upstream PR **Draft until #229 is merged and a
fresh integration rebase / CI comparison has been reviewed**.
