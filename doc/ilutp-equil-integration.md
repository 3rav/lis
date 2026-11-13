# Experimental integration: ILUTP + two-sided Equil

This branch combines two **independently optional** developments from
[ILUTP PR #75](https://github.com/3rav/lis/pull/75) and
[two-sided Equil PR #76](https://github.com/3rav/lis/pull/76)
without changing either source PR. It is based on PR #75's feature branch.

## Options

```text
-i fgmres -restart 120 -p ilut -ilutp 1 -ilutp_threshold 0.5 -ilut_equil 1 -ilut_rcm 1 -ilut_drop 0.005 -ilut_rate 12 -ilu_pivot_tol 1e-9 -omp_num_threads 1 -scale none
```

- `-ilutp 0|1`: dynamic **column** pivoting; default 0.
- `-ilutp_threshold [0,1]`: dynamic column-pivot threshold; default 0.5.
- `-ilut_equil 0|1`: one row-max pass followed by one column-max pass; default 0.
- `-ilut_rcm 0|1`: optional symmetric ordering.
- `-ilut_row_scale 0|1`: existing row-only normalization. Equil implies
  the row pass, whether this flag is zero or one.

All four main profiles are available: ordinary ILUT, ILUT + Equil, ILUTP,
and ILUTP + Equil. Row-only ILUT remains a separate reference.

## Exact factor-space mapping

With symmetric RCM permutation P, diagonal row weights R and diagonal
column weights C obtained from the **unpivoted** preconditioner matrix,
and ILUTP's dynamic column permutation Q,

```text
B = P A P^T
S = R B C
S Q ≈ L U

forward  : P^T C Q (LU)^(-1) R P
transpose: P^T R^H (LU)^(-H) Q^T C^H P
```

This ordering matters: the *column scaling vector is indexed in original
factor-space order*, before Q. During forward application the outgoing
factor-space coordinate i refers to column Q[i], so the coefficient is
C[Q[i]]. During transpose application the incoming coordinate uses
C[Q[i]]. Failing to remap C would silently give incorrect preconditioned
forward/transpose operators. The Krylov operator itself is not rescaled.

## Safety and compatibility

- Default settings preserve the existing ILUT behavior.
- Only **serial, real-double CSR** is supported for the combined path.
- ILUTP factorization and application require one OpenMP solver thread.
- The existing frozen-pattern numerical-update path is explicitly rejected
  for reordered / Equil / ILUTP preconditioners; the sequence manager must
  rebuild rather than reuse invalid patterns.
- Changing ILUTP/Equil options invalidates cached preconditioners and
  option-profile comparison includes all settings.
- MPI, complex arithmetic, FP16 factorization, and multi-threaded ILUTP
  remain unsupported; existing non-ILUTP code paths are unchanged.

## Verification

`make check` runs `testilutp` and `testilutequil`.
`testilutp` now checks all 8 combinations of RCM, row-only scaling
and Equil with dynamic pivoting, forward and transpose exact solves on
a small nonsymmetric zero-diagonal system, tau=0 equivalence, and
sequence rebuild when Equil is toggled.

Optional manufactured-solution benchmark:

```sh
./test/benchilutpequil 5 5 0.005
./test/benchilutpequil 30 2 0.005
```

The output has five profiles (`ilut`, `equil`, `ilutp`,
`ilutp_equil`, `row_only`) and reports the status, FGMRES iteration
count, **true original-system residual**, actual L+U fill, dynamic
column-pivot swap count, solver wall/setup/solve times and effective
drop/rate settings. It uses a synthetic 12x12 nonsymmetric
convection-diffusion stencil with alternating unknown units, **not**
a matrix extracted from production CalculiX or PETSc.

Do not infer a general convergence or wall-time improvement from green
unit tests. Before considering an upstream merge, compare true residual,
total setup+solve wall time, LU fill, pivot counts, and peak memory on a
larger public matrix and on contact snapshots; include negative controls.

## Initial native LIS five-profile example

The full Linux serial regression suite passed 110/110 tests and a prior
smoke benchmark completed successfully at amplitudes 5, 10 and 30:
https://github.com/3rav/lis/actions/runs/43489051029

At amplitude=30, drop=0.005, rate=5, all profiles reached the same
1e-6 solver tolerance on the **original**, unscaled A:

| Profile | FGMRES iterations | L+U nnz | True relative residual |
|---|---:|---:|---:|
| ILUT | 3 | 703 | 3.8142e-7 |
| Equil only | 3 | 1282 | 1.9900e-8 |
| ILUTP only | 11 | 907 | 4.4162e-7 |
| ILUTP + Equil | 3 | 1282 | 1.9900e-8 |
| Row-only ILUT | 3 | 703 | 3.8142e-7 |

This shows compatibility, **not unique synergy**: Equil alone also
converged in 3 iterations, and both Equil cases needed more fill. This
small manufactured matrix is primarily a numerical regression example.

