# Experimental ILUTP for difficult contact systems

This branch provides an opt-in, **serial real-double CSR** variant of ILUT with
**dynamic column pivoting**. It is a separate numerical implementation; normal
`-p ilut` remains unchanged unless `-ilutp 1` is specified.

## Example

```text
-i fgmres -restart 120 -p ilut -ilutp 1 -ilutp_threshold 0.5 -ilut_rcm 1 -ilut_row_scale 1 -ilut_drop 0.005 -ilut_rate 12 -ilu_pivot_tol 1e-9 -omp_num_threads 1 -maxiter 1000 -scale none
```

- `-ilutp 0|1`: enable dynamic column pivoting (default 0).
- `-ilutp_threshold <number>`: pivot when `abs(diagonal) < threshold * max(abs(uneliminated row))`.
  Range [0,1], default 0.5 (pivot only when the diagonal is less than half the largest remaining entry).
- `-ilu_pivot_tol`: existing *diagonal regularization*, applied **after**
  column selection. It is not itself pivoting.
- `-ilut_rcm` and `-ilut_row_scale`: optional existing preprocessors, applied
  before ILUTP.

The factorization approximates `D_row P A P^T Q = L U`, where
`Q` is the dynamically selected column permutation. Both forward and
Hermitian-transpose preconditioner application account for `Q`.
No row permutations, MPI, complex scalars, FP16, or parallel factorization
are provided in this initial version. OpenMP builds require one active solver
thread (`-omp_num_threads 1`).

The factorization follows the existing ILUT rule of dropping new fill contributions below the row-average-based threshold before insertion, and uses the existing per-row L/U fill budget.  It does not apply an extra post-accumulation or final-factor drop.  A min-heap handles dynamically discovered lower entries in O(log n) per insertion/removal rather than quadratic selection scans.  With `-ilutp_threshold 0`, ILUTP delegates numerical factorization to the existing ILUT implementation and retains an identity column permutation, providing an exact no-pivot baseline. New upper column indices are kept in original factor-space numbering
during elimination, then converted to final positions once all pivots have
been selected. ILUTP intentionally disables the frozen-pattern update path:
reusing a pattern from a different column permutation is not valid.

## Regression test

Run after regenerating Autotools files:

```sh
autoreconf -fiv
./configure --enable-test --disable-fortran --disable-f90 --disable-saamg --disable-mpi
make -j2
./test/testilutp
make -j2 check
```

`testilutp` checks a nonsymmetric 4x4 matrix requiring column pivoting and
verifies the full-fill preconditioner against exact forward and transpose
solutions, with/without RCM and row scaling.  It also compares L/U/D
factors with legacy ILUT for `tau=0` (with/without RCM/row scaling),
runs an end-to-end limited-fill FGMRES solve, checks preconditioner
sequence reuse, pivot-threshold invalidation, and frozen-pattern rebuild\nfallback, and validates option parsing. OpenMP builds additionally verify\nthat increasing the thread count after setup is explicitly rejected.
Complex LIS builds skip this real-only regression.

## Development status

**Experimental and not yet benchmarked on the 60,265-DOF CCX contact matrices.**
Column pivoting is not equivalent to SuperLU's partial *row* pivoting, and
cannot resolve every row-rank deficiency or guarantee improved convergence.
Check actual residuals and compare with the existing serial ILUT baseline.

## Native contact snapshot validation

Real 60,265-DOF CCX contact benchmarks show a strong improvement at 187 springs, but an early BREAKDOWN and larger fill at 256 springs. See [native ILUTP versus ILUT contact benchmark](ilutp_contact_ccx_v2_benchmarks.md) for full physical residuals, factor fill, timing and reproduction. Keep this option experimental.
