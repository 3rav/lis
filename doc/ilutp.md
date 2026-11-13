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
sequence reuse, pivot-threshold invalidation, and frozen-pattern rebuild
fallback, and validates option parsing. OpenMP builds additionally verify
that increasing the thread count after setup is explicitly rejected.
Complex LIS builds skip this real-only regression.

## Development status

**Experimental.** Both the original limited-fill ILUTP and the newer native complete-fill option have been benchmarked on real 60,265-DOF CCX contact snapshots. See the distinct datasets and configurations below.
Column pivoting is not equivalent to SuperLU's partial *row* pivoting, and
cannot resolve every row-rank deficiency or guarantee improved convergence.
Check actual residuals and compare with the existing serial ILUT baseline.

## Optional complete-fill native LU (no external solver dependency)

The **separate opt-in** `-ilutp_full 1` mode retains all fill and disables
coefficient dropping/diagonal perturbation inside the ILUTP factorization.
It uses LIS's own column permutations and the existing
`-ilutp_threshold` policy; it is not SuperLU, a supernodal LU, or row pivoting.
For example:

```text
-i fgmres -p ilut -ilutp 1 -ilutp_full 1 \
-ilutp_threshold 0.1 -ilut_rcm 1 -ilut_row_scale 1 \
-ilutp_full_limit 60000000 -omp_num_threads 1 -scale none
```

- `-ilutp_full 0` (default) leaves legacy and limited-fill ILUTP unchanged.
- `-ilutp_full_limit N` is a cap on **stored off-diagonal L+U nonzeros**,
  not total process resident memory; `0` uses an automatic budget.
- Reaching the cap aborts factorization with a reported error rather than
  silently dropping additional coefficients.
- RCM and row scaling act only inside the preconditioner; the Krylov operator
  and independently verified residual use the original `A,b`.
- The preconditioner is valid only for serial, real-double CSR and a single
  active OpenMP solver thread. MPI, complex and general CCX fallback are
  not implemented by this mode. It can use substantially more time and memory.

See [full-fill implementation and regression](ilutp-native-fullfill.md) and
[15-step CCX Newton-series measurements](ilutp-native-contact-v2-15steps.md).

## Native contact snapshot validation

Earlier **limited-fill ILUTP** results remain important: at 187 springs,
ILUTP converged in 76 iterations; at 256 springs, limited-fill ILUTP broke
down with true relative residual near 1. The data and negative controls
are preserved in [the original contact benchmark](ilutp_contact_ccx_v2_benchmarks.md).

A later, separately supplied private `pivot1e9_test2` archive contains
**15 original 60,265-DOF Newton A,b pairs**. With complete fill, RCM,
row-max scaling, pivot threshold 0.1 and a 60-million L+U-nonzero limit,
**15/15 linear systems** achieved independently checked original-system
relative residuals between 4.63e-15 and 2.84e-12 in one FGMRES iteration
each. The maximum observed fill was 51,108,311 L+U nonzeros and the
largest observed peak RSS was 799,024 KiB (~780 MiB). The most expensive
factorization required 65.04 CPU seconds. These are **per-linear-system**
results, not an end-to-end nonlinear CCX run.

**Do not equate this 15-step archive with the earlier labelled 187/256
spring cases:** its metadata do not reliably establish spring counts.
The archive itself remains private; only summary statistics are committed.
This PR does not alter CCX or establish Newton-step acceptance.

Both modes are experimental and default-off. No SuperLU header, library,
runtime loading, or linked symbol is required.
