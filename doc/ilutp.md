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

## Optional native bounded complete-fill LU

```text
-i fgmres -p ilut -ilutp 1 -ilutp_full 1 \
-ilutp_threshold 0.1 -ilut_rcm 1 -ilut_row_scale 1 \
-ilutp_full_limit 60000000 -omp_num_threads 1 -scale none
```

- `-ilutp_full 0` is default. The opt-in full mode retains all
  elimination fill, disables dropping and *also disables diagonal
  regularization* (`pivot_tol=0` internally). It does not isolate
  any one of these changes, and is **not** equivalent to SuperLU's
  row-pivoted or supernodal factorization.
- `-ilutp_full_limit N` bounds the total stored **off-diagonal**
  L+U nonzeros, not process RSS. The default `0` uses the automatic
  budget `max(20*nnz(A), 4*n)`. Exhaustion fails setup rather than
  silently dropping fill.
- The full mode uses the existing column selection threshold but
  skips magnitude sorting of factors, which is unnecessary without
  dropping. It supports serial real-double CSR only, with one active
  OpenMP thread. MPI, complex and half precision are unsupported.
- Both RCM and row scaling are **internal to the preconditioner**;
  the Krylov operator and physical residual check still use the
  original matrix and RHS.
- Frozen-pattern numeric refresh is not supported for changing
  pivot permutations. There is no automatic nonlinear CCX fallback
  or validated factor-cache reuse across Newton steps.

**Compatibility clarification:** while the default `-p ilut`
factorization remains the legacy implementation when `-ilutp 0`,
this PR also updates `lis_solver_set_option_replace` and
`lis_solver_set_option_replace_if_changed` so the existing
`ilut_rcm` and `ilut_row_scale` options (as well as new
ILUTP options) are copied and compared on option replacement.
That change **also affects legacy ILUT's option lifecycle**.
Do not describe all legacy ILUT option behavior as byte-for-byte
unchanged. This change is necessary to avoid stale settings in
reused solver objects, and requires independent regression review.

## Validation and limitations

`./test/testilutp` exercises forward and transpose solves,
pivot-threshold behavior, limited versus complete fill, fill-limit
failure, and algorithm compatibility checks. A 32×32 nonsymmetric
regression reduced original relative residual from roughly
`8.54e-2` with deliberately weak limited fill to `2.45e-16`
with complete fill, without changing the physical operator.
A separate *synthetic*, not CCX, 60,265-DOF test from
`test/gen_ilutp_fullfill_case.py` tests factor fill and physical
residual; neither establishes general nonlinear convergence.

The **only contact benchmark report** is
[contact outcomes, negative controls and the 15-step series](ilutp_contact_ccx_v2_benchmarks.md).
It includes the separately labelled 187- and 256-spring limited-fill
results, the 15-run/12-distinct-matrix complete-fill series,
the controlled experiment still needed to isolate pivot
regularization, and command lines for reproduction.

Existing results do **not** establish that complete fill rather
than disabled regularization was the cause of improved behavior.
They do not compare against the CCX production solver, establish
end-to-end nonlinear convergence or guarantee bounded memory on
arbitrary problems. The default-off implementation remains
experimental and is not proposed as a general replacement for a
production direct solver.
