# Native LIS complete-fill ILUTP (experimental; no SuperLU)

This experiment is implemented entirely in LIS. It does not compile, load or
link to SuperLU and adds no new library dependency. The branch contains no
SuperLU plugin. The normal incomplete ILUTP algorithm remains the default.

## Usage

```text
-i fgmres -p ilut -ilutp 1 -ilutp_full 1 -ilutp_full_limit 20000000 \
-ilut_rcm 1 -ilut_row_scale 1 -scale none -omp_num_threads 1
```

- `-ilutp_full 0` (default) preserves existing behavior.
- `-ilutp_full 1` disables dropping and diagonal perturbation, retains
  all elimination fill and uses strongest-available column pivoting (`tau=1`).
  Only the preconditioner is changed; Krylov solves the **original** A x = b.
- `-ilutp_full_limit N` controls the maximum total number of stored
  off-diagonal L+U entries. A value of zero (default) sets an automatic
  budget `max(20*nnz(A), 4*n)`. When the budget is exceeded, setup returns
  `LIS_OUT_OF_MEMORY`; the algorithm will not silently drop fill.
  This is a bound on factor nonzeros, **not** a bound on total process RSS.
- This mode requires `-ilutp 1`, serial real-double CSR and a single
  OpenMP solver thread. Existing sequence-preconditioner reuse is
  invalidated when these options change.

This is **column-pivoted** native sparse LU (A Q = L U, up to roundoff),
not SuperLU's row-pivoted or supernodal LU. It has no COLAMD ordering.
Full fill may require large memory and runtime, or fail on difficult
matrices; it is not a guaranteed substitute for a production direct solver.

## Independent numerical regression

Locally built native LIS, nonsymmetric sparse test matrix (32 x 32);
the original A and b were the same for both preconditioners:

| Variant | true relative residual |
| --- | ---: |
| ILUTP `drop=0.5, rate=0.2` | 8.54298012e-2 |
| Native full-fill ILUTP | 2.44508862e-16 |

The full factors had 290 L+U nonzeros. Ordinary and transpose solves
matched the manufactured solution, maximum absolute error 6.66e-16.
A deliberately tiny factor-nonzero budget failed cleanly; legacy ILUTP
tests including tau=0 equivalence remained green.

Run `./test/testilutp` after building LIS. The experimental GitHub Actions
workflow also checks that neither `liblis.a` nor the test executable has
unresolved/dynamic SuperLU dependencies.

**Contact matrices are not in this repository.** The 60,265-DOF 256-spring
and 187-spring CCX systems have NOT been tested with this mode. Before
claiming a fix, measure physical residual `||Ax-b||/||b||`, peak RSS,
fill size and wall time on the *original* CCX snapshots.
