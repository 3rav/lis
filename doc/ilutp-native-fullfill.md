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


## Synthetic 60,265-DOF stress (NOT contact data)

Use `test/gen_ilutp_fullfill_case.py` to generate nonsymmetric,
ill-scaled banded A (tiny diagonals every 127th row; offsets ±1, -11,
+17) with manufactured RHS b=A*x. No SuperLU is involved. The original
Ax=b operator is unchanged by the test.

```sh
python3 test/gen_ilutp_fullfill_case.py --n 60265 \
  --matrix /tmp/native-A.mtx --rhs /tmp/native-b.mtx
./test/benchilutp /tmp/native-A.mtx /tmp/native-b.mtx \
  --full 1 --full-limit 20000000 --maxiter 80 --restart 20 \
  --drop 0.05 --rate 1 --rcm 0 --row-scale 1 \
  --csv /tmp/native-60265.csv
```

Local, one-thread measurement, 301,295 input nonzeros:

| Mode | Iterations | L+U nonzeros | Original relative residual |
| --- | ---: | ---: | ---: |
| ILUT | 6 | 241030 | 1.20587e-6 |
| ILUTP tau=0.5 | 3 | 241030 | 3.72596e-7 |
| ILUTP tau=1 | 3 | 241032 | 3.72586e-7 |
| Native full-fill ILUTP | 1 | 1687037 | 4.01561e-14 |

The full mode permuted 951 column positions. A separate full-only run
used ~0.113 CPU seconds for factor setup and peak RSS 57,392 KiB on
this specific environment; these figures are not contact-model estimates.
A deliberately insufficient limit of 10,000 factor nonzeros returned
`setup_failed` rather than silently dropping coefficients.

`--full 1` appends a fourth benchmark mode. Without it, the previous
three modes remain unchanged. `--only-full 1 --full 1` runs only the full
mode. `--full-limit N` checks a nonnegative integer within LIS_INT range.

This regression demonstrates a 60,265-DOF sparse control, NOT the
original 256-spring Newton system. Original CCX A,b are unavailable in
this repository. Its fill, spectrum, memory and convergence may differ.
