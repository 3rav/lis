# Native restarted symmetric Krylov eigensolver

This module adds a self-contained, matrix-free, thick-restarted subspace
Krylov eigensolver to LIS. All iterations, Rayleigh-Ritz projection,
dense symmetric Jacobi diagonalization, metric reorthogonalization,
thick restarts, selection, and residual validation are implemented in C
within LIS. **The implementation requires no external eigensolver**;
only the standard C library and libm are used.

## Operator contract

The public header is `lis_krylov.h`; the symbol is
`lis_krylov_symmetric` in liblis. All vectors have local length n
and elements of type double; callbacks return zero on success. Matrices
may use any representation or be assembled only through matvec callbacks.

```c
lis_krylov_problem p = {
    n, LIS_KRYLOV_SHIFT_INVERT_B, sigma, context,
    apply_A, apply_B, solve_shifted
};
lis_krylov_options o = { nev, 0, 300, 1e-10, "LM" };
lis_krylov_stats stats;
int rc = lis_krylov_symmetric(&p, &o, eigenvalues,
                              eigenvectors_column_major, residuals, &stats);
```

The spectrum is that of the symmetric generalized problem
**A x = lambda B x** (B is the identity for standard mode).
Operators and the inverse action are supplied by the caller,
independent of storage format. The caller must prepare a reusable
factorization or linear-solver setup for the selected shift.

### Spectral transformations

| transform | inner product | transformed operator | eigenvalue map |
|---|---|---|---|
| `LIS_KRYLOV_STANDARD` | Euclidean | A | lambda = mu |
| `LIS_KRYLOV_SHIFT_INVERT_B` | B (SPD) | (A - sigma B)^(-1) B | lambda = sigma + 1/mu |
| `LIS_KRYLOV_SHIFT_INVERT_A` | A (SPD) | (A - sigma B)^(-1) A | lambda = sigma mu/(mu - 1) |

In SHIFT_INVERT_A, B may be indefinite and sigma must be nonzero.
No factorization is performed inside the eigensolver; `solve_shifted`
receives a right-hand side and applies the inverse of A - sigma B
using application-owned state.

The option `which` selects the **transformed** Ritz spectrum:
LM (largest magnitude; default), SM (smallest magnitude), LA
(largest algebraic), or SA (smallest algebraic).
Choose ncv > nev; the default is min(n,3*nev+8).
maxiter=0 defaults to 300 expansion cycles; tol=0 defaults to 1e-10.

## Convergence and diagnostics

The solver uses a B-orthogonal (or A-orthogonal) subspace, double-pass
modified Gram-Schmidt, exact small projected symmetric diagonalization,
and thick restarting with retained Ritz vectors and cached operator images.
All returned eigenpairs are certified against the **original** equation:

`||A*x - lambda*B*x|| / (||A*x|| + |lambda|*||B*x||)`.

Only the first stats.nconv entries are valid after partial/nonconvergent
iterations. stats also records iterations, restarts, operator applications,
metric products, linear solves and largest certified residual. There is no
static process-global solver state, and callback errors are propagated.

## Build / test

With LIS's development checkout:

```sh
autoreconf -fiv
mkdir build && cd build
../configure
make -j2
make -C test test_krylov_modes test_krylov_sparse
./test/test_krylov_modes
./test/test_krylov_sparse
```

For standalone precision and sanitizers:

```sh
cc -O2 -std=c99 -Wall -Wextra -Werror -fsanitize=address,undefined \
  -Iinclude src/esolver/lis_krylov.c test/test_krylov_modes.c -lm -o /tmp/krylov-modes
/tmp/krylov-modes
```

The test suite exercises standard symmetric and generalized cases,
positive and negative generalized eigenvalues, inverse transforms,
indefinite second operators, metric orthogonality, invalid options,
nonconvergence, callback errors, and a 160-degree-of-freedom
non-diagonal spring-chain reference with a single reused factorization.

## Scope and limitations

The initial native implementation handles serial/local,
real-double symmetric eigenproblems and 32-bit dimensions.
It has no MPI distribution or complex/quad precision support.
The method is a **thick-restarted Rayleigh-Ritz/Krylov subspace**
method, not a bit-for-bit implementation of an implicit QR-restarted
Lanczos routine. Performance depends on the supplied shifted solve
and spectral gap; large-model benchmarks and additional
eigensolver-interface integration are future validation work.
