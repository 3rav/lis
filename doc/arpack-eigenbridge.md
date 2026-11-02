# Optional ARPACK-NG bridge for CalculiX frequency and buckling

This branch adds an optional **matrix-free, real symmetric, double-precision,
serial** generalized eigensolver to LIS 2.1.13. Existing LIS eigensolvers and
default builds are unchanged. It uses ARPACK-NG's C reverse-communication API.

## Build and dependencies

Requires ARPACK-NG with arpack.pc, arpack.h, dsaupd_c, dseupd_c,
a C compiler, pkg-config and standard LIS Autotools dependencies.

~~~sh
autoreconf -fiv
mkdir build-arpack && cd build-arpack
../configure --with-arpack
make -j2
make -C test test_arpack_modes test_arpack_sparse
./test/test_arpack_modes
./test/test_arpack_sparse
~~~

The default ./configure leaves the bridge disabled. The public header is
include/lis_arpack.h. When linking to the static LIS library from downstream
code, include the ARPACK linker flags from pkg-config --libs arpack.

## API and callback contract

~~~c
#include "lis_arpack.h"
/* Callback signature:
 * int apply(void *ctx, int n, const double *x, double *y);
 */
lis_arpack_problem problem = {
    neq, 3, 0.0, ctx, apply_K, apply_M, solve_K_minus_sigma_M
};
lis_arpack_options options = { nev, 0, 300, 1e-10, "LM" };
lis_arpack_stats stats;
int rc = lis_arpack_symmetric(&problem, &options, eigenvalues,
                              eigenvectors_column_major, residuals, &stats);
~~~

All callbacks operate on contiguous local vectors of length n, return zero
on success, and receive the context pointer supplied by the caller.
The solve callback must compute the inverse action of K - sigma*M,
**using one factorization per selected shift**, prepared by the caller.
No matrix assembly, dense conversion, or internal refactorization occurs.
The callback may bridge to LIS_MATRIX_USER, CSR storage, or a CCX solver.

### Mode 3: FREQUENCY

- Generalized symmetric system: K u = lambda M u, where M is SPD.
- ARPACK B=M and OP=(K - sigma M)^-1 M.
- Frequencies f=sqrt(lambda)/(2*pi) for non-negative eigenvalues.
- Typically sigma=0 with which=LM to extract lowest modes efficiently.

### Mode 4: BUCKLE

- K u = lambda KG u, with K positive definite and KG indefinite.
- ARPACK B=K and OP=(K - sigma KG)^-1 K.
- The shift sigma must be **nonzero**, and the geometric stiffness matrix
  must be provided as the second callback.
- Preserve the existing CalculiX outer buckling shift-correction loop.
  A single shift may not capture the correct lowest positive load factor.

For either mode, ARPACK handles the implicitly restarted Lanczos iterations.
The default ncv is min(n,2*nev+8); with a chosen ncv it must satisfy
nev < ncv <= n. which can be LM, LA, SA, SM, or BE.
Check the number of converged pairs before using output arrays.

## CalculiX integration points

CalculiX's arpack.c (*FREQUENCY) and arpackbu.c (*BUCKLE) may replace their
reverse-communication eigenpair loop with lis_arpack_symmetric, forwarding
K/M/KG matvec operations and the previously factorized shifted solve.
Map computed eigenvalues/eigenvectors back to the usual CCX output pipeline.

This branch **does not modify or build a CalculiX checkout**. The CCX
linear-algebra adapter and build/link configuration must be integrated in
the target ccx source tree. The numerical tests below do not prove an
end-to-end speedup for ccx. Benchmark on identical meshes and settings,
comparing wall times, factorization counts, modal frequencies/buckling
load factors and true residual norms before claiming an improvement.

## Tests and diagnostics

- test_arpack_modes.c: analytic generalized diagonal frequency problem;
  buckling with positive and negative eigenvalues and indefinite KG;
  B-orthogonality; residuals; invalid argument / callback error cases.
- test_arpack_sparse.c: 160-DOF non-diagonal spring-chain stiffness
  matrix; five analytic low modes; only one initial factorization.
- .github/workflows/arpack-eigen.yml: CI build with ARPACK and a
  baseline build without ARPACK.

Relative residual per eigenpair is
||Kx-lambda Mx||/(||Kx|| + |lambda| ||Mx||).
The API exposes stats.nconv, stats.arpack_info, operation and solve counts.
Nonconvergence returns LIS_ARPACK_NO_CONVERGENCE; errors are propagated,
not silently accepted.

Limitations: only real double, local/serial, symmetric systems with
32-bit ARPACK indices. MPI and cyclic-symmetry eigenproblems are not
implemented by this bridge.
