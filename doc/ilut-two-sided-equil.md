# Optional two-sided ILUT equilibration

This branch adds opt-in two-sided local equilibration to serial real-double CSR ILUT, without SuperLU.

## Usage

```text
-p ilut -ilut_equil 1 -ilut_rcm 1 -ilut_drop 0.005 -ilut_rate 12
```

- `-ilut_equil 0` (default): retain existing ILUT behavior.
- `-ilut_equil 1`: row maximum-absolute-value scaling followed by column maximum-absolute-value scaling of the row-scaled matrix.
- `-ilut_row_scale 1`: existing row-only mode; when Equil is on, row scaling happens once, regardless of the row-scale option.
- `-ilut_rcm 1`: optional existing symmetric reordering, applied before scaling.

Let B = P A P^T. The factorized matrix is M = R B C, where R and C are the
diagonal row and column weights obtained by the two passes.

Forward apply: x = P^T C (LU)^(-1) R P b.
Adjoint apply: x = P^T R^H (LU)^(-H) C^H P b.

The Krylov operator is not scaled; the transformation is used only by ILUT.
No changes to legacy ILUT, ILU(k), MFAMG or other solvers when disabled.

## Constraints

- Serial real-double CSR; complex, distributed and matrix-free inputs are rejected.
- Frozen-pattern numeric refresh cannot update transformed coefficients and requires a full factor rebuild.
- The column pass costs O(nnz) time and O(n) additional memory.
- No pivoting, new reordering, or external SuperLU dependency.
- No guarantee of improved convergence; previous SciPy/SuperLU ablations showed little benefit from column scaling on the first ten FEA snapshots.

## Validation

```sh
autoreconf -fiv
./configure --enable-test --disable-fortran --disable-f90 --disable-saamg --disable-mpi
make -j2
./test/testilutequil
make -j2 check
```

`testilutequil` exercises all combinations of RCM, row-only and two-sided
scaling on a nonsymmetric ill-scaled matrix. It verifies scaled column maxima,
full-fill forward / transpose solves, option replacement, and invalid values.

On contact benchmarks compare true residual ||Ax-b||/||b||, setup time, total
wall time, peak memory and accepted Newton steps. Large-scale validation
on the motivating 60,265-DOF snapshot cases is still pending.
