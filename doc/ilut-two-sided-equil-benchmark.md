# Numerical comparison: ILUT row scaling vs two-sided Equil

This benchmark is a **manufactured convection-diffusion system**, not an
extracted production matrix and not the original PETSc example. Its
five-point stencil follows the same class of finite-difference systems
used by PETSc KSP tutorials:
https://petsc.org/release/manual/getting_started/

## How to reproduce

```sh
autoreconf -fiv
./configure --enable-test --disable-fortran --disable-f90 --disable-saamg --disable-mpi
make -j2
./test/benchequil 100 1 0.005
./test/benchequil 100 5 0.005
```

The source is `test/benchequil.c`. It constructs a 12x12 nonsymmetric
five-point advection/diffusion stencil, then multiplies each column by
alternating units `amplitude` and `1/amplitude`. The exact solution is
prescribed, so all three variants use the same original `A` and `b`.
We compute the **true residual on unscaled A** after LIS FGMRES, rather
than using the preconditioned residual. The solver options are identical
except the scaling mode:

- `none`: no explicit ILUT scaling
- `row`: `-ilut_row_scale 1 -ilut_equil 0`
- `both`: `-ilut_row_scale 1 -ilut_equil 1`

Common options: `-i fgmres -restart 50 -p ilut -ilu_pivot_tol 0.001`,
`-maxiter 300 -tol 1e-6 -scale none -omp_num_threads 1`.
CLI: `benchequil [unit_amplitude] [ilut_rate] [ilut_drop]`.
The program records the *effective* rate/drop parameters, fill budget,
actual number of nonzero L+U coefficients, iterations and true residual.

## Initial native LIS results (verified by GitHub Actions)

At `amplitude=100`, `drop=0.005`:

| Maximum fill rate | Mode | Iterations | L+U nonzeros | True relative residual |
|---:|---|---:|---:|---:|
| 1 | row | 10 | 538 | 7.1469e-7 |
| 1 | both | 13 | 549 | 8.2475e-7 |
| 5 | row | 5 | 648 | 3.4507e-7 |
| 5 | both | 3 | 1282 | 1.9957e-8 |

For the initial sweep of 7 unit contrasts x 10 fill budgets, all 210
native solves converged. Two-sided scaling reduced iteration count
in **9 of 70** direct comparisons, but **none** of those wins had
<=10% extra actual L+U fill. Thus, the raw decrease from 5 to 3
iterations at rate=5 does *not* establish an improvement at equal memory.

Native run + raw artifact:
https://github.com/3rav/lis/actions/runs/43483544224

## Scope and caveats

- The example shows an opt-in trade-off, not a general superiority of
  two-sided equilibration over row scaling.
- It uses LIS's own ILUT, not SuperLU's differently defined dropping.
- The Python approximate-ILUT model used during exploration is **not**
  interchangeable with these native results.
- Compare total setup and solve time on larger systems before declaring
  a runtime benefit. This small matrix mainly validates behavior.
- `-ilut_equil 1` remains disabled by default and the PR remains draft
  pending a convincing comparable-fill advantage on public matrices.

