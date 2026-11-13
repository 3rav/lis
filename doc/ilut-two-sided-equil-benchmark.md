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

## Extended native LIS sweep (2026-11-13)

Native GitHub Actions: https://github.com/3rav/lis/actions/runs/43484101569

Five column-unit contrasts (5, 10, 30, 100, 300), five ILUT drop
thresholds (0.0001, 0.001, 0.005, 0.02, 0.05) and ten fill-rate
settings (0.5 through 5) produce 250 parameter configurations.
Each runs three LIS variants, for **750 successful solves**.

Two-sided Equil used fewer GMRES iterations than row-only in 60/250
parameter configurations. Only **three** comparisons also had
<=10% more actual L+U nonzeros (two distinct numerical configurations;
rate=1 and 1.25 share the same integer fill budget here).

| Units | Drop | Rate | Mode | GMRES iterations | Actual L+U nonzeros | True relative residual |
|---:|---:|---:|---|---:|---:|---:|
| 100 | 0.05 | 1 | row | 10 | 533 | 4.7783e-7 |
| 100 | 0.05 | 1 | both | 9 | 544 | 6.3871e-7 |
| 5 | 0.0001 | 2 | row | 11 | 1027 | 3.5519e-7 |
| 5 | 0.0001 | 2 | both | 10 | 986 | 8.3851e-7 |
| 100 | 0.05 | 2 | row | 7 | 588 | 8.0587e-7 |
| 100 | 0.05 | 2 | both | 4 | 759 | 3.7627e-8 |
| 100 | 0.005 | 1 | row | 10 | 538 | 7.1469e-7 |
| 100 | 0.005 | 1 | both | 13 | 549 | 8.2475e-7 |

Interpretation: the 10->9 case shows a **small** improvement with
only +2.1% factor fill, while the stronger 7->4 case requires
**+29.1%** factor fill. There are also configurations where Equil
**regresses**, including 10->13 iterations for +2% fill. These
comparisons all meet the same 1e-6 original-system residual tolerance.
The experiment is a controlled *manufactured* PDE and not a public
reservoir or production FEA matrix. It establishes a possible
setup/iteration trade-off, not a general time or memory advantage.

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

