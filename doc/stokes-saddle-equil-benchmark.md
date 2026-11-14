# Native LIS: a reproducible indefinite Stokes-type saddle-point example

**Result:** Native FGMRES + ILUTP with optional two-sided Equil succeeds on
saddle-point systems where **Equil-only** and **ILUTP-only** individually fail
to reach the same iteration cap. Success is checked using BOTH the
true original-system residual and error against a known exact solution.

## Matrix construction

Use `test/gen_stokes_saddle.py`; the source is a **manufactured
Stokes-type finite-difference model** inspired by incompressible flow.
It is NOT an actual IFISS discretization and is NOT an extracted
CalculiX contact matrix.

Pressure block first in the mathematical representation:

```text
H = [ -delta*I     gamma*B ]
    [ gamma*B^T       K  ]
A = D H D
```

- `K`: two uncoupled, SPD, five-point 2-D velocity diffusion stencils,
  diagonal 4.5 and edge coupling -1.0.
- `B`: one-sided pressure/velocity divergence-like stencil. The x
  component has unit diagonal and lower-neighbor -1; hence `B` has
  **full row rank** and `H` is **nonsingular**, even when `delta=0`.
- `delta=0`: exactly zero pressure diagonal. The original `H` is
  symmetric indefinite with `2m²` positive and `m²` negative
  eigenvalues for an `m × m` grid.
- `D`: positive diagonal physical-unit transformation. For the
  selected `pressure_weak` examples pressure unknowns use `1/30`
  and velocity unknowns use `1`. Congruence preserves inertia.
- `ordering=interleaved`: pressure, x-velocity, y-velocity at each cell.
- RHS is always calculated from `b=A*x_true`, where
  `x_true[i]=cos(0.17*i)+0.25*sin(0.11*i)`.

The original A and b are written to Matrix Market files, with no
per-method modification. LIS computes the physical true residual
`||Ax-b||₂/||b||₂` and an independent **true relative solution error**
`||x-x_true||₂/||x_true||₂`. The acceptance thresholds are
`physical_relres<=1e-5` and `relative_solution_error<=1e-3` (native
FGMRES `-tol 1e-6`). The latter guards against ill-conditioning
masking solution error in the physical residual.

## Strong six-profile positive saddle examples (native LIS, 2026-11-14)

These selected cases compare six numerical methods at exactly the
same `-ilut_drop 0.005 -ilut_rate 3 -ilu_pivot_tol 0.001
-restart 70 -maxiter 180 -tol 1e-6 -ilut_rcm 0`, with
`-ilutp_threshold 0.5` for the pivoting variants.

At `--side 8 --contrast 30 --units mixed_checker
--ordering interleaved --gamma 0.3 --delta 0`,
**n=192 (64 pressure unknowns, 64 zero diagonal entries)**:

| Mode | Native status | Iterations | L+U offdiag nnz |
|---|---|---:|---:|
| Ordinary ILUT | not converged | >180 | 2,534 |
| Row-scaled ILUT | not converged | >180 | 2,534 |
| ILUTP without scaling | not converged | >180 | 2,618 |
| Row-scaled ILUTP | not converged | >180 | 2,586 |
| Two-sided Equil ILUT (no pivoting) | not converged | >180 | 2,734 |
| **ILUTP + two-sided Equil** | **converged** | **9** | **2,889** |

Combined: true physical residual **4.778e-7**, relative error
against the known exact solution **5.378e-4**, and **+5.7%**
actual factor nnz compared with the largest control.

At `--side 10 --contrast 30 --units mixed_checker
--ordering interleaved --gamma 1 --delta 0`,
**n=300 (100 pressure unknowns, 100 zero diagonal entries)**:

| Mode | Native status | Iterations | L+U offdiag nnz |
|---|---|---:|---:|
| Ordinary ILUT | not converged | >180 | 3,984 |
| Row-scaled ILUT | not converged | >180 | 3,984 |
| ILUTP without scaling | not converged | >180 | 4,244 |
| Row-scaled ILUTP | not converged | >180 | 4,245 |
| Two-sided Equil ILUT | not converged | >180 | 4,318 |
| **ILUTP + two-sided Equil** | **converged** | **11** | **4,575** |

Combined: true physical residual **3.143e-7**, relative
manufactured-solution error **1.362e-4**, **+6.0%** actual
factor nnz compared with the largest control.

The `2m²` positive and `m²` negative eigenvalues follow
analytically from K SPD, B full row rank and D positive.
The original matrices are symmetric **indefinite and nonsingular**,
even though pressure diagonal entries are exactly zero.

These are **selected manufactured examples**: evidence of a
numerical benefit from combining the methods, not proof of a
general overall speedup or success on a real Stokes/CCX FEM matrix.
The measurements use the same original A and b across all variants.

## Reproduce the n=192 case

```sh
python3 test/gen_stokes_saddle.py --side 8 --contrast 30 \\
  --units mixed_checker --ordering interleaved --gamma 0.3 \\
  --delta 0 --matrix A.mtx --rhs b.mtx

./test/benchilutp A.mtx b.mtx --maxiter 180 --restart 70 \\
  --rtol 1e-6 --drop 0.005 --rate 3 --pivot-tol 0.001 \\
  --rcm 0 --row-scale 0 --manufactured-truth 1 --equil 0 \\
  --csv no-scale.csv

./test/benchilutp A.mtx b.mtx --maxiter 180 --restart 70 \\
  --rtol 1e-6 --drop 0.005 --rate 3 --pivot-tol 0.001 \\
  --rcm 0 --row-scale 1 --manufactured-truth 1 --equil 0 \\
  --csv row-only.csv

./test/benchilutp A.mtx b.mtx --maxiter 180 --restart 70 \\
  --rtol 1e-6 --drop 0.005 --rate 3 --pivot-tol 0.001 \\
  --rcm 0 --row-scale 0 --manufactured-truth 1 --equil 1 \\
  --csv two-sided.csv
```

The `legacy_ilut` rows are ordinary or Equil-only ILUT;
`ilutp_0p5` rows provide the pivoting controls and combination.

This experimental branch carries an **optional**
`--manufactured-truth 1` flag which is not part of upstream #230.
The generic `benchilutp` in #230 supports `--equil 0|1`
and true original-system residuals. The manufactured exact
solution test additionally verifies `||x-x_true||/||x_true||`
against an independent `1e-3` acceptance threshold.

Run the full deterministic parameter sweep with

```sh
python3 test/run_stokes_saddle_sweep.py --bench ./test/benchilutp
```

## Full sweep and provenance

Native GitHub Actions (run includes CSV and Matrix Market artifacts):
https://github.com/3rav/lis/actions/runs/43585381603

- Small n=48/108: 256 configurations, **57** strict-synergy cases
  after both residual and solution-error acceptance checks.
- Larger n=192/300: 192 configurations, **70** strict-synergy cases
  after the same checks.
- Sweep parameters include variable contrasts, unit modes, orderings,
  velocity-pressure coupling, drop and fill budget. The two cited cases
  both use the same fixed options across all three methods.

The bare true-residual criterion alone produced 61 strict positives
on the small grid. The stronger solution-error criterion reduced this
to 57. This illustrates why a manufactured exact solution matters.

## Relation to published benchmarks

The operator follows the standard Stokes block saddle-point form
described in the IFISS Stokes testproblem documentation and in
general saddle-point preconditioning literature. Do not attribute
these numerical data to IFISS or to SuiteSparse.

References:
- IFISS Stokes tests: https://www.cs.umd.edu/~elman/ifiss3.6/ifiss_guide_3.6.pdf
- Example saddle point model and preconditioners:
  https://se.copernicus.org/articles/11/2031/2020/
- Public GHS_indef/stokes64s has 12,546 DOF but a two-dimensional
  numerical nullspace; it is not interchangeable with this full-rank
  manufactured test: https://sparse.tamu.edu/GHS_indef/stokes64s
