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

## Reproducible positive cases (native LIS, Linux GCC, 2026-11-14)

Common parameters: Equil flag 0/1, ILUTP flag 0/1 with `tau=0.5`,
`-ilut_rcm 0 -ilut_row_scale 0 -ilut_drop 0.005 -ilut_rate 3
-ilu_pivot_tol 0.001 -restart 70 -maxiter 180 -tol 1e-6`.
The factor size below is actual stored off-diagonal L+U, **not**
a memory budget estimate.

| Grid, unknowns | Mode | Iterations | L+U nnz | Status |
|---|---|---:|---:|---|
| 8×8, n=192 | ILUTP only | >180 | 2,646 | fails |
| 8×8, n=192 | Equil-only ILUT | >180 | 2,675 | fails |
| 8×8, n=192 | **ILUTP + Equil** | **10** | **2,883** | **converges** |
| 10×10, n=300 | ILUTP only | >180 | 4,169 | fails |
| 10×10, n=300 | Equil-only ILUT | >180 | 4,192 | fails |
| 10×10, n=300 | **ILUTP + Equil** | **11** | **4,576** | **converges** |

The combined preconditioner gives:
- n=192: true relative residual **3.909e-7**, relative solution
  error **3.365e-5**; +7.8% and +9.0% L+U fill compared with
  Equil-only and ILUTP-only.
- n=300: true relative residual **5.736e-7**, relative solution
  error **6.613e-5**; +9.2% and +9.8% L+U fill.

An independent NumPy/SciPy eigenvalue check of the *input matrices*
confirmed: n=192 has 128 positive and 64 negative eigenvalues with
64 exactly zero diagonal entries; n=300 has 200 positive and 100
negative eigenvalues with 100 exactly zero diagonal entries.
Both input matrices have full numerical rank.

These are robust **numerical synergy** examples where column pivoting
and row/column Equil together enable convergence at modest extra fill.
This does **not** show a speedup in total setup+solve time. It also
does not demonstrate the same result on a genuine Stokes FE export
or CalculiX contact matrix.

## Reproduce n=192

```sh
python3 test/gen_stokes_saddle.py --side 8 --contrast 30 \
  --units pressure_weak --ordering interleaved --gamma 1 \
  --delta 0 --matrix A.mtx --rhs b.mtx
./test/benchilutp A.mtx b.mtx --maxiter 180 --restart 70 \
  --rtol 1e-6 --drop 0.005 --rate 3 --pivot-tol 0.001 \
  --rcm 0 --row-scale 0 --manufactured-truth 1 --equil 0 \
  --csv controls.csv
./test/benchilutp A.mtx b.mtx --maxiter 180 --restart 70 \
  --rtol 1e-6 --drop 0.005 --rate 3 --pivot-tol 0.001 \
  --rcm 0 --row-scale 0 --manufactured-truth 1 --equil 1 \
  --csv scaled.csv
```

From `controls.csv` compare `variant=ilutp_0p5` (ILUTP-only).
From `scaled.csv` compare `variant=legacy_ilut` (Equil-only)
and `variant=ilutp_0p5` (combined).

This branch carries an *experimental-only* `--manufactured-truth 1`
extension to `benchilutp`. The base #230 provides `--equil 0|1`
and the physical residual but does not carry this extra
manufactured-solution audit flag.

For the complete deterministic sweep:

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
