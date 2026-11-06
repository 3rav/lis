# Contact operator snapshots: a reproducible comparison contract

`tools/contact_matrix_audit.py` is a **read-only**, standard-library Python
utility for checking matrix snapshots exported from CalculiX or another FE
code. It does not export matrices from CalculiX, run LIS or run competing
solvers. It makes the inputs and correctness checks for later comparisons
reproducible.

## Input files

All vectors and matrices use [Matrix Market](https://math.nist.gov/MatrixMarket/formats.html)
format. The matrix `A.mtx` must contain the **assembled linear system**
that LIS actually solves, **after elimination of SPC/MPC** and the insertion
of the current contact tangent. A symmetric stiffness matrix may be stored
in lower or upper triangle with the header:

```text
%%MatrixMarket matrix coordinate real symmetric
```

General nonsymmetric matrices should use `real general`. The auditor
expands symmetric coordinate matrices and sums duplicates. It supports
real/integer/pattern coordinate matrices and real/integer general arrays;
complex, Hermitian and skew-symmetric formats are intentionally rejected.

The optional right-hand side `b.mtx`, initial guess `x0.mtx`, actual
solution `x.mtx`, and candidate mode matrix `B.mtx` are Matrix Market
vectors/arrays. A mode matrix has `n` rows and one column per RBM mode.
Matrix Market array data are **column-major**.

The optional DOF map CSV must have the following header:

```csv
eq,node,dof,body
0,123,1,part-A
1,123,2,part-A
2,124,3,part-B
```

- `eq` is a **0-based active equation index**, contiguous across `0..n-1`.
- `node` is an integer FE node identifier; `dof` is displacement
  component 1, 2 or 3.
- `body` is a stable structural-part identifier **before contact**.
  Do not infer it from the graph of an operator already containing
  penalty springs.
- Export exactly one mapping row for every active equation. Constrained
  DOFs excluded by SPC/MPC do not appear. Different bodies may each have
  an incomplete set of active DOFs.

## Usage

```sh
python3 tools/contact_matrix_audit.py \
  --matrix A.mtx \
  --rhs b.mtx \
  --initial x0.mtx \
  --solution x.mtx \
  --modes B.mtx \
  --dof-map active_dofs.csv > audit.json
```

The JSON reports:

- exact file SHA-256 for `A` and `b`; size, nonzero count, diagonal
  range and number of zero diagonals
- relative Frobenius asymmetry `||A - A^T||_F / ||A||_F` and
  `||A||_infinity`
- true `||b-A*x0||_2/||b||_2` and `||b-A*x||_2/||b||_2`
  (when `b=0`, denominator is 1, so the reported measure is absolute)
- per-mode Rayleigh quotient `v^T A v / (v^T v)` and
  `||A*v||_2 / (||A||_infinity * ||v||_2)`
- counts of bodies, physical nodes, and nodes with incomplete active DOFs.

The program fails with exit status **2** on malformed inputs, including
NaN/Infinity coefficients, out-of-range matrix indices, missing entries,
non-square stiffness matrices and duplicate or incomplete DOF maps.

## Suggested capture set for nonlinear contact

Record the **identical** operator and vectors at a fixed Newton iterate,
with active contact counts `0, 1, 2, 5, 12, 53, 116` when available.
Also keep separate metadata: mesh and material version, contact law,
Newton iteration, active spring count, convergence tolerance and
precision. A comparison is valid only if all solvers consume the same
`A`, `b`, `x0` and near-nullspace `B` (where supported).

The script does **not** calculate spectral condition numbers or estimate
smallest singular values. A small Rayleigh quotient alone is not proof
of near singularity for a nonsymmetric matrix. It does not solve the
nonlinear contact complementarity problem.

## Computational constraints

The implementation uses Python standard-library dictionaries for the
sparse operator. Memory scales with the **expanded nonzero count**; a
symmetric input stored as one triangle is expanded in memory. For very
large contact operators, SciPy/PETSc-based analysis will use less memory.

The synthetic test suite runs in Linux CI using:

```sh
python3 -m unittest discover -s tools -p 'test_contact_matrix_audit.py' -v
```

This utility is independent of MFAMG PR #201, the contact-aware nodal
aggregation proposal, and all LIS solver option-profile PRs.
