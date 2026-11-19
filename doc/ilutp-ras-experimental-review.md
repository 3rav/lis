# Experimental ILUTP and RAS improvements

This branch is deliberately independent of BSR/GCRO-DR correctness fixes.
It contains two potentially performance-affecting changes and one optional
stability diagnostic. It is not a claim of faster FEA solves.

## ILUTP fill selection

The bounded heap selects the same top-`cap` entries by
`(descending |value|, ascending column)` that the previous full qsort
selected in ordinary finite-valued cases. Retained entries are re-sorted
before storage. This changes the cost of limited-fill selection from
sorting all `m` entries to maintaining `cap` best entries.

Complexity is reduced in theory when `cap << m`; wall-clock speedup and
preconditioned-iteration equivalence remain to be validated on authentic
FEA matrices. Full-fill ILUTP retains the original path.

## Optional inverse-gain diagnostic

`LIS_ILUTP_GAIN_DIAG=1`: print estimated magnitudes of the triangular
inverse majorants after successful ILUTP factorization.

`LIS_ILUTP_GAIN_LIMIT=<positive finite number>`: reject a factorization
if its conservative inverse-gain upper bound is nonfinite or exceeds the
limit. The default behavior is unchanged if the variable is absent.
An invalid value produces `LIS_ERR_ILL_ARG`.

The bound is **not** an exact condition number or proof of convergence.
It only measures possible amplification through the incomplete factors;
it cannot predict contact nonlinearity or Newton convergence.

The environment configuration is temporary/experimental and may not be
appropriate for distributed or embedded applications. A documented solver
option should supersede it before enabling a default policy. Record these
environment variables in benchmark logs for reproducibility.

## Parallel RAS scheduling

At setup time, estimate each independent subdomain's triangular solve
cost by local factor NNZ. During application, process the more expensive
domains first with dynamic OpenMP scheduling. Each task uses its own
factorization and its own scratch vectors, and writes only to its
disjoint *owner core* in the output. Any Schur/interface correction occurs
**after** all subdomain tasks have joined the barrier.

This should not change the mathematical operator. Different scheduling
can still change floating-point execution order in internal routines;
repeated-run checks with identical residuals, tolerances and thread counts
are required before making reproducibility claims.

## Acceptance evidence required before release

Run a side-by-side comparison against upstream `develop`, pinning compiler,
thread count, solver and preconditioner options, pivot threshold and matrix
ordering. For each real CalculiX/FEA matrix and RHS, record:

- Factorization setup time and apply/solve time separately;
- `nnz(L)+nnz(U)`, peak resident memory and drop/pivot diagnostics;
- FGMRES/(F)GCRODR iterations and the **true**
  `||b-A*x||_2 / ||b||_2`;
- A complete Newton/contact time sequence, including convergence/failure,
  preconditioner rebuilds and fallback decisions;
- Relative differences in output for one and multiple OpenMP threads.

CI's small synthetic matrices only establish basic numerical safety;
they do not establish a performance benefit or production readiness.
