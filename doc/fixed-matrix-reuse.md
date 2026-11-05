# Reusing a preconditioner for an unchanged matrix

This fork provides a deliberately restricted, **opt-in** session for solving multiple
right-hand sides with one prepared preconditioner.

## Supported configurations

- A fully assembled **CSR** matrix on a single MPI rank (serial or OpenMP).
- Default precision, `-scale none`, no matrix-storage conversion, additive
  Schwarz, or transpose-of-A option.
- **CG + Jacobi** (`-i cg -p jacobi`) or **BiCGSTAB + ILU**
  (`-i bicgstab -p ilu`).
- One solver and one matrix per session; the session must not be used concurrently.

Other combinations return `LIS_ERR_NOT_IMPLEMENTED`. Changes to the solver's
options or parameters, or to matrix dimensions and storage pointers, return
`LIS_ERR_ILL_ARG` until the session is rebuilt.

## C usage

```c
LIS_SOLVER solver = NULL;
LIS_SOLVE_REUSE ctx = NULL;
/* A is an existing assembled CSR matrix; b and x have matching dimensions. */
lis_solver_create(&solver);
lis_solver_set_option(
    "-i bicgstab -p ilu -scale none -tol 1.0e-11 -print none", solver);

if (lis_solve_reuse_create(A, solver, &ctx) != LIS_SUCCESS) {
    /* Report unsupported configuration or allocation error here. */
}
else {
    /* For each new right-hand side, keep A unchanged. */
    int rc = lis_solve_reuse_execute(ctx, b, x);
    /* Inspect rc and solver->retcode. */

    /* If numeric coefficients A->value changed in place without changing
       the sparsity pattern, explicitly rebuild the cached preconditioner: */
    rc = lis_solve_reuse_refresh(ctx);
    /* Check rc before another solve. */

    lis_solve_reuse_destroy(ctx);
}
lis_solver_destroy(solver);
```

**Lifetime:** Destroy the session before destroying its solver or matrix. The
session owns its `LIS_PRECON`, but borrows `A` and `solver`.

**Matrix changes:** This API **does not detect edits to existing numeric values**.
Call `lis_solve_reuse_refresh()` after any coefficient change, *before* the next
execution. If the sparsity structure, matrix type, dimensions, or storage pointers
change, destroy and recreate the session. Updating matrix coefficients without
refreshing would silently reuse stale factors.

**Solver changes:** Do not modify the solver's options or parameters until the
session is destroyed. The API deliberately refuses changed configurations.

**Failure handling:** The returned code is the kernel error, if any, otherwise
the solver's convergence status. The session retains its preconditioner across
an injected solver-work allocation failure, allowing a later retry.

## Benchmarks

Use `test/benchpreconreuse.c` to compare ordinary `lis_solve()` with a
manually prepared `LIS_PRECON` and repeated `lis_solve_kernel()` calls.
It verifies numeric equivalence and iteration counts. Results from shared CI
hardware are indicative, not guaranteed performance claims.

This API does not promise transparent reuse of SA-AMG, arbitrary preconditioners,
matrix updates without `refresh`, or multi-rank MPI calculations.
