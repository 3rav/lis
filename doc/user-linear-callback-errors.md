# Linear solver USER callback failures

This change concerns the standard-precision linear solvers and the
existing application callback ABI. It is not new MPI, QUAD or SWITCH
support.

## Status and output

A nonzero operator or preconditioner callback result must stop the
affected recurrence before the failed output is used. The callback code
is preserved in the solver's status; the existing public API convention
is unchanged: the caller must check both `lis_solve` and
`lis_solver_get_status`. The shared solve driver can still attempt a
post-solve residual check; this does not erase an earlier error.

- On failure before a usable residual is available, `resid` is NaN.
- Otherwise the last available recurrence residual estimate is retained
  where available. It is diagnostic, not a certification of the returned
  iterate. If a newly updated iterate could not be checked, the
  affected restart paths report NaN.
- `iter` describes solver progress at exit according to the individual
  method's existing iteration convention. Initial-residual failure is
  iteration zero; it is not a uniform count of completed outer steps
  across all algorithms.
- A failed solve does not promise a converged physical solution. Most
  methods retain the last available or partial physical iterate. In
  BiCGSTAB(l), a failed preconditioned correction is not a physical
  iterate, so the saved initial iterate is restored and residual is NaN.
  Do not accept any returned vector as a solution after a nonzero status.

The initial residual helper already records callback errors and the
solve driver clears the status before each execution on the current
develop base. For the default relative-to-initial-residual condition,
a finite nonzero initial norm is already normalized to one: the
additional diagnostic initialization reuses that fact rather than
performing another norm reduction. Alternate convergence conditions
may require an extra norm. No extra history entry or print operation
is added by the diagnostic initialization.

## BiCGSTAB(l) numerical fix and cleanup

The separate numerical issue included in this consolidated change is
that the initial physical iterate was copied to xp but was also left
in x, which is subsequently used for the preconditioned correction.
Final reconstruction could double-count the initial iterate. The
correction now starts at zero, and xp honors the zero-initial-guess
option. The native test `testbicgstablinit` is independent of the
callback-injection test `testuserbicgstabl`. The same change also frees
tau on the existing no-improvement MAXITER exit and callback error
exits. The correction is never reconstructed from failed preconditioner
output.

## Precision and MPI boundaries

CG (including its separately tested explicit-matrix QUAD path) and
MINRES changes were already accepted upstream in #254 and #253. This
change does not edit them or other QUAD/SWITCH implementations. Direct
matrix-free USER/OPERATOR solves are restricted to default precision
by the existing capability check.

The new callback tests are serial. This patch does not collectively
agree on rank-local callback error codes. For distributed application
callbacks, all participating ranks must follow a compatible collective
schedule and agree on success/failure before returning; a rank-local
error alone can leave peers blocked in a later collective. The same
requirement applies to distributed preconditioner callbacks. Do not
interpret this patch as support for rank-asymmetric failures.

## Regression and CI

The eleven native tests and their Automake targets are all in this
change; no dependency on the eigensolver or serial-contract PRs.
Fault injection covers selected initial and iteration applications
and selected restart/reconstruction paths, not an exhaustive proof
of every callback site or persistent recycle transition.

The workflow has two targeted PR profiles and four profiles on
matching develop pushes. The explicit paths intentionally exclude
test/Makefile.am and generated Makefile.in, following the CI budget
policy in #255. Unrelated test registrations must not launch this
specialized suite. Test sources and the workflow itself trigger it.
The existing baseline handles build-system registration changes.
OpenMP is enabled but this suite uses one thread; dedicated complex
and multithreaded testing belongs to the separate serial-contract
package.
