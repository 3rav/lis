# Reconstructed historical RCI foundation (qualification only)

Based on the supplied historical transcript (v0.42, v0.46 and v0.51).
Not recovered original source; no claim of its ABI or numerical equivalence.

- Separate size/version-checked opaque handle request with A/B, PRECON,
  shifted solve, zero and complex AXPBY; generation and token guards.
- Block dispatch, with deterministic per-column fallback and error handling.
- Non-MPI contour projection, with complex shifts and weights and caller-owned
  vector storage, using the same machine for standard and generalized cases.
  The caller implements the backend, including shifted solves. No result is
  published on error; scratch can be modified and is not rolled back.
- Synchronous nested real-double solver adapter uses existing four core methods;
  it does not implement complex shifted solves or a native LIS preconditioner.
- Tests: independent diagonal complex reference, block/column dispatch, reuse,
  error/generation/cancellation injection at each phase, nested four methods,
  version, size and aliasing checks.

No public headers installed, no configure/Autotools changes. The source change
adds an internal qualification dimension accessor, to avoid reading beyond an
incompatible child context. Keep this branch separate from the upstream PR until
the new API design and qualification have been reviewed.

Not yet shared classic solver kernels, GCRO-DR/FGCRO-DR,
stackless Krylov or extended precision LIS builds. Complex
coefficients/opaque vectors are not equivalent to full complex LIS support.

The split-complex adapter now implements (zB-A) with two real planes and a
persistent doubled-real linear context. It supports BiCGSTAB, GMRES, FGMRES
with identity preconditioning, not CG (complex shifts need not be SPD).
Standard/generalized non-diagonal 3x3 tests compare the actual nested solves
and two-RHS contour projection against independent complex Gaussian
elimination. Input and output planes must be disjoint; the caller owns them.
Callback failure preserves the current column's output; earlier columns
may already have changed. No full complex LIS ABI support is claimed.

A separate native real-double adapter borrows assembled A/B matrices and a
prepared LIS_SOLVER/preconditioner. It uses lis_matvec and lis_psolve,
checks both structure/value generations and invalidates on replacement of the
prepared preconditioner. The adapter does not prepare arbitrary native
preconditioners itself; the caller must meet their workspace and lifecycle
requirements. Qualification tests use native Jacobi and SSOR, with all four
linear methods in Serial/OpenMP and INT32/INT64. Native bridge integration
does not imply a shared classic numerical recurrence or native FEAST.

The experimental FEAST outer continuation uses the split-complex contour
projector, B-metric twice-MGS rank reduction, native projected eigensolve
callback, Ritz reconstruction and explicit independent A/B residual
verification. It requires real symmetric A and SPD B, conjugate quadrature,
and a suitable starting subspace. It returns only converged interval values
after stable count/trace; it does not certify completeness or an empty
interval, and is not integrated into the classical FEAST policy. Every outer
request has a context-local monotonic token independent of nested projector
tokens. The native test uses lis_eigen_projected_solve and actual assembled
LIS A/B matrices (known spectra 1,2,5,8); tests repeat, cancellation, callback
error, generation invalidation and independently check physical residuals.

The recycle-space continuation stores fixed GCRO U/ZU/C and flexible
FGCRO W/Z/C separately (Z authoritative), refreshes images with A and
fixed preconditioning, performs twice-MGS with identical transforms,
drops rank-deficient candidates, and commits via a two-bank transaction.
Cancellation, incomplete output, generation mismatch or callback failure
keeps the previous committed basis. Residual and direction projection
preserve A*Z=C and orthogonality to C. Tests cover two RHS, changed operator,
fixed/flexible policy, and all request-phase faults. This is a building
block, NOT complete GCRO-DR/FGCRO-DR: native harmonic Ritz selection,
automatic integration into the augmented Arnoldi/restart recurrence and
differential native qualification remain open.

The augmented linear sequence now drives the existing GMRES/FGMRES RCI
recurrence cycle-by-cycle. It projects residuals and preconditioned
directions against a compatible committed recycle basis, retains the
U/ZU/C or W/Z/C augmented Arnoldi relation, and uses the actual native
GCRO-DR small QR/extraction/inverse-iteration routines for harmonic
selection. It commits transformed physical/image bases without new A/M
calls. Selection failure is nonfatal; stale operator/preconditioner
generations are conservatively not reused. The native dense selector can
allocate and is qualification-only. This is an experimental augmented
sequence, not yet differential equivalence to classical GCRO-DR/FGCRO-DR
or shared native classic recurrence. Initial coarse projection incurs one
extra true MATVEC before child verification; no duplicated Arnoldi MATVEC.

The native symmetric Krylov continuation is generated from the exact audited
upstream source blob f275ffbb5cf3ab0a6757d5e59eaefd9dd9f9b74c. Nine
functions and 28 suspend sites are transformed into persistent C11 frames;
all synchronous numerical helpers remain the original functions. The
generator deliberately rejects any other source blob and CI checks the
committed generated output. No native callback invocation, worker thread,
stack switching or replay in the generated path. The native shift convention
is A-sigma B. Qualification compares against lis_krylov_symmetric_ex64,
including callback input/output byte traces, results, native/extended stats,
errors, cancel and budget/workspace reuse. This is not yet integrated into
the installed library or the native classical entry point.
