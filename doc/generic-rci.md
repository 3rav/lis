# Experimental reverse communication interface (proposal)

Applications with an application-owned operator or preconditioner can drive
CG, BiCGSTAB, restarted GMRES or FGMRES without embedding their callbacks in
the existing LIS solver dispatch. The optional bridge runs the same interface
with an assembled LIS matrix and identity or Jacobi preconditioning. No existing
solver is replaced. This is a proposal for maintainer review, not a release or
a commitment to API stability.

Build with `./configure --enable-rci --disable-mpi --disable-quad --disable-complex --disable-longdouble`; add `--enable-omp` and/or `--enable-longlong` if desired.
Default builds do not install RCI headers or export RCI symbols. The experimental
module requires C11 and non-MPI real-double LIS. Unsupported enabled build
combinations fail during configure.

## Usage

Include `lis_rci.h` and create a context with the dimension, iteration limit,
restart length (GMRES/FGMRES), relative and absolute tolerances and method.
Call `lis_rci_start` with b and an optional initial guess. While it returns
MATVEC or PRECON, fill every entry of the request output with A times the input
or the application of the preconditioner, respectively. Call `lis_rci_advance`
with zero on success or an application error code on failure. Do not change the
input or reuse request pointers after advancing. DONE indicates that an
independent matvec verified the residual; negative status indicates failure.
The solution and diagnostic buffers remain owned by the context. Reuse a
finished context with another right-hand side of the same dimension.
See `test/testrci.c` for an executable example.

CG requires a symmetric positive definite operator and fixed SPD preconditioner.
BiCGSTAB and GMRES require fixed linear preconditioning; FGMRES allows changes.
GMRES and FGMRES share a right-preconditioned Arnoldi implementation storing
preconditioned basis vectors. Convergence is not guaranteed for arbitrary
operators and preconditioners. The bridge supports only identity and Jacobi.
No MPI, complex, long double, quad or other platforms are qualified.
A single context/bridge must not be used concurrently. Independent generic
contexts contain no shared mutable state; LIS bridge calls inherit LIS
thread-safety constraints. Allocation occurs at creation; callbacks own
their allocations. The diagnostic residual is absolute, not relative,
and can be stale on failure.

For the bridge, compile against matching installed LIS headers and flags,
including `-D_LONG__LONG` when appropriate; check `LIS_RCI_CHECK_CLIENT_ABI()`
before creating typed LIS objects. This check screens scalar/index widths and
is not a full layout proof. Keep the borrowed assembled matrix alive until
destroying the bridge; never mutate it during a solve.

## Origin and qualification

The sources are an independent AI-assisted implementation of standard iterative
methods developed in this contribution, not imported from the unavailable prior
RC3 or copied from an external RCI package. The proposed terms are those in the
repository's COPYING; the contributor must confirm the license before merge.
No third-party matrix data, opaque source archive or external numerical
dependency is included in this PR. Historical qualification (including
CHOLMOD-assisted structural cases and failed weak-preconditioner trials) is
retained in the contributor's fork, not part of upstream's test obligations.
The single CI workflow runs deterministic tests and compares existing
ABI exports against the current merge base. See the PR for evidence.

## Extended RCI and protocol

The original two-request API remains available. The extended interface
(`lis_rci_start_ex`, `lis_rci_advance_ex`) additionally externalizes every
dot product and Euclidean norm, including RHS norm and final residual norm.
All vectors are full-length real-double vectors. DOT supplies two read-only
inputs; NORM supplies one. Supply a finite long-double scalar on completion
(nonnegative for NORM); vector requests use their existing output arrays.
The caller is responsible for computing correct scalars, just as for MATVEC.
Each request has a nonzero uint64_t token, monotonic within one context across
restarts and solves. A wrong token returns PROTOCOL without consuming the
outstanding request. Tokens do not protect writes to expired pointers.
Do not mix legacy and extended entry points within a solve. The legacy API
executes the scalar requests internally using the existing local backend.
No vectors are allocated for individual requests.

An optional, transport-independent consensus helper packs five uint64_t words:
`[token, ~token, tag, ~tag, error_word]`. The tag is an opaque numeric
identifier, compared bitwise; it is not floating point equality. Componentwise
unsigned maximum combines records. An error word is zero for success;
otherwise its high 32 bits encode the inverse rank and its low 32 bits encode
the signed error. After merging a nonempty set with unique ranks, unpack reports
token/tag disagreement and the error of the lowest failing rank separately.
The caller chooses precedence. This is a pure local primitive, not an MPI
collective or distributed solver. It does not compare arbitrary long-double
scalar values, and it does not establish compatibility with the unavailable
original RC3 protocol/API.

`test/testrciext.c` drives all four solvers with external scalar operations,
tests stale tokens and reuse, and checks 8192 deterministic consensus scenarios
and a token difference above 2^53. The installed-library CI runs both legacy
and extended APIs.
