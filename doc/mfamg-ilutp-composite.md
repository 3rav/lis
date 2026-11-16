# Experimental MFAMG + native ILUTP residual composition

This opt-in *built-in* preconditioner combines an MFAMG V-cycle with a
single global native ILUTP correction. It is selected with
`-p mfamg_ilutp` in the current development branch.

## Algorithm

For input residual `r`, the composite applies

```text
z1 = MFAMG(r)
r1 = r - A*z1             (same assembled A used by the Krylov solver)
z2 = ILUTP(r1)
z  = z1 + z2
```

Both numerical preconditioners are built during setup and reused for each
application. Failures are returned to the caller; there is **no** silent
ILUTP-only or MFAMG-only fallback. This is multiplicative residual
correction, *not* ILUTP pivoting inside MFAMG blocks and *not* recursive
`Hybrid`. The outer FGMRES still verifies convergence on the original
system.

## Initial supported scope

- Real, double-precision, non-split, single-process explicit CSR matrices.
- `-scale none`; at least one caller-supplied MFAMG near-nullspace vector.
- Global ILUTP uses a serial pivot sequence: set `OMP_NUM_THREADS=1` **before**
  creating the solver/preconditioner (even if built with OpenMP).
- No distributed MPI ranks, complex arithmetic, quad arithmetic, or
  Additive Schwarz wrapper.
- No transpose/adjoint `psolveh` implementation.
- No `-hybrid_p mfamg_ilutp` nesting. Select `-p mfamg_ilutp` directly.
- No native RAS composition. Composite ILUTP uses a single global factor;
  `-ilut_ras_*` options are not forwarded into it.
- Numeric PSD/sequence preconditioner refresh is not implemented in this
  first version; use ordinary `lis_solve` rebuilding for changing matrices,
  rather than assuming old factors remain valid after a Newton update.

## Example (LIS option string)

```text
-i fgmres -p mfamg_ilutp -scale none
-mfamg_block_size 2 -mfamg_vcycle_smooth_steps 1
-ilutp_threshold 0.5 -ilut_equil 1 -ilut_rcm 1
-ilut_drop 0.005 -ilut_rate 8 -restart 80
-maxiter 1000 -tol 1e-8
```

The composite selects native ILUTP automatically, even without an explicit
`-ilutp 1` option. The ILUTP threshold, optional complete-fill cap, RCM,
row scale, Equil, and drop/rate are inherited from the outer configuration.
MFAMG inherits its block size, fine-aggregate factor, smoothing passes,
deep-level sibling merge, near-nullspace basis, and optional FP16 storage.
ILUTP factors themselves remain in full precision. Existing `-p mfamg`
and `-p ilut` behavior is unchanged.

The host application must call `lis_solver_set_near_nullspace()` with
compatible vectors. A constant vector is useful for diffusion-like models;
contact/mechanical systems often require physically meaningful rigid-body
modes or other application-selected near-nullspace vectors. This code does
not invent those vectors.

## Validation

Run `autoreconf -fiv` and `./configure --enable-test`, then
`make -j2 && make -j2 check`. The new `testmfamgilutp` regression exercises
`lis_solve` for MFAMG and the composite, repeated solves on a deterministic
real nonsymmetric CSR system, original-system solution error, and rejection
of missing near-nullspace data.

These tests **do not** establish a contact-model convergence advantage,
iteration speedup, race freedom under arbitrary concurrent calls, or safe
reuse across changing matrix generations. Benchmark on saved *identical*
CalculiX/contact `A,b,near-nullspace` snapshots with the same outer
FGMRES settings. Compare true physical residual, iterations, factor setup
time, total wall time, and peak resident memory against standalone ILUTP
and standalone MFAMG. The second full-matrix ILUTP setup can be expensive,
so better convergence need not imply better total runtime.

## ABI / rollout

This development change appends one built-in preconditioner ID (13),
moving `LIS_PRECON_TYPE_USERDEF` from 13 to 14 and extending the internal
`LIS_PRECON_STRUCT`. Recompile code against a consistent LIS header and
binary; do not mix artifacts built against the previous internal layout.
This PR is experimental and should not be merged without CI and real
contact-matrix verification.
