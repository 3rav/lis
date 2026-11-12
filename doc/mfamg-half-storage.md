# MFAMG optional binary16 transfer storage

MFAMG now honors the existing solver `-precon_half` policy for resident
prolongation coefficients, using the same IEEE 754 binary16 converters as ILU.

* `off`: unchanged full-precision prolongation.
* `on`: attempt binary16 compression after the Galerkin hierarchy and
  terminal factorization are complete. Preserve full precision for a level if
  conversion overflows, flushes a nonzero coefficient to zero, exceeds a
  relative error of 1% (plus absolute tolerance 1e-7), or allocation fails.
* `noprobe`: attempt compression while omitting the rounding-error threshold;
  still reject overflow and nonzero-to-zero underflow.
* Complex builds retain the full-precision path.

Only prolongation coefficient arrays are compressed, not coarse matrices,
smoother factors, diagonal vectors, or dense terminal factorization.
Restriction and prolongation decode operands before accumulation in the
existing LIS scalar type. Numeric refresh with half mode enabled takes
the transactional rebuild path to avoid reusing discarded arrays.

This is an incremental implementation of FP16 storage for MFAMG, not
complete mixed-precision storage of all multigrid operators.

## Validation

Build and run `testmfamg` and `testmfamggalerkin` with the existing
repository test setup. Follow up with dedicated convergence and resident
memory tests for both `on` and `noprobe` modes.
