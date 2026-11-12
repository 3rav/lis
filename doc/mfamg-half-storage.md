# MFAMG optional FP16 preconditioner storage

MFAMG honors the existing LIS solver `-precon_half` setting. It reuses
the IEEE 754 binary16 conversion routines used by ILU.

* `off`: original full-precision storage and computation.
* `on`: compress suitable **prolongation coefficients**, **block-Jacobi
  LU factors** and **Jacobi smoother diagonals** after hierarchy
  construction and factorization, retaining full precision per array
  when any value is unsafe to encode.
* `noprobe`: skip the per-coefficient 1% relative tolerance guard, but
  still reject infinities, NaNs and nonzero-to-zero underflow.

Conversion failures and allocation failures are nonfatal: the corresponding
array stays full-precision. Compression is disabled when the build reports
FP16 unavailable, including complex builds.

Resident prolongation coefficients, block-smoother LU factors and
scalar-/block-Jacobi smoother diagonals can be compressed. Coarse CSR
matrices, pivots, terminal dense factorization and work vectors remain in
LIS scalar precision.
Restriction, prolongation and triangular solves decode coefficients and
accumulate in the solver precision. This is **not yet full mixed-precision
storage for the entire multigrid hierarchy**. Like ILU's half mode, half
storage is opt-in, separate from solver arithmetic precision.

The half-precision update path rebuilds the hierarchy transactionally
instead of attempting numeric in-place refresh of released arrays.

## Testing

`testmfamg` now tests `off`, `on`, and `noprobe` for weighted-Jacobi
and block-Jacobi smoothers. The regressions inspect resident packed arrays
and compare V-cycle correction to the full-precision reference. The GitHub
Actions Linux `half` job configures with `--enable-half`; other jobs
exercise compatibility with default builds.

Further work: full-precision residual certification under harder matrix
condition numbers, robust fallback probing, and total resident-memory
benchmarks for large problems.

## Diagonal safety and memory accounting

After the Jacobi diagonal has been computed and validated in solver
precision, it may be packed in FP16. The `on` policy guards the relative
error of its **reciprocal** (as used in the smoother); `noprobe` skips
that tolerance but still forbids zero/NaN/Inf. If a level's diagonal
cannot be packed, it remains at the original precision. This applies
to both weighted Jacobi and block-Jacobi's scalar fallback.

A read-only internal diagnostic, 
`lis_precon_mfamg_get_half_coefficient_bytes`, reports the full-
precision equivalent and current resident bytes of **only three
eligible coefficient categories** (prolongation, block LU, diagonals).
The test suite compares both counters, confirms no compression in
`off`/unsupported configurations, and verifies the PSD rebuild path.
These counters exclude matrix storage, work vectors, indices,
allocator metadata and overhead: they must **not** be interpreted as
whole-preconditioner resident memory or runtime measurements.
