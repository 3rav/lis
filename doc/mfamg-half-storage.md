# MFAMG optional FP16 preconditioner storage

MFAMG honors the existing LIS solver `-precon_half` setting. It reuses
the IEEE 754 binary16 conversion routines used by ILU.

* `off`: original full-precision storage and computation.
* `on`: compress suitable **prolongation coefficients** and **block-Jacobi
  LU factors** after hierarchy construction and factorization, retaining
  original full-precision arrays when any value is unsafe to encode.
* `noprobe`: skip the per-coefficient 1% relative tolerance guard, but
  still reject infinities, NaNs and nonzero-to-zero underflow.

Conversion failures and allocation failures are nonfatal: the corresponding
array stays full-precision. Compression is disabled when the build reports
FP16 unavailable, including complex builds.

Only resident transfer coefficients and block-smoother LU factors are
compressed. Coarse CSR matrices, scalar smoother diagonals, pivots,
terminal dense factorization and work vectors remain in LIS scalar precision.
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
