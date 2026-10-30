# Generic C API interoperability: contracts and validation

These APIs provide three optional C interfaces. Public headers are already installed
by `include/Makefile.am`.

## Nullspace preflight

`lis_solver_check_nullspace_configuration(solver, &issues)` returns
non-destructive configuration flags. The same checker is called by
`lis_solve_kernel()` to prevent drift between preflight and solve.
Zero issues do not imply the nullspace basis has the correct dimension or
that the numerical problem is solvable.

Near-nullspace acceptance is an OR of two flexible outer solvers and the
direct MFAMG preconditioner. The current kernel does not automatically treat
nested Hybrid+MFAMG as a direct MFAMG preconditioner. In this build,
`LIS_PRECISION_DEFAULT == LIS_PRECISION_DOUBLE`.

## Owned CSR

`lis_matrix_materialize_csr_owned(source, &csr)` creates a standalone,
LIS-owned CSR result. The caller must destroy it with `lis_matrix_destroy`.
The CSR does not alias borrowed source arrays or depend on a cache.
The source must support row access. Only a single-process matrix is
supported (including MPI with one rank); no distributed materialization.
If the output pointer is non-NULL, error paths set `*csr=NULL`.
Source callbacks must remain stable during the operation.

## Real64 vector transfers

`lis_vector_import_real64` and `lis_vector_export_real64` validate
the full input before writing the destination. They handle positive strides
and permit NULL data pointers for valid zero-length vectors.
NaN and Inf are rejected on import and export, including discarded imaginary
parts under ALLOW_LOSS. REJECT_LOSS refuses nonzero imaginary components and
loss of real precision. The implementation is intentionally serial for
validation-before-write semantics. It is a C interface: no Fortran wrapper
is claimed.

## AUTO sequence policy

The experimental in-place CSR value-refresh sequence reuse was **reverted**
because existing preconditioners can keep references to derived CSR
diagonal and D/L/U storage, which value refresh destroys.
Both DIRECT and CACHED_CSR therefore retain conservative reset-on-value
change. Reuse across value refresh requires a separate transactional,
preconditioner-lifetime-safe design, with benchmarks.

## Validation

Run `autoreconf -fiv`, configure and build LIS, and then run `make check`.
`testmatrixautosession` and `testvectorreal64` are in the actual
Automake `TESTS` list. For targeted regression coverage, run `testmatrixautosession`,
`testmatrixmaterialize`, `testcapabilities`, and `testvectorreal64`
directly after building. Recommended validation configurations include
Linux real/complex i4/i8, OpenMP, long-double, sanitizer, and MSYS2 UCRT64
i4/i8. CI configuration and run status are maintained outside this source PR.
