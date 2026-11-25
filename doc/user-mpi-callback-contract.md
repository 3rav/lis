# Distributed LIS_MATRIX_USER: application-owned callback contract

`LIS_MATRIX_USER` is an operator interface, **not** a distributed sparse
storage format. Explicit `user_matvec` and `user_matvech` callbacks can
operate with an MPI-distributed LIS matrix, but LIS does not infer a halo,
exchange off-rank column values, or perform a distributed transpose reduction
on behalf of a USER callback.

## Local vector and index semantics

Create an MPI-distributed matrix with `lis_matrix_create(comm,&A)` and
`lis_matrix_set_size(A,0,global_n)`. Use `lis_matrix_get_size()` and
`lis_matrix_get_range()` to discover the process's local size and global
row interval `[is,ie)`. `LIS_MATRIX_USER_MATVEC` receives pointers to
the local vector storage used by LIS. The callback must fill the **local**
output segment, return `LIS_SUCCESS` on success or a LIS error code on
failure, and not assume that remote input entries have been fetched by LIS.

For off-rank couplings, the application's callback is responsible for
communication or another mechanism to obtain remote input values before
computing the local output. In an adjoint `A^H x`, the callback is likewise
responsible for any required interprocess accumulation. A diagonal-only
operator needs no such communication and is therefore a small but useful
smoke test of the interface.

The `user_data` pointer is borrowed. Its lifetime, numerical state,
consistent updates, and threading/MPI safety are the application's
responsibility. Do not mutate it or destroy the matrix concurrently with
an in-flight operation.

## What is not yet native

- The generic row-backed USER `Ax`/`A^H x` kernels and direct
  row-source preconditioners remain restricted to serial matrix access.
- Automatic cached USER-to-CSR materialization is serial-only.
- MPI-distributed row callback indices are defined as global column IDs,
  but LIS does **not** currently build an MPI halo/reduction schedule from
  those returned rows. Installing `user_get_row` does not activate such
  a distributed implementation.
- A successful custom MPI USER matvec test is not evidence that all
  LIS storage converters or preconditioners support distributed USER.

## Regression

`test/testusermpicallback.c` creates a positive diagonal operator whose
coefficients depend on the global equation number. Every rank owns its
local rows. It checks `Ax`, `A^H x` and CG/BiCG numerical solutions
on one and two ranks. A 2-rank diagonal test intentionally does **not**
exercise off-rank communications, collective callback errors, or ghost
buffer lifetime; these need separate applications and tests.

See `.github/workflows/user-mpi-smoke.yml` for reproducible MPI,
32/64-bit integer and ASan/UBSan build settings.
