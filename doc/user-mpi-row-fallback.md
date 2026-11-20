# Native distributed USER row fallback

The `user_get_row` callback returns zero-based **global** column indices.
With two or more MPI ranks, when explicit `user_matvec` or
`user_matvech` is absent, LIS can compute `Ax` and `A^H x` directly
from the application-owned local rows. This is a correctness-first fallback:

- `Ax`: `MPI_Allgatherv` assembles the full global input vector on
  each rank. Local USER rows multiply their entries against it.
- `A^H x`: local rows accumulate into a global-length scratch vector,
  then `MPI_Allreduce` sums those contributions, and each rank extracts
  its local output segment.
- Both variants synchronize buffer/argument failures before communication,
  and query/fill callback errors afterward. A failure on one rank returns
  a nonzero LIS error on all ranks. All ranks must call the operation.
- Row query and fill count consistency and global column ranges are
  checked. Scratch storage is O(global_n + maximum_row_nnz + ranks).

**This is not a sparse halo implementation.** Every rank replicates a
full global vector and communicates O(global_n) data. This may be much
slower and more memory-intensive than CSR for large sparse distributed
systems. For scalable real-world MPI cases, applications should continue
to supply explicit communication-aware USER matvec callbacks, or use CSR
with LIS's normal sparse halo exchange.

The distributed row fallback does not enable distributed direct-row
preconditioners, automatic USER-to-CSR materialization, or diagonal
extraction from row callbacks. It does not alter the public callback ABI.
Serial and explicit-callback behavior are unchanged.

`test/testusermpirows.c` exercises cyclic off-rank couplings (so a
diagonal-only smoke test cannot pass by accident), the Hermitian
transpose, GMRES numerical recovery, and rank-local callback failures
with collective status propagation, on 1 and 2 ranks and 32/64-bit
indices with sanitizer coverage.
