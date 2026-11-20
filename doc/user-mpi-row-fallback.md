# Native distributed USER row fallback

The `user_get_row` callback returns zero-based **global** column indices.
With two or more MPI ranks, when explicit `user_matvec` or
`user_matvech` is absent, LIS can compute `Ax` and `A^H x` directly
from the application-owned local rows. This is a correctness-first fallback:

- `Ax`: each rank requests only referenced remote columns,
  and owners return their input entries using `MPI_Alltoallv`.
- `A^H x`: local row contributions are sent to the owning column ranks
  along the reverse sparse request edges.
- Both variants synchronize buffer/argument failures before communication,
  and query/fill callback errors afterward. A failure on one rank returns
  a nonzero LIS error on all ranks. All ranks must call the operation.
- Row query and fill count consistency and global column ranges are
  checked. Scratch storage is proportional to local row data, incoming requests
  and the number of MPI ranks.

**Sparse exchange is rebuilt at every matvec call.** Unlike the original
full-vector fallback of PR #127, the implementation in PR #130 transfers
O(off-rank nonzeros) values and column requests rather than replicating
the full global vector. Repeated remote columns are not yet deduplicated
and request-count metadata remains O(ranks). For maximum performance on
large systems, explicit communication-aware USER matvec callbacks or
assembled CSR with its reusable sparse halo plan may still be preferable.
See [sparse MPI exchange](user-mpi-sparse-halo.md).

The distributed row fallback does not enable distributed direct-row
ILU-style preconditioners or automatic USER-to-CSR materialization.
Jacobi diagonal extraction from row callbacks is supported in the
stacked MPI Jacobi PR #128. It does not alter the public callback ABI.
Serial and explicit-callback behavior are unchanged.

`test/testusermpirows.c` exercises cyclic off-rank couplings (so a
diagonal-only smoke test cannot pass by accident), the Hermitian
transpose, GMRES numerical recovery, and rank-local callback failures
with collective status propagation, on 1 and 2 ranks and 32/64-bit
indices with sanitizer coverage.
