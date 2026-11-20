# Sparse MPI exchange for row-only USER

This patch replaces the correctness-first global-vector implementation
in #127 and is stacked after the MPI diagonal/Jacobi support in #128.

**Contract:** Each USER row callback accepts a zero-based *local* row
number and returns zero-based *global* column IDs. Every rank must enter
`lis_matvec` or `lis_matvech` collectively on the same matrix communicator.
Explicit application-owned `user_matvec`/`user_matvech` callbacks still
take precedence; serial USER kernels are unchanged.

## Algorithm

1. Fetch and validate local rows. Reject negative/overflowing nonzero
   counts, changed query/fill counts and out-of-range global columns.
   Synchronize callback/allocator errors across MPI ranks *before*
   communication.
2. Group the off-rank global column requests by owner using the
   existing matrix row partition ranges. Exchange request counts with
   `MPI_Alltoall` and column IDs with `MPI_Alltoallv`. Validate that
   every received ID belongs to the receiving rank.
3. For `Ax`, owners return only requested entries of their local
   input vectors with `MPI_Alltoallv`. Each rank combines local
   entries with the fetched remote entries. For `A^H x`, each rank
   sends conjugated row contributions along the same request edges
   in reverse, and owners sum received contributions.
4. Save the local input vector before updating any output, allowing
   direct in-place matvec and matvech on this MPI fallback.

No full global vector is allocated, and no `MPI_Allgatherv` or
global-length `MPI_Allreduce` remains in the sparse matvec path.
Storage scales as O(local_n + local_nnz + ranks + incoming_requests),
and vector traffic scales with O(off-rank_nnz) entries rather than
O(global_n) values. MPI metadata still incurs O(ranks) costs per call.
Repeated remote indices are currently **not deduplicated**.

The request plan is rebuilt per call, intentionally supporting USER
row patterns that change without a structure-generation notification;
future work can add a versioned, owned cache. Distributed ILU and
USER-to-CSR conversion policy remain independent concerns.

Error reductions are collective at each failure-prone step. Passing
a callback error on one rank must return the same nonzero status on
all ranks rather than deadlocking in the subsequent exchange.

## Regression

`test/testusermpisparse.c` exercises 1, 2, 3 and 4 ranks, a
diagonally-dominant cyclic matrix with duplicate remote columns,
forward/transpose numerical parity, a GMRES solve, in-place products,
runtime changes to diagonal-only and empty rows, callback failures
from a single rank and malformed row data. The CI matrix covers
32-/64-bit `LIS_INT`, Release and ASan/UBSan.
