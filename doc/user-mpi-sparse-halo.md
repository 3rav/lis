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
2. Sort off-rank column occurrences by global ID, deduplicate them,
   and map every nonzero occurrence to its unique per-rank remote slot.
   The sorted IDs are naturally grouped by owning rank because the
   matrix row partition ranges are monotone. Exchange unique request
   counts with `MPI_Alltoall` and IDs with `MPI_Alltoallv`. Validate
   that every received ID belongs to the receiving rank.
3. For `Ax`, owners return each requested unique entry once using
   `MPI_Alltoallv`; every local row occurrence reads its shared remote
   slot while retaining its own matrix coefficient. For `A^H x`,
   each rank first sums conjugated row contributions to each shared
   remote slot, then sends one value per unique remote column along
   the same reverse request edges. Owners add received contributions.
4. Save the local input vector before updating any output, allowing
   direct in-place matvec and matvech on this MPI fallback.

No full global vector is allocated, and no `MPI_Allgatherv` or
global-length `MPI_Allreduce` remains in the sparse matvec path.
Storage scales as O(local_n + local_nnz + ranks + incoming_requests),
and vector traffic scales with O(unique_off_rank_columns_per_requester)
entries rather than O(global_n) values. Each requester counts a column
once even if it appears in multiple local rows. The per-call sorting
adds O(off-rank_nnz log(off-rank_nnz)) work and O(off-rank_nnz)
scratch; MPI metadata still incurs O(ranks) costs per call.

The request plan is rebuilt per call, intentionally supporting USER
row patterns that change without a structure-generation notification;
future work can add a versioned, owned cache. Distributed ILU and
USER-to-CSR conversion policy remain independent concerns.

Error reductions are collective at each failure-prone step. Passing
a callback error on one rank must return the same nonzero status on
all ranks rather than deadlocking in the subsequent exchange.

## Regression

`test/testusermpisparse.c` exercises 1, 2, 3 and 4 ranks, a
diagonally-dominant cyclic matrix with duplicate remote columns
within and across rows, forward/transpose numerical parity, a GMRES solve, in-place products,
runtime changes to diagonal-only and empty rows, callback failures
from a single rank and malformed row data. The CI matrix covers
32-/64-bit `LIS_INT`, Release and ASan/UBSan.
