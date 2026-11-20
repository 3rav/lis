# Distributed row-only USER: diagonal and Jacobi

This patch extends [the MPI row-backed USER fallback](../doc/user-mpi-row-fallback.md).
With only `user_get_row` installed, `lis_matrix_get_diagonal` reads
each rank's local rows (zero-based local row index) and selects the
diagonal via the **global** column ID `A->is + local_row`.
Duplicate diagonal entries are summed. It uses no dense matrix or
global allgather of rows; scratch is proportional to the longest local
row.

The callback query/fill pair must be stable. The code rejects negative
or overflowing counts, invalid global column IDs and inconsistent
query/fill counts. Even if a callback fails on only one rank, all
ranks return a shared nonzero error through `MPI_Allreduce`.
Accordingly, *all ranks must invoke distributed diagonal extraction
collectively*.

When Jacobi is requested for distributed row-only USER, the AUTO matrix
policy can now choose DIRECT, avoiding the unsupported MPI
USER-to-CSR materialization route. Explicit USER diagonal callbacks
retain precedence. Serial AUTO/Jacobi fallback behavior is unchanged.

The underlying distributed row-backed Ax/A^H x still uses full
global-vector collectives (O(global_n) scratch and communication per
rank) and is not a sparse-halo implementation. This feature does not
enable distributed row-source ILU/ILUT, distributed USER-to-CSR
conversion or any callback ABI change.

Regressions run 1- and 2-rank MPI solves with Jacobi and AUTO,
off-rank cyclic couplings, and one-rank callback failures, invalid
column IDs and inconsistent query/fill row lengths.
