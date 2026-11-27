# Experimental native Generic RCI

The existing array RCI API is unchanged. The opt-in `lis_rci_native.h`
facade exposes six linear methods, FEAST and Krylov under one typed context.

## Build

Run `./configure --enable-rci --disable-mpi --disable-quad --enable-shared`,
then `make -j2` and `make install`. This candidate requires Python 3,
Node.js and GCC/Clang hidden symbol visibility. Disabled builds and installed
clients do not require the generators. Generation uses pinned local sources
in build directories, without modifying the original source or accessing
the network.

## Contract

Real double, non-MPI Serial/OpenMP; linear unscaled CSR with NONE/Jacobi,
FEAST standard/generalized CSR without retries or conversion. Krylov retains
its fixed-double callback ABI. Check the client/library ABI descriptor before
typed calls. Create a context with API version and request size, start the
appropriate method and execute requests while status is positive. Complete
each request with its token and native callback status. DONE commits linear
and FEAST solution outputs; Krylov retains native partial-result semantics.
Wrong tokens do not consume requests. Keep operands and solver alive and
exclusive; stop the backend before cancellation or destruction. All calls
must use the same thread and OpenMP settings. Shift signs are explicit:
FEAST sigma B-A, Krylov A-sigma B. See test_native_installed.c in
qualification/rci-restore for a complete client.

The native continuation is reconstructed from pinned LIS source, not recovered
RC3 code. This is a candidate, not a stable ABI. Extended scalar kernels are
separately qualified; this facade does not yet enable them. MPI and double-double
solver mode are excluded. Maintainer review and contributor licensing
confirmation remain necessary.
