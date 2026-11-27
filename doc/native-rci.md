# Experimental shared native RCI

Configure with `--enable-native-rci` (default off). The legacy
`--enable-rci` flag and legacy raw-vector API remain independent.
The optional native session API uses LIS vector handles, except Krylov
which uses its original real-double array ABI. See lis_rci_native_session.h
for ownership, request versions, split-complex shift conventions and status.
The full lifecycle header is installed in real-double builds only.
Its documented restrictions are checked before native setup.
MPI, double-double solver mode, and fast/finite-only math are unsupported.
Complex and extended scalar builds support the prepared native session.
Krylov itself stays real double. Do not use private restore/test symbols.
This is experimental, not a stable new ABI. Generated continuations are
materialized from pinned native LIS source by the qualification generators;
normal builds require neither Python nor Node.

A native request can be evaluated with `lis_rci_native_execute(&request,
sizeof request)`, then its returned native callback code passed to the
session/full advance function with the unchanged token. This is optional:
applications may fill the outputs themselves, including by running a separate
child session for a shifted solve. The installed-API test covers independent
parent FEAST and child GMRES contexts and error propagation. Raw-double
Krylov requests remain application-defined.

## Scope and provenance

This is an experimental reconstruction of shared native continuations, not a
recovered copy of the original RC3. The installed API is opt-in and versioned.
The existing public LIS and legacy RCI entry points remain available.

The prepared session covers CG, BiCGSTAB, GMRES, FGMRES, GCRO-DR, FGCRO-DR,
FEAST and native Krylov. Full lifecycle setup/cleanup is available for supported
real-double linear and FEAST configurations only. The full lifecycle requires
exclusive borrowed objects, the same host thread throughout, CSR (or USER with
no preconditioner for linear solves), no scaling, storage conversion, nullspace
correction or per-solver thread override; FEAST requires CSR A/B and no automatic
subspace retry. Preconditioner construction remains synchronous. See the
installed headers for the precise contracts. A cancelled or failed solve does
not certify an iterate. Extended and complex qualification applies to prepared
sessions, not the full lifecycle; native Krylov retains its real-double ABI.

For a fresh build, use `./configure --enable-native-rci --disable-mpi --disable-quad`
with your usual LIS build options, followed by `make` and `make install`.
Run `make check` with tests enabled in a supported real-double build.
The new mode must not be used with fast-math or finite-only math.
