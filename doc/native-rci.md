# Experimental shared native RCI

Configure with \`--enable-native-rci\` (default off). The legacy
\`--enable-rci\` flag and legacy raw-vector API remain independent.
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
