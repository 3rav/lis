# Temporary native LIS qualification

This isolated branch qualifies the independent real-double reconstructed RCI, not the lost original RC3. No MPI, no private matrices, no changes to existing solvers. The tarball is byte-identical to the five source files from the local integration package; SHA-256: 95fdea9d68dbf892b98878518ac61e5ae5b39d75bad72941120ab6402a413425. Build and smoke test against actual LIS with configured headers, 32/64-bit indices and serial/OpenMP. No official ABI or release claim.
