# Generic RCI v1.0 reconstructed Serial/OpenMP release candidate

Decision (2026-11-27): freeze the tested, optional real-double implementation as a **release candidate**, not an official LIS release or a recovered original RC3. No merge to master/develop, no release tag. MPI qualification explicitly excluded.

## Frozen code and evidence

Integration baseline: 7481b46ee7da9db7ff3dcbb3bedf99c50821860f.
Industrial diagnostic: c9df5cec63418ff3c8632e01bdb6a609ad0d84b9. The four implementation/header Git blob IDs match exactly between these commits:

| Path | Git blob SHA-1 |
|---|---|
| include/lis_rci.h | 15ba43b218c4bbb86553d086825ddf3d68a6627d |
| include/lis_rci_bridge.h | 4dfb1c0fc0f32fbb1458827c4fa5623dc3a4056d |
| src/lis_rci.c | 249824ee24b7edb69d18243a3cbc66f6bb549da3 |
| src/lis_rci_bridge.c | fd02f91fef986e0b4e4058d76f95b5164ae29eb8 |

- [Final integrated CI](https://github.com/3rav/lis/actions/runs/45941553504): 7/7 jobs; 288 solver runs including 96 runs with ASan/UBSan; default-disabled build verified.
- [Final ABI/legacy regression](https://github.com/3rav/lis/actions/runs/45941553616): 2/2 jobs, INT32 and INT64; existing public headers unchanged except generated configuration; no old exported symbols removed, exactly ten new RCI exports, unchanged SONAME. Six legacy executables per index width compiled against old library and run against both versions (24 executions); new C++ consumer passed both widths. This is a targeted compatibility test, not a complete cross-platform ABI audit.
- [Structural baseline](https://github.com/3rav/lis/actions/runs/45941007032): BCSSTK14 4/4 with Jacobi after scaling, BCSSTK17 0/4 within 2000 iterations.
- [SSOR/ILUT](https://github.com/3rav/lis/actions/runs/45941235881) and [SSOR all methods](https://github.com/3rav/lis/actions/runs/45941508377): BCSSTK17 CG passed with SSOR, remaining methods did not pass at the prescribed iteration limit.
- [Independent shifted Cholesky diagnostic](https://github.com/3rav/lis/actions/runs/45941680542): BCSSTK17 4/4, fixed SPD preconditioner from CHOLMOD with shift 1e-6, same rtol 1e-12, maxiter 2000, restart 50 and independent original-coordinate residual <= 1e-10, forward error <= 1e-5.

| Method | Iterations | Original relative residual | Relative forward error |
|---|---:|---:|---:|
| CG | 6 | 1.612867937e-14 | 1.191369109e-11 |
| BiCGSTAB | 3 | 1.118375375e-13 | 1.302843789e-08 |
| GMRES | 5 | 5.478815665e-13 | 1.060258752e-08 |
| FGMRES | 5 | 5.478877631e-13 | 1.060258766e-08 |

CHOLMOD is an external diagnostic dependency only, not a built-in LIS preconditioner or a new liblis dependency. The strong factorization has a reported fill estimate of 1,043,601 entries. The above is evidence of RCI compatibility with a suitable externally supplied preconditioner, not evidence that weak built-in preconditioners were repaired or that all industrial matrices will converge.

## Build and scope

```sh
autoreconf -fiv
./configure --enable-rci --enable-shared --enable-omp --disable-mpi --disable-complex --disable-longdouble --disable-quad --disable-fortran --disable-f90 --disable-saamg
make -C src -j2
make -C include install
make -C src install
```

Omit OpenMP or use --disable-omp for serial; --enable-longlong selects 64-bit LIS indices and requires -D_LONG__LONG when compiling the client. Compile with matching generated configuration header and -DHAVE_CONFIG_H; call LIS_RCI_CHECK_CLIENT_ABI() before typed LIS calls. Use lis_rci.h for the generic core and lis_rci_bridge.h for the identity/Jacobi bridge. CG requires an SPD operator and a fixed SPD preconditioner. BiCGSTAB and GMRES assume a fixed linear preconditioner; FGMRES permits variable preconditioning. Contexts are not reentrant; caller owns callback correctness, and borrowed matrices must remain valid and unchanged while solving.

API/ABI freeze on tested Linux x86-64, real double, default C enum and structure layout: options 48 bytes, request 24, info 40; all aligned to 8 bytes. Ten exports unchanged from tested baseline. Do not reorder existing fields, change enum values or function signatures in compatible updates. No claim of equivalence to unavailable original RC2/RC3. No native MPI, extended/complex precision, other platform or broad production workload qualification. Release candidate is opt-in and has no production merge or publication approval.

Diagnostic source files are preserved under qualification/rci-reconstructed. The CHOLMOD workflow can be inspected and manually dispatched where GitHub permits it. Historical logs remain attached to the exact tested commits.
