# Retained-anchor FGMRES pilot (experimental)

This opt-in application-side driver is a reference implementation, **not an installed LIS API or a change to the default solver**. It uses native LIS FGMRES and complete-fill ILUTP, and keeps one anchor factorization for a sequence of related matrices. It is intentionally limited to unscaled, assembled serial real-double CSR and one active thread. Matrices and vectors must share the same DOF ordering. No callback, MPI, complex, QUAD, or recycling support is claimed.

## Policy

1. Factor the application-provided anchor P once; keep P alive and immutable. Published structure or value generation changes invalidate the anchor.
2. Solve the current physical system A x = b from zero with the retained P and a configured pilot limit of 20 iterations. Accept only a finite independently computed RHS-relative physical residual at most 1e-8 (absolute residual for b=0), with a native numerical status.
3. If MAXITER and the projected residual ratio r20/r10 is at most 0.25, allow a further 80 iterations from the current iterate with the same P. This continuation restarts FGMRES; the ratio is a heuristic, not a certificate. Reported native iteration counts are retained as returned by LIS (they may include a terminal iteration counter increment).
4. If the result is still unacceptable, build one fresh factorization of A in a separate temporary solver and allow 100 iterations from zero. Never rebuild the same frozen factorization automatically. Do not retry hard API errors.
5. Commit x only on physical acceptance. Preserve the native raw status in the report, including MAXITER or BREAKDOWN with a separately accepted physical residual. Report successful factor creations separately from solve status. The anchor stays allocated during fallback, so peak memory can be higher.

The driver owns its solvers and preconditioners. The application owns P, A, b and x; P must outlive the driver. A may be P; b and x must not alias one another or matrix storage. Do not use the driver concurrently or mutate any input while solving. Recreate it after a hard solve error. The fixed profile and pilot limits deliberately reproduce the studied configuration, not an automatic universal tuning policy.

## Build and regression

Build LIS in a serial real-double configuration, with OpenMP optional. From the source tree after an in-tree static build:

```sh
cc -O2 -fopenmp -DHAVE_CONFIG_H -Iinclude tools/lis_anchor_pilot.c \
   tools/test_anchor_pilot.c src/.libs/liblis.a -lm -o test-anchor-pilot
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ./test-anchor-pilot
cc -O2 -fopenmp -DHAVE_CONFIG_H -Iinclude tools/lis_anchor_pilot.c \
   tools/bench_anchor_pilot.c src/.libs/liblis.a -lm -o bench-anchor-pilot
```

Omit `-fopenmp` when linking a non-OpenMP build. For an out-of-tree build use the generated configuration include directory and library from that build. The regression exercises retained-factor use, successful continuation, fallback, fresh reference, anchor generation rejection and output preservation. The dedicated workflow builds release and ASan/UBSan profiles; an unsupported configuration is not counted as a pass.

## Reproduction on application data

Pass any ordered list of native CSR Matrix Market matrix files and real general Matrix Market array RHS files (one column). Use exactly the same order, data, build and thread settings for both modes:

```sh
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ./bench-anchor-pilot fresh A1.mtx b1.mtx A2.mtx b2.mtx
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ./bench-anchor-pilot pilot A1.mtx b1.mtx A2.mtx b2.mtx
```

Both modes independently certify the physical residual and report setup+solve time, excluding input loading. Factorization counts include the initial anchor. First pilot setup timing includes driver construction; for performance comparisons prefer total time. The executable returns nonzero on any failed solve. An application can use the functions in `tools/lis_anchor_pilot.h` directly by compiling the reference source with LIS; they are not exported by the library build.

## Prior exploratory evidence and limitations

On 15 private CalculiX contact systems (60,265 DOFs; 12 distinct matrices), an earlier application-side experiment observed 15/15 physically accepted solutions, six factorizations instead of 15, and 221.66 s versus 296.61 s in an earlier fresh-reference run (25.27% lower time; 1.34x). The maximum certified relative residual was 2.94e-10. Peak RSS increased from approximately 844 to 1136 MiB. One direct kernel result reported MAXITER despite physical acceptance. Those runs used different harnesses and single observations; they do not prove a general speedup or nonlinear CalculiX convergence. The new benchmark provides a matched fresh/pilot harness for reproducible follow-up; private inputs are not included. Do not make this heuristic a default without independent models and repeated timing comparisons.

## Matched reference-driver smoke benchmark

The new executable was compiled against the parent of this package (`18219f23`), real-double serial/OpenMP, one active thread, `-O2 -DNDEBUG`. On original snapshots 01 and 15, using the same input files, profile, executable and physical acceptance rule in both modes, fresh mode accepted 2/2 with two factorizations in 21.883959 s; pilot mode accepted 2/2 with one factorization in 12.885392 s (41.12% lower time, 1.70x observed ratio). Snapshot 15 used the continuation path, achieved relative residual 8.843346254176e-14, and required no new factorization. Fresh mode achieved 1.512875194850e-13 on that snapshot. These are sequential single measurements in a shared environment, excluding input loading; this pair was selected because previous experiments showed reuse potential. They are not an independent full-series validation.
