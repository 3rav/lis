# Retained-anchor FGMRES pilot (experimental)

This opt-in application-side driver is a reference implementation, **not an installed LIS API or a change to the default solver**. It uses native LIS FGMRES and complete-fill ILUTP, and keeps one anchor factorization for a sequence of related matrices. It is intentionally limited to unscaled, assembled serial real-double CSR and one active thread. Matrices and vectors must share the same DOF ordering. No application-defined callback, MPI, complex, QUAD, or recycling support is claimed.

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

## Optional low-memory policy

Use `lis_anchor_pilot_create_low_memory(P, &ctx)` instead of `lis_anchor_pilot_create` or select `pilot-low-memory` in the benchmark. The default retains the original anchor across fallback. The optional mode destroys the anchor solver and factor before building a temporary fresh factorization; the original immutable P is still borrowed. On the next call it checks the original P generations before rebuilding its factor lazily. Reported factor counts and setup time include lazy reconstruction. The numerical policy, physical acceptance, output ownership and hard-error handling remain unchanged. The application must keep P alive and immutable even while its factor is absent. This is a trade-off, not a hard memory cap: allocator retention, matrix storage and temporary factorization workspace also contribute to RSS.

A selected four-step contact sequence (original snapshots 01, 02, 04, 15) was run with the same executable, profile, parent build `18219f23`, one OpenMP thread per process and three separate CPU-affinity cores concurrently in a shared environment. All processes exited 0 and accepted 4/4 solutions with max relative residual below 1e-8. Input loading is excluded from total setup+solve time.

| Mode | Factorizations | Total (s) | Peak RSS (KiB) |
|---|---:|---:|---:|
| Fresh | 4 | 77.754183 | 595092 |
| Pilot | 2 | 60.951973 | 979212 |
| Pilot low-memory | 3 | 70.584457 | 637724 |

The low-memory mode reduced peak RSS by 34.87% versus the retained-anchor pilot and reduced total time by 9.22% versus fresh. The extra anchor reconstruction made it slower than the normal pilot. The two pilot modes returned identical printed physical residuals. This is one selected subsequence, not a full-series or independently replicated result. The regression additionally exercises lazy reconstruction, subsequent reuse, and invalidation while no factor is live. No new installed API or library ABI changes.

## Optional rolling-anchor policy

Use `lis_anchor_pilot_create_rolling(P, &ctx)` or benchmark mode `pilot-rolling` to retain the most recent successful fallback factor rather than rebuild the original anchor. The old factor is retired before fresh setup, limiting live factor storage to one set. Neither matrix is owned by the driver. Keep the old anchor and current A alive and unchanged throughout the synchronous call. After **every** return, use `lis_anchor_pilot_current_anchor(ctx)` to determine which matrix must remain alive and immutable. A successful accepted fallback promotes A only after the output copy succeeds; the old matrix can then be released. A rejected numerical result or a hard failure never changes the anchor identity. The old factor may be absent, and a hard solve error invalidates the context; recreate after a hard error. A numerical rejection may lazily rebuild the old factor on the next solve. This is not a hard memory cap and remains an opt-in experimental policy.

The local rolling prototype was previously measured on the full 15 original contact snapshots in the same executable against low-memory, one thread per process on separate cores, excluding input loading: both accepted 15/15, rolling required nine factor creations rather than eleven, 253.689640 s vs 267.163146 s (5.04% lower observed time), peak RSS 821964 vs 821588 KiB. Printed physical residuals matched, including a physically accepted raw MAXITER. This is a prior single comparison, not a new benchmark of the hardened implementation or independent-model validation. Standard retained-anchor and low-memory policies are unchanged.

The extended regression tests promotion, releasing the former anchor, retaining the new factor and generation invalidation. The dedicated workflow also compiles `test_anchor_pilot_faults.c` using GNU ld wrapping and runs five deterministic boundary-failure cases (fresh factorization, fresh solve, generation lookup, output copy, numerical rejection). They require unchanged anchor identity and output, release of the rejected matrix, and correct subsequent invalidation or lazy rebuild. These tests are not exhaustive internal allocation-failure coverage.

## Optional single-pass and matched-restart controls

These are opt-in reference-driver policies, not installed application APIs or
changes to the default retained-anchor algorithm. The internal synchronous
FGMRES stop hook is restricted to the serial real-double CSR/ILUT path; it
adds no solver-object fields or global callback state. Other configurations
are not supported by this hook.

- `lis_anchor_pilot_create_single_pass(P, &ctx)` / `pilot-single-pass`:
  preserve the Krylov basis across the 20-step checkpoint if the projected
  residual ratio r20/r10 <= 0.25, with the original 100-iteration hard budget
  and restart length 60. If the checkpoint rejects, reconstruct the current
  approximation before returning a numerical stop. Physical acceptance and
  fresh fallback remain unchanged.
- `lis_anchor_pilot_create_restart_matched(P, &ctx)` /
  `pilot-restart-matched`: retain the 20+80 restart structure but normalize
  the continuation tolerance against the independently measured pilot
  residual, and reset the original tolerance for every new solve. This
  isolates the restart effect from the original driver's stricter effective
  stopping tolerance.
- These two modes count completed projected anchor iterations; legacy raw
  terminal counters may include an extra increment. Fresh counters remain
  native. Never infer nonlinear application convergence from linear residuals.

The opt-in environment variable `LIS_ANCHOR_AUDIT_DIR` names an existing
directory for exclusive creation of accepted solution files
`x-001.mtx`, etc. Run audits separately from timing comparisons. The
diagonal Python regression independently recomputes physical residuals in
70-digit Decimal arithmetic, checks raw decisions/iteration counts across
paired runs, and hashes input and executable files. Use fresh output
directories for each run.

```sh
cc -O2 -fopenmp -DHAVE_CONFIG_H -Iinclude tools/lis_anchor_pilot.c \
   tools/test_anchor_single_pass.c src/.libs/liblis.a -lm -o test-anchor-single-pass
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ./test-anchor-single-pass
python3 tools/anchor_stage9_diagonal.py ./bench-anchor-pilot ./diagonal-evidence
python3 tools/anchor_stage10_contact.py ./bench-anchor-pilot ./contact-evidence \
  --width 100 --height 100 --reps 3
```

The latter generates a sparse spring network with changing penalty constraints
and nonsymmetric rank-one perturbations. It is an algebraic contact
surrogate, not a finite-element contact simulation. Both scripts use one
CPU-affinity core, separate timing/audit runs, and round-trip decimal inputs.
The sparse checker uses independent streaming binary64 products and
`math.fsum` accumulation. Neither script requires private CalculiX inputs.

### Native experimental evidence (2026-11-26)

Native Release and ASan/UBSan regressions passed on experimental commit
`3503ba04fab2bead2108e7e9c25d13e4592c1c76`:
https://github.com/3rav/lis/actions/runs/45765554090 .

The diagonal sequence (15 states, 80 DOF, three repetitions per policy/profile)
accepted all 540 timing/audit solves, with 270 independently certified audit
vectors and maximum relative residual 9.211e-9. Corrected anchor-iteration
totals were 445 / 400 / 350 for classic / matched / single-pass.
https://github.com/3rav/lis/actions/runs/45766919827 .

The sparse surrogate Release results below are medians of three repetitions,
setup plus solve excluding input loading. Each listed size independently
certified 270 audit vectors over six policies; no pilot fallback occurred.

| DOF | Fresh (s) | Classic (s) | Matched (s) | Single-pass (s) |
|---:|---:|---:|---:|---:|
| 1024 | 0.087187 | 0.059847 | 0.038133 | 0.036151 |
| 10000 | 3.102383 | 1.452082 | 0.870498 | 0.877557 |

At 1024 DOF, corrected anchor iterations were 534 / 366 / 314; at 10000
DOF, 540 / 306 / 296. The sanitizer 1024-DOF case also passed (90
independent audit vectors). Source:
https://github.com/3rav/lis/actions/runs/45767929726 .

### Extended experimental results

The 60000-DOF Release contact-surrogate sequence also passed all 90
independent audit vectors; maximum relative residual was 8.9764e-9.
This is one timing repetition, not a statistically established speedup.

| Policy | Total (s) | Setup (s) | Corrected anchor iterations | Factor creations | Peak RSS (KiB) |
|---|---:|---:|---:|---:|---:|
| Fresh | 64.540514 | 63.117776 | 0 | 15 | 318144 |
| Classic | 24.344149 | 3.998087 | 469 | 1 | 310968 |
| Low-memory | 24.096604 | 4.008687 | 469 | 1 | 311004 |
| Rolling | 23.981713 | 4.004637 | 469 | 1 | 311040 |
| Matched | 16.052933 | 4.023242 | 291 | 1 | 311052 |
| Single-pass | 16.920309 | 3.955768 | 288 | 1 | 311056 |

Fresh additionally performed 15 fresh iterations; all five pilot strategies
had zero fallbacks on this particular sequence. Source:
https://github.com/3rav/lis/actions/runs/45767929726 .

The separate fixed-budget screening tested gates 0.10, 0.25, 0.50 and 0.75
at restart 60, and restart lengths 40 and 100 at gate 0.25, on normal and
abrupt-change sequences at 1024 and 10000 DOF. It included 13 candidates
(fresh plus 12 policy/parameter combinations), three repetitions, with
2340 independent audit vectors. All candidates passed; there were no
reported failures. The 1024-DOF abrupt-change single-pass case used
4/3/2/0 fallbacks for gates 0.10/0.25/0.50/0.75 respectively, and
223/241/257/294 total capped iterations. The observed median total times
were 0.027397/0.026225/0.024605/0.021070 s. The same gate choices did
not change iterations or fallback decisions at 10000 DOF. This is
insufficient evidence for changing the fixed 0.25 gate.
https://github.com/3rav/lis/actions/runs/45768177552 .

The separate adaptive holdout tested early stagnation rejection and
projected geometric continuation prediction in addition to the existing
controls (five candidates, four size/scenario combinations, three
repetitions, 900 independent audit vectors). All passed, but neither
heuristic changed a checkpoint decision on these fixtures. Therefore
no adaptive policy is promoted and no speedup is attributed to it.
https://github.com/3rav/lis/actions/runs/45768395227 .

The initial public-matrix benchmark could not download its SuiteSparse
archives (SSL connection timeouts); no numerical result was obtained.
https://github.com/3rav/lis/actions/runs/45768298336 .
A separate retry uses NIST Matrix Market BCSSTRUC2 data. Its results
are not yet included.

The unchanged default gate, restart and budgets are retained. The integrated
optional policies passed the Release and ASan/UBSan CI, including the paired
diagonal audit:
https://github.com/3rav/lis/actions/runs/45768507080 .
Full private CalculiX revalidation remains outstanding. The results are
not grounds for making either new policy a universal default.
