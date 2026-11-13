# Experimental serial RAS-ILUT and approximate interface Schur

LIS ILUT now has a **separate opt-in** serial restricted additive Schwarz (RAS) variant, inspired by the decomposition architecture of hypre. It does not depend on hypre or use hypre code.

## Configuration

```
-i fgmres -p ilut -scale none -ilut_ras_domains 4 \\
-ilut_ras_overlap 1 -ilut_schur 1 -ilu_pivot_tol 1e-12
```

- `-ilut_ras_overlap 0|1|2` turns on RAS (default: off). Overlap is the number of sparsity-graph edges added beyond each core partition.
- `-ilut_ras_domains 2..64` selects the count of deterministic, contiguous partitions (default: 2). Each subdomain uses the normal LIS ILUT factorization on its extracted local matrix, including the outer `-precon_half` policy.
- `-ilut_schur 0|1` adds an optional **diagonal-approximate interface Schur correction** (default: off). The interface is rows coupled across the original nonoverlapping partitions.

RAS applies the local factor solve on each overlapping subdomain but injects corrections into that subdomain's **nonoverlapping core** only. For the Schur variant, LIS computes

\[\widetilde S=A_{\Gamma\Gamma}-A_{\Gamma I}\operatorname{diag}(A_{II})^{-1}A_{I\Gamma}\]

and uses a small pivoted dense LU factorization to correct the interface residual. It then applies diagonal-approximate interior back-substitution. This is a **one-step approximate Schur correction**, not hypre's distributed GMRES-on-Schur algorithm and not a general-purpose exact Schur complement.

## Scalability experiments (not yet production-ready)

The experimental successor branch `experiment/ras-schur-scalable` extends the serial implementation with:

- `-ilut_ras_domains 2..64` (previously 2..16).
- `-ilut_ras_parallel 1`: explicitly opt-in OpenMP parallelism of independent **local preconditioner applications only**. Setup/factorization remains serial because LIS currently mutates shared global dispatch and error-reporting state during construction. Serial builds execute the same loop sequentially. Per-domain errors are collected and returned after the OpenMP join.
- Bounded weak-diagonal promotion from the interior to the Schur interface before factorization, with a cap of 20% of the matrix dimension. This is not sufficient to guarantee a well-conditioned Schur operator. A larger or singular interface still fails safely.

These changes are not a claim of automatic near-linear scaling. **The current experimental successor adds a sparse CSR approximate-Schur path above 256 interface DOF**, using a bounded rowwise accumulator, `nnz(S) <= 4 * nnz(A)`, relative dropping, and one existing ILUT preconditioner application with row equilibration. Small interfaces retain dense pivoted LU. This is not an iterative Schur solve and sparse fill can still be large; exceeding the budget fails explicitly. The separate prior [3rav/lis#64](https://github.com/3rav/lis/pull/64) remains a standalone research PoC. RCM-first partitioning and parallel factorization are **not implemented**. No production contact FEA performance comparison is available.

### Explicit synthetic sparse Schur smoke test

```sh
./test/testhalfilut --sparse-schur-smoke
```

The smoke test constructs a **512-unknown scalar 3D seven-point stencil** with four contiguous partitions and an expected 384 interface DOF. It checks that sparse Schur setup succeeds, reports true residuals after one preconditioner application, and compares results with `-ilut_ras_parallel 0` and `1`. It does **not** represent three-dimensional vector elasticity, a Newton/contact KKT matrix, or a measured performance advantage. CI includes it in the Linux FP16 job.

## RCM-ranked partitioning and reproducible Matrix Market comparison

`-ilut_ras_rcm 1` optionally partitions domains along the deterministic symmetric-graph RCM ordering shared with ILUT from upstream #224. The input matrix and global variable numbering are not changed. Core injection still uses disjoint owner indices, even when they are no longer contiguous. The default remains `0` (original index ordering) to preserve prior behavior. RCM is a bandwidth-reduction heuristic, **not a guaranteed minimum-separator partitioner**; interface size and timing must be measured for each FEA mesh.

The optional extended Matrix Market benchmark compares baseline ILUT, RAS, and RAS+Schur using the same original `Ax=b` and manufactured known solution:

```sh
OMP_NUM_THREADS=8 ./test/testhalfilut --fea-ras contact.mtx 8 8 8 0.005 8 1 1 1
```

The arguments after the file are: baseline fill rate, RAS fill rate, RAS+Schur fill rate, common drop tolerance, domain count, overlap depth, RCM-ranked partitions (0/1), and OpenMP local apply (0/1). The output includes **wall-clock** setup/solve/total times, FGMRES iterations, true original-system relative residual, interface DOFs and a stored-entry proxy. Each case prints `fair_fill` only when its proxy differs by at most 5% from baseline ILUT. Because factors, sparse CSR, indices, and workspaces have different costs, the stored-entry proxy is **not peak resident memory**: measure peak RSS independently before drawing performance conclusions. If a case fails, its status and LIS error are printed rather than inventing convergence results.

The benchmark is intended for **external real FEA snapshots**. The synthetic 512-DOF stencil test is only a correctness regression and is not substituted for a 50k–200k DOF contact elasticity experiment. No production FEA speedup is claimed.

## Synthetic contact-KKT regression

`testhalfilut` now assembles an indefinite **48×48 stabilized saddle-point system** with a 36×36 spring-stiffness block, twelve displacement-difference constraints, and a negative 0.1 multiplier diagonal. It compares the legacy ILUT control with the opt-in variant, requiring convergence to a **true relative residual below 1e-7** and printing Krylov iterations and residuals. This is an algebraic contact-style regression, **not** a production nonlinear FEA benchmark and not evidence of a speedup.

## Dependency update (2026-11-13)

This branch is synchronized with the latest ILUT robustness changes from fork PR #62 (upstream #224), including the high-unit-contrast KKT test and the English user-guide option table. The fork PR #63 remains stacked cleanly on #62, with its diff restricted to RAS/Schur changes. The upstream #225 comparison against `develop` necessarily includes #224 until that dependency merges.

An entirely **separate, non-integrated** sparse CSR Schur feasibility experiment lives in fork Draft PR #64. It checks a synthetic interface larger than 256 DOF and documents OpenMP risks, but it is not a fix for the production limitations of the current dense Schur implementation.

## Local ILUT RCM ablation and numerical finding (2026-11-13)

The RAS domain-ownership option `-ilut_ras_rcm 1` only orders how the global matrix is divided into subdomains. **It is separate from** `-ilut_rcm 1`: the latter now also orders each extracted local owner+halo matrix before factorization, using the existing symmetric RCM/row-scaling implementation. Both remain opt-in. The local ordering is applied after extracting the subdomain, not merely to the owner ranks.

A controlled real structural-FEA ablation on `Boeing/bcsstk39` (46,772 DOF, 8 OpenMP threads, 8 domains, overlap=0, domain RCM=off, row scaling=off, fill rate 8) detected extreme inverse amplification **in local domain 6** when local RCM was disabled: `||M_d^-1 b_d||/||b_d||=8.16e70`, while the smallest observed factor-pivot magnitude over all domains was about `6.98e4`. Thus absolute pivot magnitude alone did not detect instability.

When local `-ilut_rcm 1` was enabled, the worst local inverse response dropped to `2.23e-4`. The complete RAS preconditioner response changed from `2.25e70` to `6.26e-5` relative to `||b||`; at 500 FGMRES iterations the true relative residual improved from `0.999963` to `0.001772`. The baseline 8-thread ILUT reached `0.004015` at the same iteration count. **None met the required convergence tolerance** and the locally ordered RAS was slower than the baseline (roughly 6.62 s versus 5.46 s setup+solve, with 19% less stored-factor-entry proxy). Adding the sparse Schur correction worsened the local-RCM residual to `0.005588`.

[Complete two-configuration ablation and worst-domain factor diagnostic](https://github.com/3rav/lis/actions/runs/43397517228).

The result supports local reordering as a way to avoid a catastrophic failure on this one matrix. It **does not prove a general speedup**, stability across contact-KKT systems, or that any particular fixed threshold/diagonal shift would be safe. A default-on factor-stability guard is **not implemented**: dimensionless tests of triangular growth and regularization policies need independent validation. This branch remains Draft.

## OpenMP ThreadSanitizer and real structural-FEA findings (2026-11-13)

The OpenMP RAS local-solve path no longer enters nested ILUT parallel regions. Local ILUT factorization and application now consistently use one thread; the sparse Schur-ILUT setup and apply do the same. The number of threads used by a factorization must not change between setup and application, since the legacy OpenMP ILUT algorithm partitions factor dependencies by thread.

The dedicated **Clang/TSan + LLVM Archer** workflow passed on the branch with packed FP16 and row-scaled RAS smoke tests: [TSan run](https://github.com/3rav/lis/actions/runs/43392959868). Archer supplies OpenMP synchronization semantics through OMPT; earlier unsynchronized libomp/runtime reports were not assumed to certify application safety. The test covers parallel *application* of independent subdomains, **not** concurrent independent calls on the same preconditioner or parallel factor construction.

A reproducible public structural-FEA benchmark was executed on **Boeing/bcsstk39**, a shuttle solid-rocket-booster stiffness matrix, with 46,772 DOF and 2,089,294 Matrix Market entries after import. It is **not** a contact-KKT system. All configurations below used FGMRES (maximum 500 iterations), 8 RAS domains, overlap 1, RCM-ranked partitions, and nominal ILUT fill rate 8.

| Threads / row scale | Baseline true residual | RAS true residual | RAS+Schur true residual |
|---|---:|---:|---:|
| 4 / 0 | 0.0655766 | 0.999915 | 0.999895 |
| 4 / 1 | 0.0498797 | 0.999959 | 0.999943 |
| 8 / 0 | 0.00401523 | 0.999918 | 0.999896 |
| 8 / 1 | 0.00506316 | 0.999958 | 0.999938 |

**None of these cases converged to the target tolerance.** The baseline ILUT was numerically closer to convergence, and RAS/Schur was slower. At 4 threads without scaling, wall-clock setup+solve totals were approximately 5.73 s (ILUT), 9.37 s (RAS), and 11.62 s (RAS+Schur). RAS factor-entry proxy exceeded baseline by 8.3%; Schur exceeded baseline by 17.7%, so neither provides a fair-fill performance comparison at these settings. These are **negative numerical results**, not a demonstrated performance benefit.

[Full public-FEA benchmark logs](https://github.com/3rav/lis/actions/runs/43392825905).

The separate [Janna/Fault_639](https://sparse.tamu.edu/Janna/Fault_639) 3D penalty-contact matrix has not been benchmarked. No indefinite Lagrange-multiplier contact-Newton KKT snapshot has been tested here.

**Keep the experimental PR Draft.** Before declaring Readiness, investigate a more robust block/interior interface correction and quantify convergence, peak memory and total time at comparable fill on multiple representative contact and structural systems. Passing test/sanitizer CI does not substitute for numerical effectiveness.

## Review scope and realistic validation

This PR **depends on the ILUT RCM/equilibration PR**. The development
branch is stacked on `feature/ilut-contact-robustness`; a GitHub PR
targeting upstream `develop` will temporarily include that dependency
until the upstream RCM PR is merged.

A production contact claim is **not** supported by the algebraic 48×48
KKT regression. A meaningful comparison must use several representative
Newton tangent snapshots and balance `nnz(L+U)`, Krylov tolerance,
setup/solve time and peak memory. The companion ILUT benchmark runner
accepts an external Matrix Market matrix without checking in private
FEA data:

```sh
./test/testhalfilut --fea-matrix path/to/contact.mtx 8 8 0.005
```

The output prints actual fill and states whether the rates are
fill-matched. RAS/Schur is **not** included in the runner's two-way
baseline-versus-RCM comparison; it must be benchmarked independently
before any performance claim.

## Limits / safeguards

- **Experimental and serial-only**: ordinary real double-precision CSR, with row indices local to the single process. MPI, USER/OPERATOR and quad-precision input are explicitly rejected when the option is used; the default ILUT remains unchanged.
- The dense Schur path is used at `n_gamma <= 256`; larger interfaces use the experimental bounded sparse CSR-ILUT path. Unsafe pivots or exceeding the `4*nnz(A)` sparse fill budget cause explicit errors, not silent fallback.
- The diagonal approximation explicitly rejects interior rows with missing or weak diagonals (relative to their largest coefficient, threshold `1e-10`) and nonfinite reciprocals. The dense Schur factorization rejects pivots smaller than `1e-12` of the initial Schur maximum coefficient. These guards prevent silent invalid arithmetic; they **do not guarantee a well-conditioned or effective contact preconditioner**.
- The regression also checks that a deliberately weak interior diagonal is rejected. These tests are algebraic safeguards, not a production FEA validation.
- The nested ILUT setup may be costly in memory and time; there is currently no automatic partitioner or communication overlap.
- The adjoint/transpose RAS application is explicitly **not implemented**. Requiring an adjoint raises an error rather than returning an incorrect result.
- Domain parts own their matrices, factors and work vectors; the outer preconditioner destroys them deterministically.
- The subdomain and Schur apply workspaces are reused per preconditioner instance. Do not call **the same instance concurrently** from multiple threads without external synchronization; independent instances are separate.
- Setup OOM paths initialize the allocated subdomain records before cleanup, so a failure in the second allocation cannot traverse uninitialized pointers.
- Use flexible GMRES or GMRES for testing and inspect true residuals, iterations, setup time, factor fill and maximum RSS. For contact KKT saddle-point matrices, success is **not guaranteed**; validate across several active contact sets and scaling regimes.

The regression `testhalfilut` executes all overlaps 0/1/2 with and without Schur, checking a substantial residual decrease and explicit rejection of unsupported adjoint applies. The current code is a foundation for later **MPI-distributed** RAS plus iterative Schur and ddPQ pivot-matching work, not a claim that these have already been delivered.

Reference: https://hypre.readthedocs.io/en/latest/solvers-ilu.html
