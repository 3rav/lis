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

## Synthetic contact-KKT regression

`testhalfilut` now assembles an indefinite **48×48 stabilized saddle-point system** with a 36×36 spring-stiffness block, twelve displacement-difference constraints, and a negative 0.1 multiplier diagonal. It compares the legacy ILUT control with the opt-in variant, requiring convergence to a **true relative residual below 1e-7** and printing Krylov iterations and residuals. This is an algebraic contact-style regression, **not** a production nonlinear FEA benchmark and not evidence of a speedup.

## Dependency update (2026-11-13)

This branch is synchronized with the latest ILUT robustness changes from fork PR #62 (upstream #224), including the high-unit-contrast KKT test and the English user-guide option table. The fork PR #63 remains stacked cleanly on #62, with its diff restricted to RAS/Schur changes. The upstream #225 comparison against `develop` necessarily includes #224 until that dependency merges.

An entirely **separate, non-integrated** sparse CSR Schur feasibility experiment lives in fork Draft PR #64. It checks a synthetic interface larger than 256 DOF and documents OpenMP risks, but it is not a fix for the production limitations of the current dense Schur implementation.

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
