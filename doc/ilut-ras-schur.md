# Experimental serial RAS-ILUT and approximate interface Schur

LIS ILUT now has a **separate opt-in** serial restricted additive Schwarz (RAS) variant, inspired by the decomposition architecture of hypre. It does not depend on hypre or use hypre code.

## Configuration

```
-i fgmres -p ilut -scale none -ilut_ras_domains 4 \\
-ilut_ras_overlap 1 -ilut_schur 1 -ilu_pivot_tol 1e-12
```

- `-ilut_ras_overlap 0|1|2` turns on RAS (default: off). Overlap is the number of sparsity-graph edges added beyond each core partition.
- `-ilut_ras_domains 2..16` selects the count of deterministic, contiguous partitions (default: 2). Each subdomain uses the normal LIS ILUT factorization on its extracted local matrix, including the outer `-precon_half` policy.
- `-ilut_schur 0|1` adds an optional **diagonal-approximate interface Schur correction** (default: off). The interface is rows coupled across the original nonoverlapping partitions.

RAS applies the local factor solve on each overlapping subdomain but injects corrections into that subdomain's **nonoverlapping core** only. For the Schur variant, LIS computes

\[\widetilde S=A_{\Gamma\Gamma}-A_{\Gamma I}\operatorname{diag}(A_{II})^{-1}A_{I\Gamma}\]

and uses a small pivoted dense LU factorization to correct the interface residual. It then applies diagonal-approximate interior back-substitution. This is a **one-step approximate Schur correction**, not hypre's distributed GMRES-on-Schur algorithm and not a general-purpose exact Schur complement.

## Synthetic contact-KKT regression

`testhalfilut` now assembles an indefinite **48×48 stabilized saddle-point system** with a 36×36 spring-stiffness block, twelve displacement-difference constraints, and a negative 0.1 multiplier diagonal. It compares the legacy ILUT control with the opt-in variant, requiring convergence to a **true relative residual below 1e-7** and printing Krylov iterations and residuals. This is an algebraic contact-style regression, **not** a production nonlinear FEA benchmark and not evidence of a speedup.

## Limits / safeguards

- **Experimental and serial-only**: ordinary real double-precision CSR, with row indices local to the single process. MPI, USER/OPERATOR and quad-precision input are explicitly rejected when the option is used; the default ILUT remains unchanged.
- Maximum interface size for the dense Schur option is 256 DOF. A larger interface or singular approximate Schur leads to an explicit error rather than an undocumented silent fallback.
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
