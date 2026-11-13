# Experimental sparse interface Schur: standalone proof of concept

**Status: research-only, not integrated into LIS.** This branch is stacked on the experimental serial RAS/Schur branch (fork PR #63 / upstream #225) and deliberately does not modify that implementation.

## Scope

`test/experimental_sparse_schur.c` implements a standalone, C99, CSR-only approximation of

```text
S_tilde = A_GammaGamma - A_GammaI diag(A_II)^(-1) A_IG
```

using sparse row accumulation with marker arrays rather than a dense `ngamma × ngamma` allocation. It exercises an **illustrative 384-DOF interface** (larger than the current 256-DOF dense limit), checks the sparse operator action independently using the original block rows, and solves a deliberately diagonally dominant synthetic example with Jacobi-preconditioned Richardson.

This is **not** a robust general Schur solver, an indefinite KKT preconditioner, or evidence of scalability to real 3D contact. The row accumulator takes `O(ngamma)` transient memory and the assembled operator takes `O(nnz(S_tilde))` memory. The actual Schur fill can still become dense and the nested accumulation cost can be significant.

Run it as:

```sh
cc -std=c99 -O2 -Wall -Wextra -Werror -pedantic \
  test/experimental_sparse_schur.c -lm -o sparse-schur-poc
./sparse-schur-poc
```

A dedicated workflow also compiles it under ASan and UBSan.

## Gate for any LIS integration

Before any solver-level integration, implement validated sparse storage with overflow checks and bounded fill, a robust iterative interface solve (e.g. FGMRES/ILUT with finite-precision safeguards), termination/failure propagation, independent normal and adjoint behavior (or explicit rejection), and tests on highly indefinite/zero-diagonal KKT systems and real FEA contact matrices. Sparse storage alone does not fix the failure of `diag(A_II)^(-1)` on multiplier interiors. Evaluate interface partition quality and factorization costs separately.

## OpenMP assessment for serial RAS

The current RAS apply path loops over independent subdomains. Each subdomain has its own `p->solver`, `p->b`, `p->x` and writes to disjoint owner-core entries in output `X`. **This suggests a potential per-domain OpenMP parallel region**, followed by a barrier before the global Schur residual/correction. The following must be verified before adding parallelism:

1. Thread safety and reentrancy of `lis_psolve`, local solver/preconditioner state, and error reporting; avoid nested OpenMP oversubscription if local ILUT uses threads.
2. Deterministic error handling: workers must report individual error codes and must not write to shared error variables without synchronization.
3. Read-only shared `B` and `s->parts`, disjoint writes to `X`, and privately owned per-domain work vectors; clean shutdown of workers on errors.
4. Factorization setup is currently serial; parallelizing only apply does not accelerate setup. Measure setup, apply, memory and complete FGMRES solve separately.
5. Test non-divisible partitions, overlap 0/1/2, weak diagonals, FP16, OOM, ASan/UBSan/TSan where available, and numerical results across 1/2/4/8 threads.

**No nearly-linear OpenMP speedup or thread safety is claimed.** Do not modify upstream PR #225 until the architecture and numerical tests justify it.

## Branch relationship

- Base: `feature/ilut-ras-schur` in `3rav/lis` (fork PR #63).
- Experimental standalone additions only: C99 prototype, documentation and CI.
- Upstream PR #225 remains Draft.
