# Native LIS ILUTP validation on CCX contact Newton matrices

**Status: mixed outcome; keep ILUTP experimental and opt-in.** These results are from an actual serial LIS executable built by GitHub Actions, not SciPy or SuperLU.

## Inputs and methodology

- The private CCX 2.24 snapshot archive `lis_snapshots_v2.7z` contains real physical Newton matrices and RHS vectors. Only statistics are published here, **not the underlying snapshots**.
- Snapshot `03` (256 contact springs): 60,265 x 60,265 and 3,639,839 stored entries.
- Snapshot `05` (187 contact springs): 60,265 x 60,265 and 3,634,871 stored entries.
- Binary: `test/benchilutp`, compiled from this PR in [GitHub Actions run 43480090281](https://github.com/3rav/lis/actions/runs/43480090281).
- Fixed settings: serial CSR, FGMRES, one thread, restart 120, maxiter 240, tol 1e-5, no LIS global scaling, RCM, row-max scaling, ILUT drop 0.005, rate 12, pivot regularization 1e-9.
- Actual residual: `||A*x-b||_2 / ||b||_2` on the original physical A and b. Reported times are single-thread CPU seconds from the benchmark tool.
- Moved columns count the nonidentity entries of the **final** column permutation (not the number of swaps).

## 187 contact springs — ILUTP helps

| Variant | LIS status | Iter | Physical relres | Setup CPU s | Solve CPU s | L nnz | U nnz | Moved cols |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Legacy ILUT | MAXITER | 241 | 1.01477492352054 | 4.8241 | 24.2147 | 5,697,332 | 5,498,058 | 0 |
| **ILUTP τ=0.5** | **SUCCESS** | **76** | **9.5305e-6** | **5.0921** | **2.9702** | 5,415,236 | 5,087,917 | 40 |
| **ILUTP τ=1.0** | **SUCCESS** | **76** | **9.4012e-6** | **5.2295** | **3.0036** | 5,419,048 | 5,107,171 | 212 |

With tau=0.5, the preconditioner changed only 40 final permutation positions, yet it reached the target residual with approximately the same factor setup cost and lower fill than legacy ILUT. LIS may count maxiter+1 at the iteration limit.

## 256 contact springs — ILUTP does not fix breakdown

| Variant | LIS status | Iter | Physical relres | Setup CPU s | Solve CPU s | L nnz | U nnz | Moved cols |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Legacy ILUT | BREAKDOWN | 1 | 1 | 14.1174 | 14.2416 | 7,287,173 | 10,132,181 | 0 |
| ILUTP τ=0.5 | BREAKDOWN | 2 | 1 | 34.6607 | 33.3866 | 13,273,971 | 10,355,890 | 20,547 |
| ILUTP τ=1.0 | BREAKDOWN | 2 | 1 | 35.8344 | 35.7726 | 13,410,461 | 10,584,738 | 32,897 |

The true residual remained approximately 1 in all three cases. ILUTP caused more fill and more than doubled factor setup cost, while still breaking down. This is a *negative result* and must be visible alongside the positive 187-spring result.

Both input matrices have nonzero diagonal entries. Minimum `abs(A[i,i]) / max_j abs(A[i,j])`: 0.0815 (256 springs), 0.5858 (187 springs). The pathological pivots therefore arise during elimination rather than from simply zero diagonal entries of the original matrices.

## Reproduce

Build this PR with serial real-double LIS and test targets. Keep the original CCX Matrix Market A and b pair and run:

    OMP_NUM_THREADS=1 ./test/benchilutp \
      /path/to/pivot1e9_test2_00000000000000000005_newton_A.mtx \
      /path/to/pivot1e9_test2_00000000000000000005_newton_b.mtx \
      --maxiter 240 --restart 120 --rtol 1e-5 \
      --drop 0.005 --rate 12 --pivot-tol 1e-9 \
      --rcm 1 --row-scale 1 --csv ccx_187.csv

Change snapshot 05 to 03 for the 256-spring experiment. The CLI reports setup/solve CPU time, original true residual, L/U fill and the final permutation moved count. The CLI was separately smoke-tested against a deterministic 100x100 Matrix Market input in CI.

## Recommendation

**Keep upstream PR Draft and ILUTP opt-in**; the 187-spring success is strong evidence for value on some contact systems, but the harder 256-spring case is still a breakdown with significant extra fill. Before promotion: diagnose early FGMRES/preconditioner breakdown, measure pivot growth and memory, and benchmark against matching-based ordering (such as MC64) and SuperLU row pivoting.

These are direct linear-solver snapshots, not a nonlinear CCX increment acceptance test. No claim is made that the overall CCX nonlinear analysis now converges.
