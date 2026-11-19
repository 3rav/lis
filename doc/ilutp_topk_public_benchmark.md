# Public ILUTP limited-fill A/B benchmark

This is an independent, synthetic benchmark supporting [fork PR #101](https://github.com/3rav/lis/pull/101). No private CalculiX A/b matrices are included, and a positive synthetic result does not imply a speedup for contact mechanics.

## Pinned comparison

- Baseline: `3d027903ab079158d99bfc91f8ab6867c212dc20` (upstream `develop`).
- Candidate: `78e694a82fe15575b22df11977b7277617e9b039` (bounded top-k, #101).
- Deterministically generated sparse, nonsymmetric, strictly row-diagonally dominant Matrix Market A and manufactured b=A*x; nearby and distant couplings.
- `test/benchilutp` from both revisions: separately measures CPU setup and CPU solve, computes independent physical residual, factor nnz, column permutations and Krylov iterations.
- Sweep `-ilut_rate 1,2,4,8`; fixed one thread, FGMRES, no RCM/row scaling/equilibration, drop=0, tau=0.5/1.0; legacy ILUT as control.

## Local reproduction

Build each *pinned SHA*, not a moving branch. Standard tools: Python 3, GCC, Autotools, libtool.

```sh
git clone https://github.com/3rav/lis.git lis-base
git clone https://github.com/3rav/lis.git lis-topk
git -C lis-base checkout --detach 3d027903ab079158d99bfc91f8ab6867c212dc20
git -C lis-topk checkout --detach 78e694a82fe15575b22df11977b7277617e9b039
for tree in lis-base lis-topk; do
  (cd "$tree" && autoreconf -fiv && \
    ./configure --enable-test --disable-mpi --disable-fortran \
      --disable-f90 --disable-saamg CFLAGS='-O2 -g' && \
    make -j2 && make -C test benchilutp)
done
python3 tools/bench_ilutp_topk_public.py \
  --baseline-bin "$PWD/lis-base/test/benchilutp" \
  --candidate-bin "$PWD/lis-topk/test/benchilutp" \
  --n 4000 --far 8 --seed 20261119 \
  --rates 1,2,4,8 --warmups 1 --repeats 7 \
  --out ./ilutp-topk-results
```

The harness resides on the benchmark branch and can run from any directory if the binary paths are adjusted. For a smoke test, use `--n 1200 --repeats 3`; the GitHub Actions workflow instead runs three named profiles on pinned clean checkouts, each with seven paired A/B repetitions: `sparse-control` (n=1200, far=1), `medium-selection` (n=1200, far=8) and `large-selection` (n=4000, far=8).

## Outputs and acceptance criteria

`inputs/synthetic_A.mtx`, `inputs/synthetic_b.mtx`, and `manifest.json` record the public data, SHA-256 hashes, seed, refs, host and options. `raw.csv` retains each paired observation. `summary.csv` includes per-rate median setup and solve CPU times and `(candidate/baseline-1)*100%` changes, plus the median, minimum and maximum of **paired A/B** percentage changes (matched by repetition index). Per-run stderr/stdout and `validation.txt` are also stored.

Pairs alternate baseline/candidate order, use separate processes and one thread; warmups are excluded. For each rate and variant, **status, iterations, factor L/U counts, moved-column count, API/solver statuses must be identical**, true residuals close (abs tol 5e-12, relative tol 1e-8) and no greater than 1e-7, and all cases must converge.

**The legacy ILUT variant is a negative timing control, and the `far=1` profile probes when top-k has little work to do. Timing never gates CI.** Shared GitHub runners exhibit noisy timing and CPU frequency variation. Publish 7+ paired repetitions on an otherwise idle dedicated machine before claiming performance improvement; report neutral/high-fill results too.

## Limits

These synthetic matrices are deliberately diagonal dominant and not pivot-sensitive nonlinear FEA contact problems. This test cannot establish Newton convergence, peak memory, production speedups, or generality. It does not replace independent code review.
