# CSR non-split transpose: one-worker scratch bypass benchmark

Date: 2026-11-20. This is a **microbenchmark**, not an end-to-end solver performance claim.

The [upstream submission branch](https://github.com/3rav/lis/tree/upstream/perf-csr-transpose-serial-20261120) skips the OpenMP scratch allocation, dense zeroing and parallel region when `omp_get_max_threads()==1` for **non-split CSR**. It uses the existing row-order serial accumulation and leaves all other cases unchanged.

## Matched A/B measurements

The optional diagnostic branch contains `tools/bench_csr_transpose_oneworker.c` and the paired GitHub Actions workflow. Each job builds upstream develop at commit `096478e4483f7f60e46b46d1a9faf1262dcdbe72` and the optimized fork branch using the same Ubuntu runner, compiler flags (`-O2`) and OpenMP configuration. For each case, the test checks the true transpose result numerically before measuring, runs repeated calls in batches, alternates baseline and optimized measurements, and checks the numerical checksums across versions. The workflow gates **correctness only**, not performance thresholds.

Median nanoseconds per transposed SpMV:

| Size | A/B run #44714040939 baseline | Optimized | Speedup | A/B run #44714699995 baseline | Optimized | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| 256 | 1668.8 | 784.3 | 2.13x | 1670.3 | 862.9 | 1.94x |
| 1024 | 5727.5 | 4248.8 | 1.35x | 5679.4 | 4910.8 | 1.16x |
| 4096 | 21491.2 | 17345.6 | 1.24x | 21761.0 | 17998.9 | 1.21x |

Both runs validated 84 batches with matching numerical checksums.

- [First full A/B run and raw CSV](https://github.com/3rav/lis/actions/runs/44714040939)
- [Independent repeat and raw CSV](https://github.com/3rav/lis/actions/runs/44714699995)

**Limitations:** single OpenMP worker, a synthetic banded CSR matrix with seven nonzero entries per row, short absolute runtimes, shared CI runner, isolated `lis_matvech_csr` calls. These results support reducing one-worker transpose overhead but **do not demonstrate comparable end-to-end LIS, MPI, multi-thread, matrix assembly or nonlinear contact speedups**. Some differences between the two runs indicate scheduler/noise sensitivity. No changes to production CSR code were made by this benchmark branch beyond its base.

The production-ready CSR PR intentionally excludes the benchmark harness and its heavy two-build CI workflow. They remain here for reproducibility.
