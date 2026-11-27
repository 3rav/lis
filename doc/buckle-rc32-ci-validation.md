# RC32 recovery: verified port results

Verified head: [532d2b604166e6eceed2891e86c3277f3455350d](https://github.com/3rav/lis/commit/532d2b604166e6eceed2891e86c3277f3455350d)
Upstream baseline: 6a4ae9c11bb03b028a6477455cb4f9f0a790132d.
[GitHub Actions run 45992236387](https://github.com/3rav/lis/actions/runs/45992236387): completed, success, all eight jobs.

| Configuration | Standard make check PASS | SKIP | RC32 driver cases PASS |
| --- | ---: | ---: | ---: |
| Serial | 143 | 2 | 62 |
| OpenMP (1 and 4 threads for RC32) | 143 | 2 | 124 |
| Serial ASan/UBSan | 143 | 2 | 62 |
| Half | 149 | 2 | not run |
| Wide index | 143 | 2 | not run |
| Complex | 124 | 21 | not run |
| Complex + wide index | 124 | 21 | not run |
| MPI single rank | targeted smoke tests passed | — | not run |
| Total counted | 969 | 52 | 248 |

Zero failures/errors in the seven full legacy test suites. The dedicated driver checks require both physical residual acceptance and explicit symbolic-LU / refresh selection markers. Positive NEV/inner-option test markers were observed in all six real configurations including the single-rank MPI smoke tests. Complex builds test legacy compatibility and documented real-only rejection, not experimental complex LU.

The public include/lis.h is unchanged relative to baseline. The five preconditioner files remain byte-identical to the initial recovered port commit 50757b160d3bd77fd48ab1c7a97f740ea8fb523e (Git blob comparison):
- src/precon/lis_precon_ilutp.inc: a32fb3d666c0b2a9f5438e53c5b284b5ec8b90ae
- src/precon/lis_ilutp_refactor.inc: 853826573ca203c12baa09ec5ad38d0b9489f3a5
- src/precon/lis_ilut_min_degree.inc: 15cf723f37c4cdef802b39df63673a2e8320a2cd
- src/precon/lis_precon_ilut.c: e96179e245fcb247fc8bd9d65f445e0585f119dd
- src/precon/lis_precon.c: fc3a804a61a0a8492359974f025ec4ab7ffa3f7d

Original consolidated patch: [archived unchanged](https://github.com/3rav/lis/blob/archive/buckle-rc32-original/recovered/RC32/LIS_2.1.25_RC32_consolidated.patch), Git blob 5bdd49e4a29243b9e873cc8e99fdb3adc9c37823.

Disposition: suitable for experimental draft review, not production acceptance. The ABI-changing near-nullspace prototype is excluded from the port but preserved in the archive. Original RC32 parallel numeric factorization is not exercised by the recommended serial-numeric tests; no TSan or production-scale OOM guarantee. Real CCX matrices and reference vectors are missing, so historical speed/RSS results have not been independently reproduced. The experimental environment options remain opt-in; legacy compatibility runs do not imply feature support for MPI, wide-index, half or complex.
