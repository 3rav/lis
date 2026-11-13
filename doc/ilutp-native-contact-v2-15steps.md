# Native LIS on a 15-step CCX Newton snapshot series (2026-11-13)

These results are from **original Matrix Market A,b pairs** supplied as a
private user-uploaded 7z archive `lis_snapshots_v2(1).7z`.
The raw 15 matrices, RHS vectors and DOF mappings are **not** in GitHub.
The archive's naming prefix was `pivot1e9_test2`, and the
matrices had 60,265 unknowns and 3.62–3.65 million nonzeros each. This
archive does **not** reliably identify a 187-spring or 256-spring model.

Implementation: exclusively native LIS ILUTP with **full fill**,
RCM symmetric reordering, row-max scaling, `-ilutp_threshold 0.1`,
`-ilutp_full_limit 60000000`, FGMRES(40), solver
rtol=1e-6. No SuperLU code, API, plugin or linking was used. Each RHS
was the original exported vector; there was no manufactured b or diagonal
shift to the linear system. Every accepted solution was independently
checked against the **original** `||Ax-b||_2/||b||_2`,
acceptance threshold 1e-5. Original input files were SHA-256 checked
before and after every solve.

Each row below is a separate factorization (even when a byte-identical
A appears in multiple steps). Factors were not cached between cases.
All 15 solves converged in **one FGMRES iteration**. Steps 1–15 had
a 45-second per-case timeout except step 3, which first timed out at
45 seconds and then passed in a separate 150-second run with exactly
the same numerical options.

| Step | Factor CPU s | L+U nnz | Peak RSS KiB | True relative residual |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 9.59 | 23533512 | 431836 | 4.631e-15 |
| 2 | 38.30 | 39107753 | 638232 | 2.028e-14 |
| 3 | 65.04 | 51108311 | 799024 | 1.664e-12 |
| 4 | 9.11 | 23533512 | 431868 | 4.631e-15 |
| 5 | 28.61 | 33338232 | 561388 | 3.816e-14 |
| 6 | 9.79 | 23533512 | 431704 | 4.631e-15 |
| 7 | 12.99 | 26710622 | 473212 | 4.943e-14 |
| 8 | 12.13 | 26021990 | 464356 | 5.722e-14 |
| 9 | 10.55 | 24611888 | 444312 | 2.837e-12 |
| 10 | 10.32 | 24756200 | 446168 | 5.145e-14 |
| 11 | 11.31 | 25737599 | 459092 | 5.952e-13 |
| 12 | 9.50 | 23533512 | 431676 | 4.631e-15 |
| 13 | 10.88 | 24986252 | 449248 | 8.641e-14 |
| 14 | 10.29 | 24786986 | 446804 | 4.747e-13 |
| 15 | 10.10 | 24440936 | 442352 | 1.513e-13 |

Totals / boundaries:
- **15/15 passes** at original residual tolerance 1e-5, max actual
  residual 2.837e-12.
- **12 distinct matrices A**: the A files in steps 1, 4, 6 and 12
  have identical SHA-256 hashes. Reusing factors for identical A
  is a potential improvement not covered by the current run.
- Highest fill **51,108,311** off-diagonal factor nonzeros,
  peak RSS **799,024 KiB (~780 MiB)** on the step-3 process,
  step-3 factorization **65.04 CPU seconds**.
- Sum of reported factor CPU times: **258.5 seconds**. This is not
  the wall time of a CCX nonlinear analysis and excludes orchestration.
- Initial full-fill baseline at RCM=0 and forced `tau=1`
  failed with a 20-million fill cap and timed out at 150 s
  with a 60-million cap on step 1. The revised **RCM=1 / tau=0.1**
  combination was successful on the same original A,b.
  This comparison changes two parameters and is **not** a controlled
  attribution of improvement to RCM or tau individually.

The current branch only changes LIS and its tests. It does **not**
automatically enable fallback in CCX and does not prove end-to-end
Newton convergence, factor-cache correctness, behavior under MPI,
or across other FEA models. The full-fill mode must remain opt-in.

To reproduce with the original private files, build native LIS and use:
```sh
python3 tools/run_native_ilutp_contact.py \
  --bench ./test/benchilutp \
  --case newton01 /path/pivot1e9_test2_00000000000000000001_newton_A.mtx \
         /path/pivot1e9_test2_00000000000000000001_newton_b.mtx \
  --full-limit 60000000 --full-tau 0.1 \
  --rcm 1 --row-scale 1 --restart 40 --maxiter 80 \
  --accept-residual 1e-5 --rtol 1e-6 \
  --max-rss-mib 3072 --max-vmem-mib 3584 --timeout 150 \
  --report native-ccx-audit.json --csv native-ccx-audit.csv
```

The GitHub CI checks smaller synthetic and numerical regressions, but
cannot rerun the private CCX series without the original user archive.
