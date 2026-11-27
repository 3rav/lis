# BUCKLE recovery ledger (2026-11-27)

Source: user-provided text export of the earlier BUCKLE work (178 messages). The export contains historical reports and artifact names, not the original patch bytes or runnable CCX matrix. Historical measurements below are reported evidence, not reproduced in this audit. No private matrix or conversation transcript is published.

## Existing LIS: historical RC31/RC32 patches

| ID | Historical change | Original decision | Recovery status |
| --- | --- | --- | --- |
| 01 | Block GLI requested NEV convergence | keep | Fresh reconstruction for transformed GLI on fix/buckle-block-gli-recovered, tests added, not run |
| 02 | Block GLI inner linear options | keep | Fresh reconstruction on same branch, tests added, not run |
| 03 | Selective shifted CSR materialization | keep | Existing upstream scalar rebuild helper for ILU/ILUT/ILUC reused by reconstructed block path; MFAMG-specific path not verified |
| 04 | Inner near-nullspace propagation API | prototype, ABI change | Existing solver coarse-basis API is available; eigensolver propagation not recovered |
| 05 | Continue after MAXITER without clearing iterate | separate prototype | Not recovered; audit retry policy and ensure physical residual certification |
| 06 | Singleton absorption | reject | Do not import rejected experiment |
| 07–08 | Experimental ILUTP dropping variants | reject | Do not import rejected experiments |
| 09 | Symbolic complete-fill LU | keep | Original source missing; upstream native ILUTP uses a serial dynamic-column-pivot numerical path, exact historical specialization not established |
| 10 | Paired LU updates | keep | Source missing |
| 11 | Grouped LU8 updates | keep | Source missing |
| 12 | ILUTP serial apply with OpenMP-enabled build | keep | Source missing; verify serial pivot factor/solve dispatch separately from legacy OpenMP ILUT |
| 13 | Grouped LU16 | reject | Do not import rejected variant |
| 14 | Parallel numerical LU | opt-in, off | Source missing; race audit required |
| 15 | DOT4 triangular solve | keep | Source missing; benchmark separately |
| 16 | Integrated Minimum Degree | keep | Source missing; RC1 contains a separate prototype ordering, not historical LU integration |
| 17 | Symbolic structure cache | keep | Existing frozen-pattern ILUT (#154) is related but not equivalent to complete-fill LU symbolic caching |
| 18 | Direct CSR scatter during numerical update | keep | Source missing |
| 19 | Integrated numerical refactorization in sequence and eigensolver | keep | Source missing; existing sequence/session facilities (#154, #262) are related but not equivalent |
| 20 | Shared numerical LU kernel, cache without OpenMP (RC32) | keep | Source missing |
| 21 | Direct trial before FGMRES | reject | Historical benchmark found regression |
| 22 | Group8 SIMD | reject tested variant | Historical benchmark found no benefit |

## New LDLT: historical stages 35–54 (patches 23–43)

| ID | Historical change | Original decision / disposition |
| --- | --- | --- |
| 23 | Left-looking LDLT | develop |
| 24 | Single factor LDLT | retain |
| 25 | Grouped forward solve | retain |
| 26 | Grouped numeric LDLT | tested |
| 27 | LDLT symbolic cache | tested |
| 28 | Grouped backward LDLT | tested |
| 29 | Supernodal LDLT, native C | tested |
| 30 | DGEMM/OpenBLAS path | tested historically, excluded by user requirement |
| 31 | Native register kernel | tested, no BLAS |
| 32 | Parallel supernodal LDLT | tested historically; race audit required |
| 33 | Parallel backward LDLT | tested historically |
| 34 | Forward gather | reject: slower on CCX |
| 35 | Parallel vector mapping | reject: no repeatable apply gain |
| 36 | Packed LDLT in-place | tested |
| 37 | Index pool compression | tested; historical ABI change |
| 38 | Supernode dependency graph replacing temporary L | tested; refactorization benefit |
| 39 | Supernode ready queue | tested experimentally |
| 40 | Critical-path queue | optional, inconclusive overall timing |
| 41 | Supernode column tasks | stage 52 recommended |
| 42 | Shared stripe pack | optional, inconclusive |
| 43 | Adjacent 2x2 pivot LDLT | experimental, tested on small cases; no 2x2 used on CCX |

PR #269 contains standalone RC1 stages 55–70 with related features. It is NOT the historical code for patches 23–43 and does not establish numerical equivalence, CCX speed or native preconditioner integration. Preserve no-BLAS requirement. OpenMP and persistent cache remain experimental; do not enable by default.

## Reported historical evidence

- 152,592 DOF, six buckling modes.
- RC31, stage 31: 12.35 s, 879 MiB, 6/6, max original residual 4.57e-11; standard regression 125 PASS, 1 SKIP.
- RC32, stage 32: 11.85 s control, 879 MiB, 6/6.
- Stage 52: 9.70 s in the stage 54 comparison, 635 MiB, 6/6, max residual 3.032e-11. Level scheduler plus column tasks; no BLAS.
- Stage 54: 11.29 s, 698 MiB, 6/6, max residual 3.703e-11. Patch 43 regression 232/232 reported; 2x2 pivots were not required by the real CCX case.

## Missing artifacts (exact names where reported)

- LIS_BUCKLE_stage10_patches.zip
- LIS_BUCKLE_RC31.zip
- LIS_BUCKLE_RC32.zip
- Stage 52 source package, patch 41
- Stage 54 source package, patch 43
- Actual CCX K/KG and RHS or input model for independent benchmark reproduction

The export strips attachment links and does not contain those archives. The source code must be recovered from the actual archives or reconstructed and independently tested; never label a reconstruction as a recovered original. No benchmark on the actual CCX matrix was run in this audit.
