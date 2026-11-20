# GCRO-DR orthogonalization: double MGS vs adaptive DGKS vs C-only pass

Date: 2026-11-20. **Diagnostic-only experiment. Not production-ready.**

## Original code and hypothesis

GCRO-DR performs two Modified Gram–Schmidt passes on both the
recycle-image basis C and the current Arnoldi basis V. The extra
projections and global reductions may cost more than necessary
when the first projection is stable. However, simply removing
second-pass orthogonalization can be numerically unsafe.

This work is on a separate experimental branch atop the diagnostic
dynamic-operator PR #109; it does **not** change upstream LIS nor
fork production PRs #107/#108.

An explicitly compiled opt-in macro,
LIS_EXPERIMENT_GCRODR_ADAPTIVE_MGS, enables these run-time
diagnostic modes through LIS_EXPERIMENT_GCRODR_ADAPTIVE_MGS:

- full: original two MGS passes on C and V (the default, unchanged);
- dgks: a second pass on either basis only if the sum of first-pass
  squared projection coefficients exceeds 0.945 of squared remaining
  vector norm, or if the norm is nonfinite. This is a DGKS-style
  heuristic, not a proof of orthogonality;
- one_c: always project once against C, but retain the original
  unconditional **two passes against V**.

Unlike one_c, dgks needs an extra norm reduction after the C
projection. A dedicated diagnostic counter records actual second
pass decisions. No runtime environment switches are compiled into
normal LIS builds.

## Experiment and validation

- The same six consecutive numeric-value changes to a persistent,
  deterministic nonsymmetric CSR system as in PR #109, with the
  sparsity pattern held fixed.
- Jacobi/ILUTP preconditioners; n=784 and 1600; restart=12/24;
  k=4/8; Ritz stride=1/4; four trials per method, first as warmup.
- Full baseline, DGKS, one_c, then repeated full baseline on the
  **same runner** to account for drift.
- All 8448 Release solves (4 modes x 2112) converged with
  independently recomputed true residual <=1e-8.
- Release and AddressSanitizer/UndefinedBehaviorSanitizer validation
  passed for all three orthogonalization modes.
- Four existing GCRO-DR regressions ran with full/dgks/one_c.
  The strengthened stride regression independently checks that
  the **stored recycle subspace** satisfies C^H C approximately I,
  and A*ZU approximately C, using 1e-8 absolute tolerances.
  All tests passed.
- The exercise is serial (one OpenMP worker). **No MPI adaptive
  orthogonalization test or scaling claim is made.**

[Exact-source workflow and artifacts](https://github.com/3rav/lis/actions/runs/44696566340).

## Measured result

Across 32 positive-k GCRO-DR configurations (two preconditioners,
sizes, restart lengths, recycle dimensions and Ritz strides):

| Metric | Adaptive DGKS | Single pass against C |
|---|---:|---:|
| Median speedup relative to full two-pass MGS | **1.033x** | **1.041x** |
| Configurations faster in matched experiment | **32/32** | **32/32** |
| Median change in Krylov iteration count | **0** | **0** |

From instrumented DGKS runs:

| Basis | Second pass performed | Second pass skipped |
|---|---:|---:|
| C (recycle image) | 0 | 19,780 |
| V (Arnoldi) | 19,780 | 1,792 |

For the one_c variant, the C second pass was skipped in all
19,780 eligible cases. The V path retains two passes; the
diagnostic counter is not collected in this mode.

These are **short synthetic shared-runner CPU measurements**,
not high-confidence production performance estimates. The
FGMRES baseline remained considerably faster than either
recycling algorithm in prior tests.

## Interpretation and next steps

1. In these workloads, the **first projection on C was sufficient
   empirically**, and both conditional approaches yielded modest
   CPU reductions. The second V pass was much more often warranted
   under the projection-norm criterion.
2. The median gain of only 3–4% is considerably smaller than the
   benefit of reducing Ritz UPDATE frequency with cheap Jacobi.
   This favors prioritizing recycle selection/transport overhead
   or representative application matrices over immediately making
   the orthogonalization change a production default.
3. A single C projection is potentially unsafe for nearly
   dependent directions, highly nonnormal/ill-conditioned matrices,
   or low precision. Before a production PR, require targeted
   numerical conditioning stress, more accurate orthogonality
   tests, MPI branch-consistency tests, and real changed-operator
   A/b snapshots from a contact FEA sequence.
4. No evidence yet justifies enabling one_c or dgks by default
   or claiming overall superiority over FGMRES.

The experiment intentionally does not replace the previous
diagnostic PR #109 nor alter fork PR #107/#108.
