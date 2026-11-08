# MFAMG v8: flexible Krylov recycling on contact snapshots

## Scope and hypothesis

The 31-spring snapshot has 60265 equations. Baseline 12-mode FGMRES(180)+MFAMG ends at 1.832690842e-3 after 1000 iterations. Six additional offline eigenmodes reduce that residual to 2.244405261e-4 with factor=4, but those modes were computed from the known matrix and cannot be assumed available during CCX.

LIS already has FGCRO-DR with an A*Z=C recycle relation and harmonic Ritz extraction. This experimental branch tests whether ONLINE Krylov recycling improves the fixed 31-spring snapshot without injecting precomputed negative eigenvectors. Flexible recycling is the important distinction from standard GCRO-DR: right-preconditioned search directions may change with the preconditioner.

## Inspiration (algorithmic ideas only; no external code copied)

- PETSc KSPHPDDM: https://petsc.org/release/manualpages/KSP/KSPHPDDM/
- PETSc KSPDGMRES adaptive deflation: https://petsc.org/main/manualpages/KSP/KSPDGMRES/
- Trilinos Belos GCRO-DR and orthogonalization: https://trilinos.github.io/belos.html
- Ginkgo flexible GMRES and restart/memory tradeoffs: https://gko-project.org/1.11/user-guide/concepts/solvers/gmres.html
- hypre BoomerAMG near-nullspace interpolation and coarsening: https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html
- Fabulous block deflated restarting: https://solverstack.gitlabpages.inria.fr/fabulous/

These motivate separate investigations of harmonic Ritz recycling and automatic local AMG mode enrichment. FGCRO-DR recycling does not directly alter MFAMG prolongation.

## Experiment controls

The new replay-only options are:

    LIS_MFAMG_REPLAY_METHOD=fgmres|fgcrodr     (default fgmres)
    LIS_MFAMG_REPLAY_RECYCLE_DIM=0..64        (default 8; only with fgcrodr)
    LIS_MFAMG_REPLAY_RESTART=2..600           (default 180)
    LIS_MFAMG_REPLAY_MAXITER=100..3000        (default 1000)

Recycle dimension must be strictly less than restart. The MFAMG library and production CCX solver are unchanged. Without new settings replay retains exactly the original FGMRES defaults.

Within-solve recycling is an additional, **default-off** experimental switch:

    LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1
    LIS_FGCRODR_INTRA_TRACE=1

Only with this switch does the newly built harmonic Ritz recycle relation become
active at the next Arnoldi cycle *within the same lis_solve call*. After each
cycle the active rank and B-column stride are refreshed, the true residual
is reprojected through C and Z, and the inherited FGMRES logic checks the
true residual at each restart. The previous cross-solve logic remains
unchanged when unset.

The replay can additionally show stable recycling diagnostics:

    LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1

This checks C^T*C approximately equals I and A*Z approximately equals C,
with a strict 1e-8 maximum error in the opt-in validation. Validation
adds at most k matrix-vector products and is never implicitly enabled.

The source-side test/test script is test/mfamg_replay_fgcrodr_sweep.sh. It uses the same A, b, and modes for all method comparisons and checks the replay method selector before loading large snapshots. Its default data directory is $HOME/mfamg_replay and it accepts MFAMG_REPLAY_DATA and MFAMG_REPLAY_EXE overrides. Set MFAMG_V8_FULL=1 for extra controls.

Review REPLAY_RESULT status, iter, true_relative and wall_seconds. A solve_rc of zero does not imply convergence when status=4. The new REPLAY_RECYCLE line reports the built recycling subspace when available.

## CI evidence and remaining validation gate

GitHub Actions Windows/UCRT64 regression succeeded:
https://github.com/3rav/lis/actions/runs/42495507072

- Legacy FGMRES(2) and FGCRODR(k=0) on a heterogeneous 240-row RHS
  match exactly: 4 iterations, true relative residual 1.117716414192e-6.
- FGCRODR with k=1 but with intra-solve recycling OFF also matches.
- With intra-solve recycling ON, an actual recycled basis was activated
  at iter=2 (k=1), and the true relative residual after 4 iterations
  was 9.570565480928e-7 (approximately 14.4% lower).
- Windows/UCRT64 compilation and testmfamgmodes also passed.

Algebraic validation passed in CI (orthogonality 4.44e-16, A*Z=C
relation 1.43e-16): https://github.com/3rav/lis/actions/runs/42495830822

The experimental Windows UCRT64 replay EXE was uploaded successfully:
https://github.com/3rav/lis/actions/runs/42496327918

**Still unverified:** the 60,265-DOF contact snapshots and any production
CalculiX integration. The small CI test is not evidence of a contact
acceleration.

## Required validation on the large snapshots

1. Compile replay on MSYS2 UCRT64 and verify the default fgmres reproduces 31-spring residual 1.832690841870e-3.
2. Compare fgcrodr recycling k=0, 8, 16 at identical restart=180 and maxiter=1000.
3. Check true residuals, status, iteration counts, actual recycle dimension, runtimes, and possible breakdowns.
4. Confirm the 116-spring control does not regress from 592 converged iterations.
5. Do NOT merge to develop or modify CalculiX without verified numerical evidence.

This work is an isolated experimental proposal, not a demonstrated acceleration yet.


## Symmetric-indefinite pilot (separate from contact evidence)

The test/mfamg_replay_indefinite_ci.sh generates a deterministic 480x480
real symmetric, strictly diagonally dominant and nonsingular matrix
with positive and negative diagonal blocks (hence indefinite), a
heterogeneous RHS and two near-nullspace modes. No private CCX data.

With restart=4 and maxiter=300, compare FGMRES, FGCRO-DR(k=0), and
intra-solve FGCRO-DR(k=1,2) on identical inputs. Check the exact k=0
fallback, a real restart, actual recycle use and C^T*C / AZ=C algebra.
Report independently computed true relative residuals and statuses.
The pilot passed in two Windows UCRT64 CI runs:
https://github.com/3rav/lis/actions/runs/42497389335
https://github.com/3rav/lis/actions/runs/42497393714

FGMRES and k=0 both converged in 7 iterations to 8.977625403445e-6.
FGCRO-DR k=1 converged in 7 iterations to 8.330184788348e-6;
k=2 in 7 iterations to 7.011277277581e-6. Both activated recycle
at iter=4. Maximum orthogonality error was 1.22e-15,
maximum AZ=C error 3.03e-16.

This change promotes the synthetic-case probe into a GATING regression:
solver convergence, independently measured true residual below 1e-5,
exact k=0 fallback, actual recycle use and valid recycle-space algebra.
It does NOT assert speedup; wall times around 0.002 s are uninformative.
Do not merge without testing the real contact snapshots.
