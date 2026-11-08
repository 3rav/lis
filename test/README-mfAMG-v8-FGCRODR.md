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
    LIS_MFAMG_REPLAY_MAXITER=100..6000        (default 1000)

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


## Near-singular symmetric-indefinite pilot (non-gating)

The separate test/mfamg_replay_illconditioned_ci.sh constructs a PUBLIC,
deterministic 480x480 real symmetric system A=diag(T,-T), where each
240x240 block T is tridiagonal with diagonal 2+1e-5 and off-diagonal -1.
Thus its eigenvalues are +/-[2+1e-5-2*cos(j*pi/241)], j=1..240,
and its exact 2-norm condition number is about 22230.49.

Two independent Windows/UCRT64 workflow runs passed:
https://github.com/3rav/lis/actions/runs/42505597543
https://github.com/3rav/lis/actions/runs/42505602628

With restart=4, maxiter=300 and equal preconditioner settings:
- FGMRES and FGCRO-DR(k=0): 25 iterations, true_relative 9.690077339196e-6
- FGCRO-DR(k=1): 22 iterations, true_relative 9.668199352259e-6
- FGCRO-DR(k=2): 22 iterations, true_relative 8.337573750868e-6
All solvers converged with status=0. Active recycling occurred after
iterations 4, 8, 12, 16 and 20, with builds=1, updates=5, rejects=0.
Orthogonality errors were below 1e-15; final A*Z=C relation errors
were below 2e-13. k=0 gave an identical numeric fallback to FGMRES.

This follow-up promotes the case to a GATING numerical CI regression:
success/status=0, independent true-relative residual <1e-5, exact k=0
fallback, at least two actual recycled restarts, multiple harmonic
updates and valid recycle-space algebra. It does NOT assert a speedup.
Measured wall times of 0.002-0.004 seconds are not meaningful as a
performance comparison, and this synthetic SPD/negative-SPD block
system is not a substitute for a contact Jacobian.

No production code is changed. The private 31/116-spring replay remains
the essential gate before any production use or merge.


## Local contact-snapshot CSV summary (no private data in CI)

The local runner test/mfamg_replay_fgcrodr_sweep.sh now writes its case
logs under $MFAMG_V8_RESULTS_DIR (default: $MFAMG_REPLAY_DATA/v8_results).
After completing the sweep, it generates v8_summary.csv via the portable
test/mfamg_replay_v8_report.awk parser. The compact CSV includes: system
31/116, mode count, factor, solver method, recycle k, restart/maxiter,
process exit, solve return code, solver status, iteration count, both
residuals, runtime, recycle-space size/stats, true number of intra-solve
recycle *cycle activations* from trace lines, orthogonality error and AZ=C
relation error. Missing optional diagnostics are marked NA, not zero.

The standard run now includes **both** FGMRES and FGCRO-DR k=8 for the
116-spring control as well as the previous 31-spring cases. The optional
MFAMG_V8_FULL=1 includes k=0 and k=16 controls and 18-mode FGMRES.

Do not mistake process_exit=0 or solve_rc=0 for successful convergence:
the solver status and independent true_relative must be considered.
In particular, the known 31-spring baseline ran to maxiter=1000 with
true_relative around 1.832690841870e-3, not the 1e-5 tolerance.
The report deliberately keeps this nonconverged baseline as evidence.

A self-contained synthetic log-parser regression
test/mfamg_replay_v8_report_smoke.sh is run in Windows UCRT64 CI. This
test checks that a MAXITER result and its nonzero process exit are
preserved, active recycle cycles are counted independently of uses,
and missing REPLAY_RESULT is rejected. It does not run private matrices.

Logs and CSV remain on the local user's computer; do not upload private
contact snapshots or associated diagnostic logs to the public repo.


## End-to-end private sweep reporting regression

The standard sweep includes both 31-spring FGMRES and k=0 controls,
plus the 116-spring FGMRES baseline. MFAMG_V8_FULL=1 runs nine
configurations, including k=16 and offline 18-mode comparisons.
All cases are attempted even if one crashes, and the available logs
are combined into the CSV. A missing REPLAY_RESULT produces a nonzero
runner exit, while a valid MAXITER result is retained as data.

The CI-required test/mfamg_replay_v8_sweep_smoke.sh uses placeholders
and fake executables to exercise the full runner, preserve status=4,
count recycle cycles, and test handling of a crash before REPLAY_RESULT.
It never reads private data and cannot prove contact convergence.


## Automatic decision screen from private local replay

The local sweep enables recycle-space algebra validation for k>0 and
writes v8_assessment.txt through test/mfamg_replay_v8_assess.awk.
The screen requires the 31-spring FGMRES baseline (60,265 equations,
1000 iterations and true relative residual near 1.832690841870e-3),
exact k=0 parity, useful k=8 residual reduction or convergence with
valid AZ=C algebra, and preservation of convergence on the 116-spring
control (historically 592 FGMRES iterations). It compares only matching
matrix, modes, factor, restart and iteration-limit settings.
An incomplete CSV causes a clear error. A failed numerical comparison
produces FAIL_OR_INCONCLUSIVE rather than falsely suggesting success.
PASS_CONTACT_REPLAY_ONLY means the one-snapshot numerical screen passed;
PRODUCTION_GATE is always BLOCKED_NEEDS_REPEATED_TIMING_AND_NEWTON_SEQUENCE.
This is deliberately conservative: timings require repeated large runs,
and Newton-step operator changes still need validation. Do not merge.


## Extended FGMRES time-to-convergence control

The first real 31-spring 60,265-DOF replay with the OFFLINE 18-mode
factor=4 preconditioner gave FGMRES status=4 and true_relative
4.442649280794e-5 at iteration 3001 (maxiter=3000). The opt-in
FGCRO-DR(k=8) reached status=0 at iteration 2772 with true_relative
9.986595263773e-6 and a validated recycle relation (AZ=C error
7.687679465279e-9, close to the existing 1e-8 rejection threshold).

To compare time-to-the-SAME-tolerance instead of comparing converged
FGCRO-DR to a stopped FGMRES, the replay-only maximum allowable
LIS_MFAMG_REPLAY_MAXITER is raised from 3000 to 6000. This DOES NOT
change the default 1000 or the production solver, the MFAMG library,
or CalculiX. The restart memory bound remains unchanged. A CI option
smoke checks that 6000 parses and 6001 rejects before any real matrix
is opened. An updated EXE is required to use the extended limit.

The six additional negative-eigenvalue modes were built OFFLINE;
the result does not demonstrate online mode discovery. Do not merge
until production-feasible enrichment, repeated timing and changed-
operator Newton sequence validation have been demonstrated.
