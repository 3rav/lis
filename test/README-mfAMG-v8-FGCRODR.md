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


## Extended real contact baseline: 6000 FGMRES iterations (2026-11-08)

The user ran the Windows/UCRT64 v8 replay built with the 6000-iteration
bound on the local, private **31-spring, 60,265-DOF** snapshot, using
the same OFFLINE 18 modes (12 baseline + 6 negative-eigenvalue modes),
MFAMG fine-aggregation factor=4, smooth_steps=3, projected-relative
criterion and restart=180 as in the prior converged k=8 run.

| Method | maxiter | Solver status | iterations | independently measured true-relative residual | seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| FGMRES | 3000 | 4 (MAXITER) | 3001 | 4.442649280794e-5 | 207.988 |
| **FGMRES (extended)** | **6000** | **4 (MAXITER)** | **6001** | **3.833810842564e-5** | **399.776** |
| FGCRO-DR k=8 (same preconditioner) | 3000 | **0 (converged)** | **2772** | **9.986595263773e-6** | **226.086** |

After an additional 3000 FGMRES iterations the residual improved
only 13.70%, and remained 3.834 times the convergence tolerance of
1e-5. The measured FGMRES run took 173.690 seconds longer than the
converged FGCRO-DR run yet did NOT achieve the same accuracy.

**Scientific claim:** in these single local runs, FGCRO-DR reached the
requested tolerance within 3000 iterations whereas FGMRES did not
reach it even within 6000 iterations. This is strong evidence of
improved *convergence* on this one fixed real contact snapshot and
preconditioner combination, not a statistically validated speedup
at matched target residual: the FGMRES time-to-tolerance is unknown.

**Remaining limits:** the extra six near-negative modes were selected
OFFLINE (not demonstrated in online production CalculiX), wall times
have no repeated-run uncertainty estimates, the changed-operator
Newton sequence is untested, and the FGCRO-DR final AZ=C validation
error of 7.687679465279e-9 is relatively close to the 1e-8 threshold.
Neither the private contact files nor raw logs are published.
Keep PR #47 Draft; no production changes or merge.


## Factorial isolation: 12/18 OFFLINE modes × MFAMG factor=2/4 (2026-11-08)

The user ran the missing two combinations on a private 31-spring,
60,265-DOF real contact snapshot, with restart=180, maxiter=1000 and
identical solver/preconditioner options. Results below integrate the
earlier real-snapshot runs, and are **single-run measurements**:

| Modes | MFAMG fine factor | FGMRES true_relative | FGMRES seconds | FGCRO-DR(k=8) true_relative | FGCRO-DR seconds |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 2 | 1.832690841867e-3 | 33.995 | 4.518330164961e-4 | 45.963 |
| 12 | 4 | 1.924037018838e-3 | 29.754 | 4.852053168278e-4 | 41.307 |
| 18 | 2 | 2.361246837060e-4 | 96.903 | 1.778285835728e-4 | 108.061 |
| 18 | 4 | 2.275074776069e-4 | 70.134 | 1.677387800578e-4 | 82.006 |

All eight solver runs stopped at status=4, iter=1001; do not
misinterpret the lower residual as convergence to tol=1e-5.
Every FGCRO-DR run with k=8 had successful validation; the two
new runs had AZ=C relation errors 3.133714268681e-11 (12/f4)
and 1.051279499587e-10 (18/f2), and orthogonality errors below
2e-14.

**Isolated mode effect:** at factor=2, 12→18 offline modes reduce
FGMRES residual 7.76× and FGCRO-DR(k=8) residual 2.54×. At
factor=4, the factors are 8.46× and 2.89×, respectively. The extra
six offline negative-eigenvalue modes dominate the improvement.

**Isolated aggregation effect:** at 12 modes, factor=2→4 actually
increases residual by 4.98% (FGMRES) and 7.39% (k=8); at 18 modes,
it reduces residual by only 3.65% and 5.67%, respectively.
For 18 modes, factor=4 single-run wall times were about 27.6%
lower for FGMRES and 24.1% lower for k=8, but these timings are
not repeated benchmark measurements.

**Completed diagnostic:** rerun 18 modes / factor=2 / FGCRO-DR(k=8)
at maxiter=3000, to test whether factor=4 was needed for
convergence (both factor=2 and factor=4 now converge at this limit;
see the result below). This does not replace
the central need to create the six extra modes online: the 18-mode
snapshot remains an OFFLINE diagnostic, not a demonstrated production
recipe. Keep PR Draft / NO MERGE, preserve v7 and CalculiX.


## Factor=2 versus factor=4 at actual 18-mode convergence (2026-11-08)

The user completed the 31-spring/60,265-DOF real snapshot with
18 OFFLINE modes, FGCRO-DR(k=8), restart=180, maxiter=3000,
and the same MFAMG settings except fine-aggregation factor:

| Metric | factor=2 | factor=4 |
| --- | ---: | ---: |
| Solver status | 0 (converged) | 0 (converged) |
| Iterations | 2634 | 2772 |
| true_relative | 9.868651998311e-6 | 9.986595263773e-6 |
| wall_seconds | 279.149 | 226.086 |
| Harmonic-Ritz updates | 14 | 15 |
| Recycle rejections | 0 | 0 |
| orth_error | 2.198241588758e-14 | 2.575717417130e-14 |
| relation_error (AZ=C) | 7.574320088498e-9 | 7.687679465279e-9 |

**Interpretation:** factor=4 is *not necessary for convergence*.
With factor=2 the method converges in 138 fewer iterations
(4.98% fewer than factor=4), but spends 53.063 seconds more
(23.47% longer relative to factor=4). This demonstrates the
difference between iteration count and full preconditioned cost:
factor=4 is the faster one of these two *single* converged
FGCRO-DR measurements. Both AZ=C errors pass the 1e-8 cutoff,
though they are close enough to merit longer-run stability checks.

The extra six modes were selected OFFLINE and remain unproven
in an automated CalculiX workflow. The remaining research gates
are (1) deterministic online negative-mode discovery/reconstruction,
(2) repeat-run timing variability, and (3) changed-operator Newton
sequence tests. No production integration or merge is justified.

## Read-only spectral audit of OFFLINE added modes

Before building an online negative-mode enrichment algorithm, analyze
the six known extra modes on the private fixed 31-spring matrix. The
standalone script `test/mfamg_mode_audit.py` uses only Python 3
standard library; it streams the sparse Matrix Market coordinate file,
so it does not retain the ~3.6M sparse matrix entries in memory.

From `$HOME/lis-v8` in MSYS2 UCRT64 (requires `python3`):

```bash
mkdir -p "$HOME/mfamg_replay/v8_mode_audit"
python3 test/mfamg_mode_audit.py \
  "$HOME/mfamg_replay/challenge_31/A04.mtx" \
  "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
  "$HOME/mfamg_replay/challenge_31/modes04_12plus6negative.mtx" \
  --csv-out "$HOME/mfamg_replay/v8_mode_audit/modes.csv" \
  | tee "$HOME/mfamg_replay/v8_mode_audit/summary.txt"
```

The report checks exact first-12 correspondence (relative tolerance
1e-10), original mode norms and Rayleigh quotients, normalized
full-operator eigenvector residual estimates, independence from the
baseline 12-dimensional span, rank and orthogonality after two-pass
modified Gram-Schmidt, and the six-dimensional Ritz values of the
**symmetric part of the projected operator** on the extra directions.
The projected symmetry diagnostic does not certify that the entire
CCX matrix is symmetric. Rayleigh/short Ritz results are descriptive
and do not prove the extra vectors are actual eigenvectors or that
they can be constructed online; further algorithm design and real
Newton-step validation remain necessary. Only the user-controlled
local CSV and stdout summary are written; no snapshot file is modified.

Public independent CI job `mode-audit-smoke` tests a synthetic
five-row indefinite system, negative projected Ritz values,
non-eigenvector residual, baseline mismatch and unsupported input.
No private CCX matrices are used or uploaded. PR remains Draft.


## OFFLINE projected-Ritz selection experiment (5 vs 6 extra modes)

The first real 31-spring mode audit gave six independent extra
directions with NEGATIVE individual Rayleigh quotients but only
**five negative** and **one positive** Ritz values in the symmetric
part of their six-dimensional projected operator:
`-925.680248743,-463.766476123,-160.55219208,-74.6704595587,-5.84714873376,+101.892963799`.
The extra modes are nearly orthogonal to the original 12, but their
large full-operator normalized eigen-residuals (0.58..0.88)
mean they must NOT be called exact eigenvectors of the full A.

The new `test/mfamg_mode_ritz_selection.py` accepts that SAME
offline 18-mode file and rotates its 6 extra directions into
Ritz vectors of the symmetric projected operator. It writes two
new local Matrix Market arrays without overwriting existing files:

- `modes_negative_ritz.mtx`: original 12 + 5 negative-Ritz directions
- `modes_all_ritz.mtx`: original 12 + all 6 rotated directions
  (same 18-dimensional span as original input)

This is a controlled experiment: the second output distinguishes
any effects caused by rotating/orthogonalizing the basis itself
from effects caused by excluding the positive Ritz direction.
The rotation does not recover modes online, and negative curvature
of the projection does not certify global negative eigenvectors.
Never replace the production near-nullspace basis with these
research-only arrays.

From an up-to-date `$HOME/lis-v8` checkout, using Python 3:

```bash
mkdir -p "$HOME/mfamg_replay/v8_mode_audit/ritz"
python3 test/mfamg_mode_ritz_selection.py \
  "$HOME/mfamg_replay/challenge_31/A04.mtx" \
  "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
  "$HOME/mfamg_replay/challenge_31/modes04_12plus6negative.mtx" \
  --output-dir "$HOME/mfamg_replay/v8_mode_audit/ritz"
```

The two real-snapshot runs described below were completed locally
with the same MFAMG factor=4, k=8, restart=180 and maxiter=3000.
Their interpretation remains diagnostic, not evidence for
replacing offline enrichment with online mode generation.
An independent Linux CI job checks both files on a synthetic
indefinite example and rejects existing-output overwrites.


## Real 31-spring Ritz selection experiment: 17 vs 18 modes (2026-11-08)

A real Windows UCRT64 replay on the private 60,265-DOF, 3,623,639-entry
31-spring A04 operator tested the new Ritz-selected basis and full
Ritz-rotated control against the previously recorded unrotated
18-mode enrichment (same MFAMG factor=4, FGCRO-DR(k=8),
restart=180, maxiter=3000, tolerance=1e-5):

| Mode basis | Status | Iterations | true_relative | Wall time, s |
| --- | ---: | ---: | ---: | ---: |
| Original 12 + 6 offline modes | **0 (converged)** | 2772 | 9.986595263773e-6 | 226.086 |
| Original 12 + 5 negative Ritz directions | **4 (MAXITER)** | 3001 | 4.193222729068e-5 | 235.031 |
| Original 12 + all 6 rotated Ritz directions | **0 (converged)** | 2805 | 9.950907090924e-6 | 227.719 |

For the 17-mode case, recycling was active (builds=1, updates=16,
uses=1, rejects=0), with recycle validation rc=0,
orth_error=9.769962616701e-15 and relation_error=2.061032868769e-9.
Despite valid algebra the solver did **not** reach the tolerance.

For the rotated 18-mode control, recycling was active (builds=1,
updates=15, uses=1, rejects=0), validation rc=0,
orth_error=1.587618925214e-14 and relation_error=7.848487326205e-9.
It converged after only 33 more iterations / 1.633 s more than the
unrotated 18-mode case, supporting robustness to this basis rotation
at this one snapshot. Note that AZ=C error is close to the 1e-8 gate.

**Important inference:** retaining the extra sixth Ritz direction
(the positive-curvature direction in the 6-dimensional symmetric
projected operator) matters for convergence *within 3000 iterations*
in this setup. Selecting only negative projected Ritz directions
fails here, although each of the six original offline vectors had a
negative individual Rayleigh quotient. Neither individual Rayleigh
sign nor projected Ritz sign alone is a sufficient universal rule
for enrichment. This is not proof that the sixth direction is
indispensable at every iteration cap or for every contact operator.

**Research status:** all six source directions are still produced
OFFLINE. The selector only rotates and drops directions; it cannot
generate them from A and the original 12 modes alone. No production
MFAMG, CalculiX, or v7 change is justified. No private data files
or complete logs are uploaded to the repository. PR remains Draft
and blocked pending online enrichment, repeated time measurements,
and a changed-operator Newton sequence experiment.


## Additional historical repeat of Ritz negative-only vs all-18 (private, earlier log)

The user subsequently supplied terminal logs for another fixed 31-spring
A04/b04 replay, using the already present Ritz-selected mode files. The
attempted `mfamg_mode_ritz_selection.py` invocation printed
`RITZ_SELECTION_ERROR=output files exist; choose a new empty output directory`.
This is the intended non-overwrite safeguard: **those files were NOT
regenerated** during the pasted session. Their contents/binary hashes
were not independently re-verified from the terminal transcript.

| Basis | Earlier wall_s | Additional wall_s | Status | Iter | True relative residual |
| --- | ---: | ---: | --- | ---: | ---: |
| 17 modes = baseline + five NEGATIVE Ritz; **omit direction 06** (+101.89) | 235.031 | 217.365 | 4 MAXITER (both) | 3001 | 4.193222729068e-5 |
| 18 modes = all six rotated Ritz | 227.719 | 219.470 | 0 converged (both) | 2805 | 9.950907090924e-6 |

The additional negative-only run had recycle builds=1, updates=16,
uses=1, rejects=0; orth_error=9.769962616701e-15,
relation_error=2.061032868769e-9. The additional 18-mode run had
builds=1, updates=15, uses=1, rejects=0;
orth_error=1.587618925214e-14,
relation_error=7.848487326205e-9. Both passed algebra validation.
The additional run returned process exit=1 for the 17-mode MAXITER
case and exit=0 for the converged 18-mode case, consistent with
solver status.

Across these two reported historical runs, wall time changed by
17.666 s (~7.5%) for the unconverged 17-mode case and 8.249 s
(~3.6%) for the converged 18-mode case, while the recorded iteration
counts and true residuals matched to the shown precision. The
nonconverged run must never be ranked as a faster valid solution.
The measured timing variation motivates matched, alternated-order
benchmarking rather than any performance claim.

**DO NOT CONFLATE** this 17-mode negative-only/omit-06 space with the
separate 17-mode **omit-01** space, which converged at iteration 2833
in a different single replay. The additional historical logs here
do **not** satisfy the new omit01-vs-all18 matched three-pair gate:
different 17-mode basis, no pair-order alternation, missing exact
binary/mode provenance. No production or Newton-sequence evidence.

## Offline six-direction leave-one-out selection: five remaining cases

The six independent extra Ritz directions on the 31-spring snapshot
are ordered by projected eigenvalues -925.68, -463.77, -160.55,
-74.67, -5.85, +101.89. The prior 17-mode negative-only variant
omitted direction 06 and failed to converge at maxiter=3000
(status=4, true_relative=4.193222729068e-5). The full rotated
18-mode control converged after 2805 iterations.

The new offline generator test/mfamg_mode_ritz_leave_one_out.py
creates six local 17-mode files named modes_omit_ritz_01.mtx
through modes_omit_ritz_06.mtx by deleting exactly one column
from the already rotated 18-mode file. The original first 12
columns are unchanged. It checks that case 06 exactly equals the
previously evaluated negative-only 17-mode control.

In MSYS2 UCRT64, after updating the research checkout:

    cd "$HOME/lis-v8"
    git pull --ff-only
    mkdir -p "$HOME/mfamg_replay/v8_mode_audit/loo"
    python3 test/mfamg_mode_ritz_leave_one_out.py \
      --baseline "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
      --rotated "$HOME/mfamg_replay/v8_mode_audit/ritz/modes_all_ritz.mtx" \
      --negative-control "$HOME/mfamg_replay/v8_mode_audit/ritz/modes_negative_ritz.mtx" \
      --output-dir "$HOME/mfamg_replay/v8_mode_audit/loo"

    MFAMG_REPLAY_DATA="$HOME/mfamg_replay" \
    MFAMG_REPLAY_EXE="$HOME/lis-v8/bin6000/mfamgreplay.exe" \
    bash test/mfamg_replay_ritz_leave_one_out.sh

The replay runner tests ONLY the five new cases 01..05, fixing
MFAMG factor=4, FGCRO-DR k=8, restart=180 and maxiter=3000.
It logs every result under the local directory
$HOME/mfamg_replay/v8_mode_audit/loo/logs and refuses overwrite.
The synthetic CI smoke validates basis deletion, equality of
the sixth case to the previous control, unchanged baseline,
mock nonconvergence handling and no overwrite.

This is OFFLINE vector deletion, not online mode generation
from the actual operator. It cannot establish production
readiness, Newton-step behavior or repeat-run timing. Never
upload private matrices or raw contact logs; keep PR Draft.


## Real 31-spring Ritz leave-one-out outcomes (2026-11-08)

The Windows MSYS2 UCRT64 user ran all five previously missing
17-mode offline leave-one-out variants against their private
31-spring fixed A04/b04 snapshot (60,265 DOF, 3,623,639 input
matrix entries). Method FGCRO-DR(k=8), MFAMG factor=4,
restart=180, maxiter=3000, tol=1e-5, one timing per case.

Direction order is ascending by the eigenvalues of the symmetric
part of the six-dimensional projected operator, not certified
eigenvalues of the full nonsymmetric Matrix Market general A04.

| Omitted Ritz direction | Projected value (approx) | Status | Iter | true_relative | wall_s | AZ=C relation_error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 | -925.68 | **0 converged** | 2833 | 9.967925888120e-6 | **207.627** | 6.665478388731e-9 |
| 02 | -463.77 | 4 MAXITER | 3001 | 1.015222798738e-5 | 220.828 | 7.106887924046e-9 |
| 03 | -160.55 | **0 converged** | 2998 | 9.934484468416e-6 | 222.154 | 7.455637638979e-9 |
| 04 | -74.67 | 4 MAXITER | 3001 | 3.313684461675e-5 | 220.258 | 4.779255519919e-9 |
| 05 | -5.85 | 4 MAXITER | 3001 | 2.994747563559e-5 | 224.589 | 5.094805373938e-9 |
| 06 | +101.89 | 4 MAXITER | 3001 | 4.193222729068e-5 | 235.031 | 2.061032868769e-9 |
| Full 18-mode Ritz control | all six | **0 converged** | 2805 | 9.950907090924e-6 | 227.719 | 7.848487326205e-9 |
| Original unrotated 18-mode control | all six | **0 converged** | 2772 | 9.986595263773e-6 | 226.086 | 7.687679465279e-9 |

Every new case 01..05 had recycle builds=1, uses=1,
rejects=0, and REPLAY_RECYCLE_CHECK rc=0, with orthogonality
errors 1.376676550535e-14 to 2.642330798608e-14. Cases with
status=4 must NOT be counted as converged even though
solve_rc=0; in particular case 02 was only 1.522% above
the target residual, and case 03 converged only at iteration
2998, almost at the configured ceiling.

**Result:** 17 total modes *can* converge: omitting the most
negative projected Ritz direction (01) converged at iteration
2833 with 207.627 s, 20.092 s / 8.82% shorter single-run wall
time than the rotated 18-mode control (227.719 s). Omitting
direction 03 also converged but near maxiter. Omitting 02,
04, 05, or the only positive-curvature direction 06 did not
converge by 3000. Thus a blanket "keep only negative Ritz
directions" rule fails, and neither magnitude nor sign alone
predicts a direction's effect on convergence.

**Limitations:** these are different offline-selected subspaces,
so this is not proof of universal necessity of any one direction,
nor a demonstrated statistical speedup or online mode selection
algorithm. All results are single run, on one fixed operator
with offline 18-mode input. Next useful validation is repeated
matched runs of 17-mode omit-01 versus rotated 18-mode control,
preferably with alternating test order to reduce machine-state
bias. Larger-scale Newton-step/changed-A testing and online
mode discovery are still mandatory before integration.
The recycle relation errors remain below 1e-8 but approach
that threshold; monitor for numerical instability.
NO MERGE: PR remains Draft; production LIS and CalculiX untouched.

## v8 bridge toward merge readiness (2026-11-09)

PR #47 remains Draft. The independent standalone FGCRO-DR feature
in anishida/lis#214 and its prerequisite MAXITER leak fix in
anishida/lis#213 were **both merged into upstream `develop` (not `master`)** on
2026-11-09 (merge commits 3d3b31ec and 030aa748). The experimental
v8 branch predates the merged per-solver option/profile logic,
so do NOT cherry-pick overlapping numerical solver diffs.
The solver file now matches accepted upstream `develop` #214 byte-for-byte. The remaining gate is verification of the research-only replay option bridge and the final upstream integration path, not waiting for upstream review.

### Matched replay: independent 17-vs-18 mode benchmark

The single historical 17-mode omit01 run (207.627 s) versus rotated
18-mode control (227.719 s) was exploratory, **not** repeatable evidence.
The public test/mfamg_replay_v8_repeats.py now runs matched AB/BA
pairs using the same environment, changing only the mode file:

    python3 test/mfamg_replay_v8_repeats.py \
      --exe "$HOME/mfamg_replay/build-v8-upstream214-20261109-winfix1/test/.libs/mfamgreplay.exe" \
      --matrix "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      --rhs "$HOME/mfamg_replay/challenge_31/b04.mtx" \
      --modes-17 "$HOME/mfamg_replay/v8_mode_audit/loo/modes_omit_ritz_01.mtx" \
      --modes-18 "$HOME/mfamg_replay/v8_mode_audit/ritz/modes_all_ritz.mtx" \
      --out-dir "$HOME/mfamg_replay/v8_mode_audit/repeats_17v18_NEW" \
      --pairs 3 --require-converged

The order is 17/18, 18/17, 17/18. Outputs stay local,
in a newly created output directory: per-run logs, matched_runs.csv
and assessment.json. The parser rejects malformed/missing results,
nonfinite values and contradictions between exit code/status.
It also verifies the replay header (17 or 18 modes, FGCRO-DR k=8,
restart=180, maxiter=3000), stable matrix dimensions/NNZ across runs,
and **actual recycle activity** (builds, updates and uses >= 1).
Passing an AZ=C check with an empty/unused recycle basis is not enough.
Only a positive wall time, status=0, true residual <= tolerance and
passing recycle algebra count as a converged benchmark observation.

The CSV is written incrementally after each completed solve. If a later
replay crashes, has malformed output or is interrupted, earlier local
rows and all completed per-run logs remain in the dedicated output
directory. Do not resume by overwriting it; create a new directory
for a new comparable measurement series.

In addition to separate 17/18 median wall times, assessment.json
contains the per-pair 17-vs-18 relative wall reductions and their
median; these paired statistics are defined **only when all paired
runs pass**. An INCONCLUSIVE result never counts as a timing
win. The paired estimator reduces some order bias but is not
a confidence interval or evidence of a statistically significant
speedup with three pairs. None of these enhancements establishes
changed-operator/Newton-step behavior.

### Observed private 31-spring paired replay (2026-11-09, UCRT64)

**Measured, not predicted:** on the same 60,265-row A04/b04 contact
snapshot, the verified 17-mode *omit Ritz direction 01* basis and
the rotated all-18 Ritz basis each converged in all three paired trials
(MATCHED_CONVERGED, process exit 0). The three matched pairs used
17/18, 18/17, 17/18 order, freshly compiled v8 HEAD `bd798df`,
fixed FGCRO-DR k=8/restart=180/MAXITER=3000/tol=1e-5, with actual
intra-solve recycling and algebra/residual checks enforced by
`mfamg_replay_v8_repeats.py`.

| Paired trial | 17 omit01 (s) | 18 all Ritz (s) | 17 vs 18 wall reduction |
| --- | ---: | ---: | ---: |
| 1 | 249.347 | 221.322 | -12.663% |
| 2 | 196.469 | 196.182 | -0.146% |
| 3 | 241.138 | 200.488 | -20.276% |

- **All six runs converged**; iterations 2839 (17) versus
  2792 (18), constant across repeats of each basis.
- Median wall time: **241.138 s (17) versus 200.488 s (18)**.
  The difference of medians is **-20.276% wall reduction** for 17
  (i.e. 17 is 20.276% slower relative to the all-18 median).
- The median of the *three within-pair percentage reductions* is
  **-12.663%**. This paired statistic is not the same estimator as
  the difference of marginal medians; do not substitute one for another.
- All three pairs favored 18. The earlier isolated 17-vs-18 result
  (207.627 s vs 227.719 s) **must not be cited as a reliable speedup**.
- This is a small, exploratory six-run sample on **one** private
  contact snapshot, not a proof of universal superiority, statistical
  significance, or Newton-step recycling. The offline 18-mode set is
  the preferred benchmark comparator for this snapshot until further
  case studies.

The benchmark's manifest, six detailed logs, `matched_runs.csv` and
`assessment.json` remain on the local Windows machine; **no private
A/b/mode files or raw contact logs were uploaded to the repository**.
These aggregate measurements are the only published results.

UCRT64 reproducibility caveat: the native Windows `python3.exe`
may interpret Bash `$HOME` differently (`C:/Users/...` rather than
`C:/msys64/home/...`). Pass absolute paths to Python via
`"$(cygpath -m "$PATH_TO_INPUT")"`, including `--exe`, `--matrix`,
`--rhs`, both mode paths, and `--out-dir`. Do not rely on
`pathlib.Path.home()` to find MSYS2 private data. On MSYS2/Windows,
`mkdir -m 700` can also fail at chmod even after creating the
directory; use `mkdir` for the build directory without deleting it,
and keep each benchmark output path unique.

### One-command private 31-spring fresh upstream #214 benchmark

The recommended **next private evidence** after source reconciliation
is a paired fixed-operator comparison of FGMRES vs FGCRO-DR on
exactly A04/b04 and all-18 offline Ritz modes. Do **not** reuse the
old `build-v8-abba-20261109` EXE from the pre-upstream solver;
all earlier 17/18 wall-time numbers were measured with that older core.

`test/mfamg_v8_31_upstream214_one_shot.sh` is a guarded UCRT64
procedure. It never checks out, resets, cleans, merges, or overwrites
the user's existing checkout (including untracked `bin6000/`).
It fetches only the remote-tracking experimental Draft branch,
checks that the solver file's Git blob is exactly the accepted
upstream `develop` #214 blob, archives the public source to a new
private source directory, uses fresh private build/results directories,
compiles its own EXE, and runs 3 alternating FGMRES/FGCRODR pairs
with convergence required. It prints only aggregate results;
matrices, modes and per-run logs stay on Windows.

Run this command **locally** in MSYS2 UCRT64 (no other data needed):

    (
      set -euo pipefail
      cd "$HOME/lis-v8"
      BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
      git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
      git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_31_upstream214_one_shot.sh" | bash
    )

Windows/MSYS2 recovery (2026-11-09): the first attempt stopped
during Git archive extraction because tracked `COPYING` is a symlink
to `LICENSE` and the Windows tar could not create it. It did **not**
compile or run benchmarks. The script now excludes the `COPYING`
archive entry and creates a regular `COPYING` file from the complete
`LICENSE` text; the UCRT64 CI explicitly regression-tests this.
The partially extracted original directory is **not removed, modified
or overwritten**. The retry exclusively uses three **new** paths:
`src-v8-upstream214-20261109-winfix1`,
`build-v8-upstream214-20261109-winfix1` and
`v8_fgmres_vs_fgcrodr_31_upstream214_20261109-winfix1`, all beneath
`$HOME/mfamg_replay`. The script fails closed if any of these exist;
no `rm`, `git clean`, `reset` or existing-output overwrite is used.
To review its source without running, replace `| bash` with
`> /some/new/path` (only if that output path does not already exist).

For result triage, report just `V8_31_SOLVER_BENCHMARK_EXIT`,
`V8_SOLVER_PAIR_ASSESSMENT` and the aggregate `assessment.json`.
Private inputs/logs and the manifest with their paths/hashes remain local.
Both solver modes are in the **same NEW v8 executable**, so this is
the narrow algorithm comparison, not a separate historical-v7
binary comparison nor an actual Newton/contact trajectory.

### Real A04/b04 31-spring FGMRES vs FGCRO-DR upstream #214 result (2026-11-09)

The one-shot UCRT64 `-winfix1` benchmark successfully compiled the
byte-pinned upstream #214 FGCRO-DR core and ran three matched AB/BA
pairs on the same private A04/b04 snapshot, **full 18 Ritz modes**,
restart=180, maxiter=3000, k=8 for FGCRODR, and tol=1e-5.

| Pair | FGMRES iterations / status / wall (s) | FGCRO-DR iterations / status / wall (s) |
| --- | --- | --- |
| 1 | 3001 / 4 MAXITER / 129.952 | 2792 / 0 converged / 233.491 |
| 2 | 3001 / 4 MAXITER / 148.620 | 2792 / 0 converged / 239.993 |
| 3 | 3001 / 4 MAXITER / 148.381 | 2792 / 0 converged / 239.901 |

**FGMRES 0/3 converged; FGCRO-DR 3/3 converged.** Median true
relative residuals: 4.289042843312e-5 (FGMRES) and
9.989979821954e-6 (FGCRO-DR). As designed, the matched runner
returned `INCONCLUSIVE_NONCONVERGENCE` and no numerical speedup
estimate (all speedup fields are null). The shorter FGMRES wall times
are due to the early MAXITER stop, **not faster successful solving**.
The result is valuable convergence evidence on **one fixed-A**
private snapshot only; it establishes neither a speedup nor any
Newton/contact-trajectory or production CalculiX claim.

A controlled next step is **the same solver pair with a common larger
iteration limit**, rather than tuning only the FGMRES budget.
The runner `mfamg_v8_solver_pair_compare.py` now supports
`--maxiter 6000` for both variants and strictly checks that the
executable's REPLAY headers actually report that value. The 6000 cap
is within the existing `mfamgreplay` supported 100..6000 range.
Run the guarded local one-pair probe with the *already compiled*
`-winfix1` EXE (no new full build):

    (
      set -euo pipefail
      cd "$HOME/lis-v8"
      BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
      git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
      git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_31_maxiter6000_probe.sh" | bash
    )

This new probe creates only new
`$HOME/mfamg_replay/v8_pair_6000_python_20261109` and
`$HOME/mfamg_replay/v8_pair_6000_probe_31_20261109` locations;
it refuses to overwrite any existing path and retains the previous
3000-iteration benchmark as-is. It fetches **only the three public
Python parser/driver source files**, not contact inputs; every matrix,
mode, JSON manifest and full log remains local. One pair is just a
convergence feasibility probe, **not a reliable timing speedup claim**.
If both converge at maxiter=6000, conduct a separately named repeated
AB/BA study using the same budget. Do not publish private logs.

### Real A04/b04 follow-up: same-binary common MAXITER=6000 probe (2026-11-10)

The user ran the one-pair private UCRT64 probe with the already built
accepted upstream #214 `-winfix1` EXE, fixed A04/b04, the identical
all-18 Ritz basis, restart 180, tolerance 1e-5 and **common
`LIS_MFAMG_REPLAY_MAXITER=6000` for both solvers**.

| Method | Status | Iterations | True relative residual | Wall time |
| --- | --- | ---: | ---: | ---: |
| FGMRES control | 4 (MAXITER) | 6001 | 3.65628547731e-05 | 256.392 s |
| FGCRO-DR v8, intra k=8 | 0 (converged) | 2792 | 9.989979821954e-06 | 197.002 s |

`MAXITER6000_PROBE_EXIT=1` and
`V8_SOLVER_PAIR_ASSESSMENT=INCONCLUSIVE_NONCONVERGENCE` are **correct**.
The benchmark intentionally reports no speedup because FGMRES is
nonconverged; shorter/longer wall time on a failed solve is NOT a
valid speed comparison. A single pair adds feasibility evidence, not
repeatable timing statistics. The previously recorded 3-pair
maxiter=3000 test also showed FGMRES 0/3 converged versus FGCRO-DR
3/3 converged, with the same 2792 FGCRO-DR iterations. At maxiter
6000, FGMRES residual remains above the specified 1e-5 tolerance.

**Interpretation:** On this one fixed real A04/b04 contact snapshot
and the same 18-mode basis, the accepted upstream FGCRO-DR reaches
the requested accuracy within the iteration budget and FGMRES does
not, including with twice the original MAXITER. This is a
**convergence / iteration-budget advantage**, not a proven wall-time
speedup, broad contact robustness, or CalculiX production result.

**Next private gate:** do not keep raising MAXITER solely to force a
speedup comparison. Prefer a separate 116-spring validation and
authentic consecutive Newton operators (A1/b1 → A2/b2 → A3/b3)
if available; keep exact status/true residual/iter and recycling
algebra checks. The private inputs, manifests, per-run logs and
output folders remain local. PR #47 remains Draft / NO MERGE.

### Separate 116-spring real-contact convergence/regression gate (prepared)

The **116-spring** contact snapshot must not reuse the 31-spring
numerical configuration. The historical 116 control used
`control_116/A02.mtx`, `control_116/b02.mtx`,
`control_116/modes02.mtx` (12 columns), `MFAMG_FINE_AGG_FACTOR=2`,
`restart=180`, `MAXITER=1000`, `tol=1e-5`; the previously reported
FGMRES baseline converged in **592 iterations**. This historical
number is a reference to check, not a new upstream-#214 measurement.

The paired runner now accepts explicit `--fine-agg-factor 2`
(default remains 4 for the separate 31-spring experiments) and
`--maxiter 1000` (default remains 3000). Both are applied to FGMRES
and FGCRO-DR, logged in the local input hash manifest and enforced
through the replay header's iteration-limit check.

`test/mfamg_v8_116_upstream214_one_shot.sh` performs three alternating
FGMRES/FGCRODR pairs using the **same previously built `-winfix1`
upstream-#214 executable**. It checks source Git blob identity, checks
the three local private input files, fetches only three public Python
test scripts into a **new private** scripts folder, and refuses
preexisting output folders. It never modifies the existing checkout,
recompiles the solver, deletes files, or uploads private data.

Run only in MSYS2 UCRT64:

```bash
(
set -euo pipefail
cd "$HOME/lis-v8"
BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_116_upstream214_one_shot.sh" | bash
)
```

The new local-only directories are
`$HOME/mfamg_replay/v8_116_python_upstream214_20261110` and
`$HOME/mfamg_replay/v8_116_fgmres_vs_fgcrodr_upstream214_20261110`.
`V8_116_BENCHMARK_EXIT` and aggregate `assessment.json` can be
shared without private matrices or run logs. A result is classified
`MATCHED_CONVERGED` **only if both methods converge in all pairs**
and FGCRO-DR's reuse and `AZ=C` checks pass. A solver MAXITER,
nonconvergence or invalid recycling yields an inconclusive
assessment with no speedup statistic; it is still useful evidence.
Compare FGMRES iteration counts with the historical 592 baseline,
but do not assume exact parity before the real test.

This gate concerns one fixed 116-spring A/b snapshot, not a
Newton sequence or production CalculiX. Keep the PR Draft.

### Real 116-spring A02/b02 matched result: convergence preserved; extra work loses (2026-11-10)

The user executed `test/mfamg_v8_116_upstream214_one_shot.sh`
with the same already-built accepted upstream develop #214 UCRT64 EXE
on private `control_116/A02.mtx`, `b02.mtx`, and `modes02.mtx`
(12 modes). Identical settings for both algorithms: MFAMG
`factor=2`, `restart=180`, `maxiter=1000`, `tol=1e-5`;
FGCRO-DR intra-solve uses `k=8` and the runner validates its
recycling algebra. All three paired AB/BA controls converged.

| Pair | FGMRES iter/status/wall (s) | FGCRO-DR iter/status/wall (s) |
| --- | --- | --- |
| 1 | 592 / 0 / 15.148 | 455 / 0 / 22.488 |
| 2 | 592 / 0 / 15.203 | 455 / 0 / 22.522 |
| 3 | 592 / 0 / 15.229 | 455 / 0 / 27.261 |

`V8_116_BENCHMARK_EXIT=0` and `MATCHED_CONVERGED`, 6/6 solves.
Median relative residuals were `9.99874391912e-06` (FGMRES) and
`9.980683447558e-06` (FGCRO-DR). The historical **592-iteration**
116-spring FGMRES control was reproduced exactly. The FGCRO-DR solver
uses 137 fewer iterations (23.1%), but wall time is LONGER in every
pair. Median FGMRES wall=15.203 s; FGCRO-DR wall=22.522 s.
The ratio of separate medians gives FGCRO-DR a **48.14% wall-time
overhead** relative to FGMRES; median within-pair overhead=48.46%.
The 27.261 s FGCRO-DR trial increases the time spread; three pairs
are exploratory and do not establish general performance statistics.

**Joint interpretation with 31 springs:** on difficult fixed A04/b04
(18 modes, factor=4), FGCRO-DR converged at 2792 iterations, whereas
FGMRES failed at both common MAXITER 3000 and 6000. That comparison
measures a convergence/iteration-budget advantage only because
FGMRES was not converged. On the easier 116-spring control with
factor=2, both methods converge, and FGMRES is faster. Thus the
accepted intra-solve recycle feature should stay **optional and OFF
by default**; never imply universal speedup, change production
CalculiX or automatically switch all MFAMG solves to FGCRO-DR.

Both data sets are single fixed-operator replays, not a genuine
consecutive Newton/contact sequence. All matrices, modes, private
manifests and per-run logs remain on the local computer. The remaining
real-data gates are authentic A1/b1→A2/b2→A3/b3 and the separately
planned candidate 12/14/18 audit. PR stays Draft / NO MERGE.

### Paired FGMRES control vs FGCRO-DR (same A/b/modes)

To isolate the algorithmic opt-in, the new
`test/mfamg_v8_solver_pair_compare.py` compares `fgmres` against
`fgcrodr` **on the same compiled v8 executable**, using the exact
same A04/b04 and all-18 Ritz modes. No historical-v7 binary is implied
by this default experiment. The runner also supports an explicit
`--exe-v7` for a separately compiled historical v7 FGMRES control,
which it labels separately and checks for the genuinely older log
format. An EXE with v8 `method=` header is rejected when supplied
as the alleged historical v7 baseline.

For manual reruns **only after the one-shot above successfully builds
the new upstream-#214 executable**, use UCRT64 and a unique NEW
output directory; do not use the earlier pre-upstream EXE:

    cd "$HOME/lis-v8"
    ROOT="$HOME/mfamg_replay"
    python3 test/mfamg_v8_solver_pair_compare.py \
      --exe-v8 "$(cygpath -m "$ROOT/build-v8-upstream214-20261109-winfix1/test/.libs/mfamgreplay.exe")" \
      --matrix "$(cygpath -m "$ROOT/challenge_31/A04.mtx")" \
      --rhs "$(cygpath -m "$ROOT/challenge_31/b04.mtx")" \
      --modes "$(cygpath -m "$ROOT/v8_mode_audit/ritz/modes_all_ritz.mtx")" \
      --expected-modes 18 --pairs 3 --require-both-converged \
      --out-dir "$(cygpath -m "$ROOT/v8_pair_fgmres_vs_fgcrodr_NEW")"

Order is FGMRES/FGCRODR, FGCRODR/FGMRES,
FGMRES/FGCRODR. Each run starts a fresh process. The runner cleans
inherited `LIS_*` variables and fixes the shared MFAMG environment;
only the method-specific recycling settings differ. FGCRO-DR uses
k=8, restart=180, maxiter=3000, tol=1e-5, active intra-solve
recycling and checked relation `AZ=C`; FGMRES controls use no FGCRO-DR
recycling settings. Common-mode configuration, actual input and EXE
SHA-256/size, fresh-process records, and per-run logs remain private
under the new output directory.

For an **actual historical v7 vs v8** test, first compile an independent
v7 `mfamgreplay.exe` from the exact original v7 source and the same
toolchain/configuration, then add
`--exe-v7 "$(cygpath -m "$ROOT/path_to_separate_v7/mfamgreplay.exe")"`.
The v7 and v8 binaries must be different files. Differences between
these EXEs may reflect more than the solver algorithm, hence the
same-binary comparison should be examined first.

A positive speedup is reported **only when every paired run of both
methods converged**, has true-relative residual <=1e-5, consistent
method/configuration, and the FGCRO-DR recycling counters/algebra
actually validate. The percentage from separately computed medians
and the median of within-pair percentages are reported separately;
they are different estimators. A MAXITER/control failure yields an
`INCONCLUSIVE_NONCONVERGENCE` assessment with no speedup number.
Unexpected headers, missing replay data, changed inputs, and
preexisting output directories fail closed.

No real FGMRES-vs-FGCRODR comparison has been measured yet.
Passing mock CI asserts behavior of the benchmark runner, not its
physical performance. A fixed A/b replay is not an actual CCX Newton
trajectory or a production integration test. Do not publish any private
matrix, mode or raw contact logs. Windows output directories use
Windows ACLs; do not assume POSIX chmod 0700 there.

### Automatic candidate mode experiment, WITHOUT offline enrichment

test/mfamg_mode_auto_candidates.py is a research-only prototype.
Given the matrix and baseline modes, its default --target small-singular
creates candidate vectors using seeded projected gradient filtering
on the normal operator A^T*A. A second opt-in
--target negative-curvature filter uses the symmetric part
(A + A^T)/2, with projected gradient steps and a mandatory
negative Rayleigh quotient check (v^T A v < 0 after deflation).
This is a curvature diagnostic for a possibly NONSYMMETRIC A, not
an eigenpair certificate, proof of physical negative-energy modes,
or a guarantee of contact convergence. Both targets start from
A and baseline modes ONLY; neither reads the six offline vectors.
The first baseline columns are preserved exactly, and all new
columns are explicitly filtered against the baseline and earlier
accepted vectors. It checks finite operator products, reduction
of the selected objective (A-norm for small-singular; Rayleigh
quotient for negative-curvature), rank, input dimensions and
exclusive output creation.

Start with a **small diagnostic experiment**, not six expensive
full-size modes at once:

    python3 test/mfamg_mode_auto_candidates.py \
      "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
      --output "$HOME/mfamg_replay/v8_mode_audit/auto_candidate_2.mtx" \
      --count 2 --steps 40 --seed 20261109

To compare the two filters on a *private* 31-spring snapshot,
write to a distinct, local output file and audit the generated basis:

    python3 test/mfamg_mode_auto_candidates.py \
      "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
      --target negative-curvature --count 2 --steps 40 --seed 20261109 \
      --output "$HOME/mfamg_replay/v8_mode_audit/auto_negative_2.mtx"

    python3 test/mfamg_mode_audit.py \
      "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
      "$HOME/mfamg_replay/v8_mode_audit/auto_negative_2.mtx" \
      --csv-out "$HOME/mfamg_replay/v8_mode_audit/auto_negative_2_audit.csv"

A filter that cannot verify independent negative curvature fails
rather than emitting fake negative-curvature modes. Do not run a
costly contact replay or change production interpolation until the
audit and a controlled baseline-vs-generated comparison justify it.
The public CI tests two extracted negative directions on a known
symmetric-indefinite diagonal matrix, deterministic bytes and
rejection of an SPD input. These simple tests do not validate CCX.

The script holds sparse coordinates in memory and its pure Python
sparse operator iterations can be expensive for millions of entries.
The resulting directions are **not** certified near-nullspace modes,
not FGCRO-DR harmonic Ritz eigenvectors, and are not integrated into
MFAMG's production interpolation. A lower A-norm alone does not
prove convergence improvement on a nonsymmetric/indefinite contact
operator. Validate rank and run the same fixed replay comparison
before any claim. Synthetic smoke checks mode novelty and output
non-overwrite; it does NOT constitute contact evidence.

### Guarded upstream #214 small-singular 12/14/18 feasibility stage

An optional **one-shot** for the 31-spring private fixed A04/b04
control is now available: `test/mfamg_v8_auto14_upstream214_one_shot.sh`.
It uses the previously compiled `-winfix1` **accepted upstream #214**
EXE, without rebuilding or modifying the user's LIS worktree.

Before doing any expensive Python iterations, it verifies that
`challenge_31/modes04.mtx` has exactly 12 columns and that the
previously verified Ritz `v8_mode_audit/ritz/modes_all_ritz.mtx`
has 18 columns with exactly the **same first 12**, preventing
accidental comparisons against an unrelated baseline. The generator
receives **ONLY A04 and the original 12-mode baseline**:
`--target small-singular --count 2 --steps 60 --seed 20261109`.
It does **not** use the offline six Ritz vectors to construct
the two new candidates. (An earlier negative-curvature real-data
attempt was inconclusive; small-singular is a separate experiment.)

The script puts the generated 14-mode file, archived public Python
scripts and comparator logs in new, separate private paths; it
refuses to overwrite anything. The comparison is one fresh solve
each at 12/14/18 modes, all with FGCRO-DR, k=8, factor=4,
restart=180, maxiter=3000 and independently measured true residual.
The offline 18-mode control is given ONLY to the comparator, not
the generator. `CONTROL_CONVERGED` is not enough to claim automatic
14-mode success: check `metrics.generated14.converged` and
`metrics.generated14.all_converged` before drawing conclusions.

If the generator cannot establish independence and A-norm decrease,
the script exits nonzero and **does not invent modes** or start
the expensive solver comparison. The generator's pure-Python
3.6-million-entry matrix filtering can be resource-intensive.

To run later on the user's own MSYS2 UCRT64 machine:

```bash
(
set -euo pipefail
cd "$HOME/lis-v8"
BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_auto14_upstream214_one_shot.sh" | bash
)
```

Only aggregate `V8_AUTO14_GENERATOR_EXIT`, `V8_AUTO14_COMPARE_EXIT`
and the local `assessment.json` are suitable for later sharing.
All raw A/b/modes and detailed solver logs remain private.
This is **not** a completed benchmark, speedup claim, validation
of 12+2 for contact, or production CalculiX integration.

### Local controlled replay of automatic candidates (12/14/18)

A fixed-operator comparison runner is now available:
`test/mfamg_replay_v8_candidate_compare.py`. It accepts:

- **baseline12** = original 12 modes (unmodified).
- **generated14** = those 12 plus 2 **automatic** candidate directions.
- **offline18** = original 12 plus the 6 offline directions,
  ideally the previously verified full 18-mode rotated Ritz control.

Start with the first two automatically generated directions using
`--target negative-curvature` (if that fails to identify two negative
directions, retain the diagnostic failure; do not silently substitute
the offline vectors). The alternative `--target small-singular` can
be evaluated as a separate experimental series.

All files stay local in MSYS2 UCRT64. **Choose NEW paths for outputs**;
the commands refuse overwriting existing results:

    cd "$HOME/lis-v8"
    # Legacy manual variant; prefer the guarded one-shot above.
    # Do not run git pull/reset/clean on a worktree with local experiments.
    python3 test/mfamg_mode_auto_candidates.py \
      "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
      --target negative-curvature --count 2 --steps 40 --seed 20261109 \
      --output "$HOME/mfamg_replay/v8_mode_audit/auto_negative_2_NEW.mtx"

    python3 test/mfamg_replay_v8_candidate_compare.py \
      --exe "$HOME/lis-v8/bin6000/mfamgreplay.exe" \
      --matrix "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      --rhs "$HOME/mfamg_replay/challenge_31/b04.mtx" \
      --baseline "$HOME/mfamg_replay/challenge_31/modes04.mtx" \
      --generated "$HOME/mfamg_replay/v8_mode_audit/auto_negative_2_NEW.mtx" \
      --offline-control "$HOME/mfamg_replay/v8_mode_audit/ritz/modes_all_ritz.mtx" \
      --out-dir "$HOME/mfamg_replay/v8_mode_audit/candidate_compare_NEW" \
      --repetitions 1 --require-offline-control-converged

The preflight strictly verifies column counts, exact unchanged
baseline columns, and linear independence of new columns. All
variants use the same environment, FGCRO-DR k=8, MFAMG factor=4,
restart=180, maxiter=3000 and true-relative tolerance 1e-5.
It records local SHA-256 fingerprints of **every input and EXE**,
detects file mutation during replay, and journals logs/CSV after each
completed run. Results: `manifest.json`, `candidate_runs.csv`,
`assessment.json` and per-run logs, all private.

The first pass executes 12 → 14 → 18 modes. With
`--repetitions 2`, the second pass reverses to 18 → 14 → 12;
these are full **fresh-process** solves on one fixed A, **not a Newton
sequence**. MAXITER statuses are recorded and not treated as success.
`CONTROL_CONVERGED` means ONLY the offline18 control reached tolerance
in every pass, **not** that automatic candidates worked: inspect
`metrics.generated14.all_converged` and its measured residual.
The tool reports comparative wall medians only if both generated14
and offline18 converge in all repetitions. One or two passes are
exploratory and are not statistically validated speedup evidence.

An independent public synthetic CI smoke requires the correct
12/14/18 rank/baseline invariants, nonconverged baseline and
generated branches, missing replay result handling, a converged
offline control, output non-overwrite and manifest generation.
No private snapshots are read in CI.

The separate 17-vs-18 omit01 runner now **verifies the actual
basis relationship** before execution: 17-mode columns must equal
all-18 rotated Ritz columns with the *first extra* direction
(column 13 overall) removed. Optional `--baseline-modes` verifies
the original 12 columns too. Its `manifest.json` stores the
local input/EXE hashes and catches in-run mutation. This prevents
accidentally treating the old **omit06 negative-only** 17-mode data
as the promising **omit01** 17-mode experiment.

### Cross-operator recycling: matrix identity vs true value-generation reuse

**Important clarification:** two calls to `lis_solve` on the same
`LIS_SOLVER` with two separately allocated matrices A1 and A2 do
**not** establish cross-operator recycling. FGCRO-DR intentionally
**REJECTS** the old `Z,C` vectors when the matrix identity or
structure-generation changes. The preconditioner/recycle can then
rebuild during the A2 solve. A successful A2 solve and valid A2 Z=C
relation do not prove the A1 vectors were transported.

The replay now prints `REPLAY_NEWTON_POLICY` with the **delta**
of lifetime counters for the second solve:
`cross_uses_delta`, `keeps_delta`,
`transports_delta`, `rejects_delta`. When A2 is a separate
matrix object, the CI explicitly requires
`cross_uses_delta=0` and positive `rejects_delta`.

A separate *synthetic diagnostic*, intended only for a fixed-sparsity
diagonal perturbation, accepts
`LIS_MFAMG_REPLAY_NEXT_DIAG_SHIFT=-1`: after solving A1, it
modifies **the existing matrix object A1 in place** by applying
`lis_matrix_shift_diagonal(A1,-1)` and then explicitly marks
`lis_matrix_mark_values_changed(A1)`. A1 is thereby changed to
A1+I without changing its identity or graph. The second solve uses
that same matrix object, with the new RHS, and the separate
`NEXT_A` file is retained as an independent reference for the
second true residual. The CI generator constructs exactly
`NEXT_A=A1+I`. This path MUST demonstrate positive
`cross_uses_delta` and valid A2 Z=C, or the regression fails.

### Exact local graph preflight before running private Newton A1/A2 replay

A **read-only, no-solve** preflight is available in the experimental
`mfamgreplay` binary. It uses the **same exact CSR compatibility
check** as the opt-in value replacement, on matrices *after*
`lis_input_matrix` imports them. Unlike comparisons of the raw
Matrix Market coordinate lines, it checks the precise `ptr[]` and
`index[]` layout that FGCRO-DR will see.

    "$HOME/lis-v8/bin6000/mfamgreplay.exe" --compare-csr \
      "$HOME/mfamg_replay/challenge_31/A04.mtx" \
      "/PRIVATE/REAL/A_next.mtx"

The second path is a **placeholder**; use a locally existing,
authentic next Newton/contact matrix. This mode requires only
two matrix files: no RHS, modes or solver configuration. It does
not solve equations, copy matrix values, modify files or upload
any content. Positive stdout is formatted:

    REPLAY_CSR_COMPARE compatible=1 n=N nnz=NNZ changed_nnz=K

The command returns **0 only when** both imported operators are
owned, unscaled, unsplit, finite real CSR with exactly identical
dimensions, number of entries, row pointers and column indices
(including column ordering), and at least **one changed value**.
Failure produces a diagnostic and a **nonzero** exit code; missing
files, changed sparsity, or an unchanged A2 must not pass. This
is a strict sufficient test, not a claim that the matrices are
mathematically dissimilar when CSR indices differ merely in order.
The actual replay will repeat the same guard before updating values.

Run this preflight *before* choosing
`LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1`; if it fails, do **not**
force a keep/transport by pretending graphs match. Use the
separate-object/rebuild control. Even a successful preflight is
only structural eligibility, not proof of FGCRO-DR cross-step
convergence or an (AZ=C) bound. Those require the subsequent
independent A2 residual and recycle-validation checks.

The Linux and Windows synthetic CI require a positive preflight;
Linux additionally verifies a same-size mismatched column graph,
unchanged A, and missing A2 all cause errors. All private graph
diagnostics and inputs remain **local only**.

### General same-graph CSR A1→A2 value replacement (opt-in)

The experimental replay now also accepts:

    LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1

together with BOTH `LIS_MFAMG_REPLAY_NEXT_A` and
`LIS_MFAMG_REPLAY_NEXT_B`. This is a more general
same-matrix-object test than the diagonal-shift example: it uses the
**entire numerical value array from independently imported A2**,
including changed off-diagonal coefficients, as the second operator
in the *existing* A1 matrix object.

This is deliberately **fail-closed**. Before modifying any value, the
runner requires both matrices to be owned, unsplit, unscaled, assembled
real CSR, the same dimension and exact nnz, with **identical CSR row
pointers AND column indices in the same order**. Every numerical
coefficient must be finite, and at least one coefficient must change.
If ANY part of the graph or order differs, it refuses transfer rather
than silently treating the matrices as interchangeable. No graph
relabeling, row sorting, or nonzero insertion is attempted. After
preflight it copies the numeric values into A1 and explicitly calls
`lis_matrix_mark_values_changed(A1)` so that FGCRO-DR considers its
`Z` for `KEEP / TRANSPORT / REJECT` on the new value generation.

The second input A2 remains an **independent reference** for the
`b2-A2*x2` true residual and the post-solve `AZ=C` validation.
`REPLAY_NEWTON_MATRIX mode=inplace_replace_values changed_nnz=N`
proves that the second call used this checked value-replacement path.
`REPLAY_NEWTON_POLICY cross_uses_delta` proves (or disproves)
whether old recycling was used. Neither success nor the existence
of the same solver object is sufficient alone.

Example on two *actually existing* private consecutive snapshots,
run **locally only** in MSYS2 UCRT64 with the regular v8 solver
environment already documented below (never upload private inputs):

    cd "$HOME/mfamg_replay"
    export LIS_MFAMG_REPLAY_NEXT_A="/PRIVATE/REAL/A_next.mtx"
    export LIS_MFAMG_REPLAY_NEXT_B="/PRIVATE/REAL/b_next.mtx"
    export LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES=1
    unset LIS_MFAMG_REPLAY_NEXT_DIAG_SHIFT
    "$HOME/lis-v8/bin6000/mfamgreplay.exe" \
      challenge_31/A04.mtx challenge_31/b04.mtx \
      challenge_31/modes04_12plus6negative.mtx

The A_next/b_next paths above are **placeholders**, not files created
by the project. For a graph-changing Newton/contact step, this mode
**must fail**; use the ordinary separate-object step to obtain the
rebuild control. Neither path validates multiple Newton iterations,
a nonlinear convergence criterion, or changing contact topology.

New public synthetic regressions verify ordinary separate-object
rejection, diagonal KEEP, and genuine **off-diagonal-nonsymmetric**
value replacement with a positive cross-use counter; they also reject
same-sized but graph-incompatible matrices, unchanged operators, and
conflicting mutation options. Linux builds the actual C code under
AddressSanitizer; Windows UCRT64 also exercises the opt-in replacement.

For any general contact A2, **do not supply a guessed diagonal
shift**. This option only represents A2 accurately when the entire
difference A2-A1 equals the declared diagonal shift. It is not a
general matrix replacement, nor does a two-step synthetic result
prove contact-Newton robustness. The general separate-object path
continues to be a rejection-and-rebuild control.

### A1-to-A2 in-process operator-change smoke (same solver)

The opt-in replay executable now supports a **second operator in the
same process and the same LIS_SOLVER object**, for synthetic and
eventually real consecutive snapshots. This is distinct from the
fresh-process 12/14/18 comparison above. Supply BOTH:

- `LIS_MFAMG_REPLAY_NEXT_A`: second Matrix Market matrix A2;
- `LIS_MFAMG_REPLAY_NEXT_B`: corresponding second RHS b2.

The second solve starts with x2=0. It keeps the original MFAMG mode
vectors and solver object, but replaces the matrix and RHS. Both
operators must have identical row counts/distributions. Unless
explicitly running the documented diagonal-shift diagnostic, A2
is a different matrix object and the recycle policy rejects A1's
previously stored directions. The original
first solve still prints `REPLAY_RESULT`; after the second solve,
the program prints `REPLAY_NEWTON_STEP`, recycling statistics and
AZ=C validation, and `REPLAY_NEWTON_RESULT` with an independently
recomputed true relative residual. Failure of either solve, algebraic
validation, explicit A*x or residual tolerance causes a nonzero exit.
It never uploads matrices and never changes production CalculiX.

A public mandatory UCRT64 regression creates a synthetic A2 by
modifying the diagonal of A1, then runs:

- A1/b1 and A2/b2 using the **same** LIS solver (k=1, intra-solve on),
  but with a **separate A2 object**, verifying rejection of old Z;
- A1/b1 followed by an **in-place A1+I** change with the same
  solver and the same A object, checking positive cross-use counters;
- A2/b2 as an independent fresh-solver reference;
- both explicit true-residual checks and second-space AZ=C validation;
- rejection of a missing half of the A2/b2 argument pair.

The original v8 regression still verifies k=0 / FGMRES parity and
single-operator behavior. This is a **two-linear-system operator-change
smoke**, not an actual sequence of Newton contact iterations, nor
evidence of robust transport across arbitrary operator changes.

When (and only when) two PRIVATE consecutive contact snapshots are
available locally, the user may run the experimental replay as follows
using a newly built `mfamgreplay.exe` from this branch:

    cd "$HOME/mfamg_replay"
    export LIS_MFAMG_EXPERIMENT_FINE_AGG=1
    export LIS_MFAMG_EXPERIMENT_DEEP_MERGE=1
    export LIS_MFAMG_FINE_AGG_FACTOR=4
    export LIS_MFAMG_STAR_SCHUR=0
    export LIS_MFAMG_VCYCLE_SMOOTH_STEPS=3
    export LIS_FGMRES_PROJECTED_RELATIVE=1
    export LIS_MFAMG_REPLAY_METHOD=fgcrodr
    export LIS_MFAMG_REPLAY_RECYCLE_DIM=8
    export LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1
    export LIS_MFAMG_REPLAY_VALIDATE_RECYCLE=1
    export LIS_MFAMG_REPLAY_RESTART=180
    export LIS_MFAMG_REPLAY_MAXITER=3000
    # Substitute actual existing local step-2 paths; DO NOT invent data:
    export LIS_MFAMG_REPLAY_NEXT_A="$HOME/PRIVATE_STEP2_A.mtx"
    export LIS_MFAMG_REPLAY_NEXT_B="$HOME/PRIVATE_STEP2_b.mtx"
    "$HOME/lis-v8/bin6000/mfamgreplay.exe" \
      challenge_31/A04.mtx challenge_31/b04.mtx \
      challenge_31/modes04_12plus6negative.mtx

The command is only an EXAMPLE and requires two actual consecutive
snapshots; do not run it with nonexistent placeholders. A meaningful
real Newton validation also requires numerical comparison with an
independent solution at step 2 and several successive operator
updates, which the current two-step smoke does not cover.

### Local adjacent-snapshot sequence runner (automatic private validation)

The new `test/mfamg_v8_snapshot_sequence.py` turns a **user-held** list of
real consecutive contact/Newton A/b snapshots into a reproducible
**pairwise** validation. It does not need files sent to GitHub or to
the assistant. It uses the same fixed FGCRO-DR/MFAMG environment as the
paired 17-vs-18 runner, sanitizes inherited LIS variables and
fingerprints every input, the modes, the CSV and the executable.

Create a local CSV **outside the repository**. For example, after
substituting actual existing files and their correct Newton order:

```csv
step,matrix,rhs
n01,/PRIVATE/STEP01/A.mtx,/PRIVATE/STEP01/b.mtx
n02,/PRIVATE/STEP02/A.mtx,/PRIVATE/STEP02/b.mtx
n03,/PRIVATE/STEP03/A.mtx,/PRIVATE/STEP03/b.mtx
```

Only the header and the format are real. The A/b locations above
are **placeholders** and must not be used as if they existed.
Relative paths in the CSV resolve relative to the CSV directory.
Each step ID must be unique, and the table must have 2–50 rows.

Run with an actual serial, real-valued, locally built v8 replay binary:

```bash
python3 test/mfamg_v8_snapshot_sequence.py \
  --exe "$HOME/lis-v8/bin6000/mfamgreplay.exe" \
  --modes "/PRIVATE/REAL/12plus6_modes.mtx" \
  --snapshots "/PRIVATE/REAL/newton_steps.csv" \
  --out-dir "/PRIVATE/REAL/v8_sequence_NEW" \
  --fresh-control
```

**Every transition runs three checks** (with `--fresh-control`):
(1) exact read-only `--compare-csr A_i A_{i+1}`,
(2) same-`LIS_SOLVER` A_i→A_{i+1} *within that pair*,
and (3) an independent fresh-solver A_{i+1} control.
The preflight selects `NEXT_REPLACE_VALUES=1` **only** when
CSR graphs/index ordering match, every coefficient is finite and A
has actually changed. If the graph differs (or A is unchanged),
it selects the safe separate-object/rebuild path. Invalid imports,
missing preflight markers or unexpected errors **never silently
become a graph-mismatch/rebuild classification**.

Each pair checks the initial and next true residual (≤1e-5),
FGCRO-DR status, k=8, restart=180, maxiter=3000, the
reported policy deltas, and second-step `AZ=C`/orthogonality
errors (≤1e-8). A fresh solver, if requested, must also converge.
A graph-changing pair MUST NOT report cross-operator reuse.
A same-graph pair may legitimately **REJECT** its recycled space
if the solver considers it unsuitable; convergence and diagnostics,
not mandatory KEEP, are the pass criteria. Results are journaled
after each pair so a failure does not erase previous work.
`--continue-on-failure` allows additional diagnostic pairs
after failure; otherwise the first failed pair ends the test.

All outputs stay local in the **new** private directory:
`manifest.json`, `transitions.csv`, `summary.json`,
plus per-pair CSR/replay/fresh logs. The runner never uploads them,
does not overwrite an existing directory, and detects modified
inputs between paired solves. `PAIRWISE_PASS` means all adjacent
windows met their measured gates, not a speedup or solver production
approval. It contains no invented time medians.

**Important boundary:** every adjacent transition starts a NEW
process, keeping the solver only across A_i→A_{i+1} *inside that
transition*. This is **not** one solver persisted across an entire
nonlinear Newton/contact trajectory. It is a local structural and
numerical screening step before actual in-application validation.

Public CI uses a fully mock, synthetic three-snapshot CSV and covers
same-graph transfer, changed-graph REJECT, fresh controls, MAXITER,
malformed preflight, forbidden input changes, and non-overwriting.

### Guarded LOCAL Newton three-step preflight (never invent chronology)

**A filename sequence is not a Newton trajectory.** No authentic
three-consecutive-step user-provenance CSV has yet been supplied for
this real-data gate. Do not silently treat `A01`, `A02`, `A03` or
matrix file timestamps as proof of consecutive Newton iterations.
The true prerequisite is a **privately curated** `step,matrix,rhs` CSV
whose author independently confirms all rows came from one nonlinear
contact solve in the correct adjacent Newton order. It must have at
least three declared rows; paths may be relative to the CSV. A
physically appropriate baseline mode array must be chosen separately.

The new `test/mfamg_v8_newton_preflight.py` is **strictly read-only**:
it checks the CSV's required columns/unique rows, existence of all
inputs, Matrix Market matrix/RHS/mode dimensions, actual LIS importer
`--compare-csr` outcomes and immutable input/EXE hashes before and
after the CSR checks. It refuses an all-unchanged operator sequence
and mismatching or malformed preflight headers. It prints only
aggregate graph classifications, and explicitly labels CSV order
`CSV_ORDER_DECLARED_NOT_INDEPENDENTLY_PROVEN`.

The UCRT64 wrapper `test/mfamg_v8_newton_triplet_one_shot.sh`
uses the already built upstream #214 `-winfix1` EXE, pinned exact
solver blob and separately archived public Python scripts. It is
**preflight-only by default**; if there is no curated CSV among
`$HOME/mfamg_replay/newton_steps.csv`,
`$HOME/mfamg_replay/newton/newton_steps.csv`, or
`$HOME/mfamg_replay/challenge_31/newton_steps.csv`, it prints
`V8_NEWTON_GATE=NEEDS_EXPLICIT_CURATED_STEPS_CSV` and does nothing.
Nothing is generated from numbered A/b filenames.

To inspect the local gate without performing an expensive solve, run
in MSYS2 UCRT64:

```bash
(
  set -euo pipefail
  cd "$HOME/lis-v8"
  BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
  git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
  git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_newton_triplet_one_shot.sh" | bash
)
```

To actually run the first three-step persistent solver window,
the operator's chronology must FIRST be independently confirmed,
and its exact private CSV/modes explicitly selected. For example,
within that same one-shot invocation, set
`MFAMG_V8_NEWTON_STEPS_CSV` to the correct local ordered CSV,
`MFAMG_V8_NEWTON_MODES` to the matching baseline 12-mode file,
`MFAMG_V8_NEWTON_FINE_AGG_FACTOR=4` (or 2 for a separately verified
116-spring series), `MFAMG_V8_NEWTON_EXECUTE=1` and
`MFAMG_V8_NEWTON_CONSECUTIVE_CONFIRMED=1`. Without the explicit
confirmation flag, even a successful structural preflight does
**not** execute the Newton replay. Both the read-only preflight
and the actual replay fail closed on malformed or changed inputs.

The runner executes **one LIS_SOLVER for three consecutive declared
systems** and optionally a separate fresh third solve. New result
directories are private and never overwritten; the preflight and
triplet use different NEW local Python source directories so the
first inspection cannot cause a future execution overwrite.
Python on Windows receives `cygpath -m` paths; Windows ACLs are
respected instead of assuming POSIX chmod. This is still only an
offline simulation of the declared Newton windows, not in-CalculiX
continuous state, an automatic optimization or production evidence.

CI smoke uses only synthetic mock matrices to verify a legitimate
same-CSR/change-graph triplet, rejected duplicate/missing rows,
all-unchanged or malformed preflight and mutated inputs. Until the
real user's chronological CSV exists, the Newton gate stays open.

### Three genuinely consecutive linear systems in the SAME solver

The opt-in C replay also accepts a third operator and RHS via
`LIS_MFAMG_REPLAY_THIRD_A` and `LIS_MFAMG_REPLAY_THIRD_B` (both
required). It uses **one LIS_SOLVER through all THREE solves**
A1→A2→A3, not separate processes for the second/third steps.
If third A3 has the same exactly checked CSR graph as the matrix
actually used on step 2, opt in to
`LIS_MFAMG_REPLAY_THIRD_REPLACE_VALUES=1` to update that SAME
matrix object in place (including changed off-diagonals) and mark
its numerical value generation. With a different graph, omit this
flag: A3 becomes a separate LIS_MATRIX, allowing FGCRO-DR's
policy to explicitly REJECT stale Z/C. **Never override a negative
CSR graph preflight.**

The program prints `REPLAY_SEQUENCE_MATRIX step=3`,
`REPLAY_SEQUENCE_STEP step=3`,
`REPLAY_SEQUENCE_RECYCLE_CHECK step=3`,
`REPLAY_SEQUENCE_POLICY step=3` and
`REPLAY_SEQUENCE_RESULT step=3`. The third true residual
is computed against an independently imported A3. Any failure
of a preceding solve, the third solver, AZ=C check, finite
residual or tolerance makes the entire invocation fail.
Synthetic Linux ASan and Windows UCRT64 checks exercise
both 3-step KEEP and graph-change REJECT.

A second LOCAL Python runner automates this over a genuine private
Newton snapshot CSV (same 3-column format as above):

```bash
python3 test/mfamg_v8_triplet_windows.py \
  --exe "$HOME/lis-v8/bin6000/mfamgreplay.exe" \
  --modes "/PRIVATE/REAL/12plus6_modes.mtx" \
  --snapshots "/PRIVATE/REAL/newton_steps.csv" \
  --out-dir "/PRIVATE/REAL/v8_triplets_NEW" \
  --fresh-third
```

Each adjacent three-row window: checks CSR for BOTH graph
transitions using the *real LIS importer*, chooses separate-object
REJECT or in-place value replacement separately for steps 2 and 3,
executes **one persistent three-step solver**, verifies each
step's status and independent true residual, checks the step-2
and step-3 recycle relation and policy, and optionally solves
A3 afresh in a separate process for validation. Logs and SHA-256
provenance stay local; failed windows are journaled in
`triplets.csv` and `summary.json`. It refuses existing output
directories, detects changed input files, and normally stops on
first failure. Use `--continue-on-failure` for diagnostics and
`--max-windows N` to limit costly runs (default: all windows).

`PERSISTENT_TRIPLETS_PASS` means only that every measured
THREE-step window met its numerical checks. **Each overlapping
window still starts a NEW process**; it is not an unbroken N-step
trajectory through the actual nonlinear CCX/Newton state, nor
evidence of speedup. In particular, real step ordering and changed
contact topology must be verified before drawing conclusions.

### Synthetic cross-operator KEEP / REJECT / TRANSPORT coverage

The 240-row in-process `mfamg_v8_cross_operator_smoke.sh` replay tests
**KEEP** (identical CSR matrix object, changed values) and **REJECT**
(different objects or graph), including three successive solves,
AZ=C algebra and independent residual checks. Its CI label correctly
remains KEEP/REJECT; this smoke does not itself force TRANSPORT.

The additional standalone `testmfamgv8transport.c` uses the accepted
upstream #214 FGCRO-DR solver options directly. It starts with a
120-row nonsymmetric diagonally dominant CSR system, builds a valid
rank-one recycle basis, then applies a *nonuniform diagonal value*
perturbation in the **same matrix object/CSR graph**, marks the
value generation and solves a new RHS with the same `LIS_SOLVER`.
For a deterministic branch test only, `-fgcrodr_keep 1` and
`-fgcrodr_reject 0` direct a nonparallel but finite recycled
subspace to TRANSPORT instead of KEEP/REJECT. This is intentionally
NOT a production recommendation for the threshold policy.

Linux CI [run 42703842134](https://github.com/3rav/lis/actions/runs/42703842134)
verified `FGCRODR_TRANSPORT_POLICY_OK`: `cross_uses=0->1`,
`transports=0->1`, `true_relative=2.967084384051e-06`,
`orth=2.220446e-16`, `AZ=C relation=2.200247e-16`.
The test ran under AddressSanitizer/LeakSanitizer and also checks
changed generation, nonempty valid recycle state and no policy
reject/keep count changes. A separate Windows UCRT64 invocation
is also part of the full workflow; check its latest completion.

This closes the **synthetic branch-coverage gap** only. It does
not show that an authentic changing-contact Newton operator
would pass the same policy gate or that TRANSPORT is faster.

### Pinned v7 vs v8 default FGMRES numerical parity (public CI)

A new Linux CI job `linux-v7-v8-fgmres-parity` independently
checks out the exact historical v7 base commit
`c4bc64eb55cd3ca0ff0a00db5472aafbafb7a122` and the v8 HEAD,
configures and compiles both from their OWN source trees with matching
flags, and runs the actual compiled offline replay executable on
the **same newly generated 240-row synthetic A/b and two near-nullspace
modes**. Both select FGMRES, no FGCRO-DR recycling. It compares solver
status, return code, iteration count, imported dimensions, fixed options
and explicit true-relative residual with numeric tolerance, but **never
compares wall-clock times**. The test cannot use private contact data
and is intentionally separate from the real-data v8 FGCRODR benchmark.

A PASS would protect a narrow historical-v7 baseline parity contract;
it does not demonstrate any contact speedup, general solver behavior
on nonsymmetric contact graphs, or production CalculiX integration.
The actual CI result must be checked on the latest HEAD before claiming
the gate passed. This job adds no merge or deployment action.

### v8 MAXITER/LSan and zero-recycle regression safeguards

A separate `test/testmfamgv8maxiter.c` target (based on the
standalone regression in upstream anishida/lis#213, but named separately
to avoid duplicating the already-merged upstream test target) forces two solves
in the same LIS solver: first a converged solve that builds a nonempty
recycling basis, then a deliberate ordinary MAXITER exit with that
basis active. It verifies active use, nonzero remaining residual and
recycle algebra. The Linux synthetic C job runs this executable under
**ASan + LeakSanitizer with `detect_leaks=1`**, so a leaked `bcoef`
workspace fails CI. The Windows UCRT64 job also runs the same C
functional regression, without claiming Windows LeakSanitizer support.

A separate zero-recycle sentinel tests
`LIS_MFAMG_REPLAY_METHOD=fgcrodr`,
`LIS_MFAMG_REPLAY_RECYCLE_DIM=0`,
`LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1` on a synthetic indefinite
operator against normal FGMRES/k=0. The intra-solve code now treats
k=0 as a no-op rather than allocating zero bytes for its B
coefficient workspace (which is allocator-dependent). The complete
`REPLAY_RESULT` numerical signature, other than wall time, must
match FGMRES. The accepted upstream solver now supplies the per-solver implementation; a separate offline-replay compatibility bridge translates the old test flags.

These synthetic checks do not establish a contact speedup or substitute
for true consecutive Newton data. PR remains Draft until the separate
real-data and upstream review gates are addressed.

### Upstream FGCRO-DR reconciliation (merged `develop`; PR #47 stays Draft)

Official anishida/lis PRs **#213** (active-recycle MAXITER leak fix) and
**#214** (per-solver opt-in intra-solve harmonic Ritz recycling) were
merged on 2026-11-09 into upstream **`develop`, NOT `master`**.

The historical v7 solver source is byte-identical to the base used by
both upstream PRs (blob `0358908c841f2f242cde8bb3d06bca1c41ec0b85`).
PR #47's solver file was replaced with the exact accepted #214 source
at merge commit `3d3b31ec40812ec6b04a82c9192cf9d69093d474`
(blob `815a73ada1d18d970ec9b6d24306d23ff631fb35`).
**No competing implementation was cherry-picked or combined.**

`test/mfamgreplay.c` is the ONLY compatibility bridge for historical
research environment switches. In FGCRO-DR mode it maps
`LIS_FGCRODR_INTRA_SOLVE_RECYCLE`, `LIS_FGCRODR_INTRA_TRACE` and
`LIS_FGMRES_PROJECTED_RELATIVE` to per-solver options
`-fgcrodr_intra`, `-fgcrodr_intra_trace` and
`-fgcrodr_projected_relative`. The upstream source itself reads
no experimental process-global environment flags.

Before marking the source gate complete, check the entire latest CI
after the replacement and verify that the bridge still exercises
real reuse, AZ=C, true residual, no-leak MAXITER, k=0 and
same-object value-generation/changed-operator behavior.
Then repeat private real benchmarks on the NEW executable, because
old 31-spring times were from the replaced experimental core.
Upstream `develop` integration and production CalculiX deployment
are separate decisions; neither Ready nor merge is automatic.

Green synthetic CI confirms only those tests. It does NOT prove
Newton/contact speedup, TRANSPORT branch coverage or production safety.

### Outstanding release gates

1. Confirm upstream #213/#214 source-blob reconciliation, the replay
   option bridge and latest CI; keep the accepted upstream interface.
2. Repeat private 31-spring and 116-spring benchmarks with true
   residuals, paired timing medians and any MAXITER failures recorded.
3. Compare automatically generated candidates against baseline
   (12 modes) and offline-enriched control (18 modes); reject
   candidates if algebraic validation or convergence fails.
4. Test a real consecutive Newton/contact operator sequence and
   recycling invalidation or transport under changing A, if available.
5. Keep public CI independent of private snapshots, perform final
   integration/review, then consider Ready for review. NO MERGE yet.
