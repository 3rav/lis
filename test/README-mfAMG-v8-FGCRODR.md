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

PR #47 remains Draft. A clean standalone FGCRO-DR feature is under
upstream review in anishida/lis#214 and its independent MAXITER fix
in anishida/lis#213. The experimental v8 branch predates the new
per-solver option/profile logic, so do NOT simply merge duplicate
solver implementation diffs. Reconcile once the upstream code is ready.

### Matched replay: independent 17-vs-18 mode benchmark

The single historical 17-mode omit01 run (207.627 s) versus rotated
18-mode control (227.719 s) was exploratory, **not** repeatable evidence.
The public test/mfamg_replay_v8_repeats.py now runs matched AB/BA
pairs using the same environment, changing only the mode file:

    python3 test/mfamg_replay_v8_repeats.py \
      --exe "$HOME/lis-v8/bin6000/mfamgreplay.exe" \
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

### Outstanding release gates

1. Reconcile the solver implementation with upstream #213/#214;
   retain default behavior and all new regression safeguards.
2. Repeat private 31-spring and 116-spring benchmarks with true
   residuals, paired timing medians and any MAXITER failures recorded.
3. Compare automatically generated candidates against baseline
   (12 modes) and offline-enriched control (18 modes); reject
   candidates if algebraic validation or convergence fails.
4. Test a real consecutive Newton/contact operator sequence and
   recycling invalidation or transport under changing A, if available.
5. Keep public CI independent of private snapshots, perform final
   integration/review, then consider Ready for review. NO MERGE yet.
