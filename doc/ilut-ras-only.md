# Experimental opt-in RAS ILUT (without Schur)

## Scope

This narrow review branch isolates **Restricted Additive Schwarz (RAS)** from the earlier RAS + Schur proof of concept. It intentionally **does not contain** the approximate dense or sparse interface Schur operator, diagonal interior elimination, promotion of weak-diagonal KKT rows, or any claimed two-level correction.

Options (serial, double-precision CSR only; MPI and adjoint `psolveh` are **not supported**):

| Option | Default | Meaning |
|---|---|---|
| `-ilut_ras_overlap 0..2` | -1 (off) | Enable RAS with 0–2 edge-overlap layers |
| `-ilut_ras_domains 2..64` | 2 | Domain count |
| `-ilut_ras_parallel 0\|1` | 0 | Enable OpenMP outer-domain *apply*, **not setup** |
| `-ilut_ras_rcm 0\|1` | 0 | Partition domain owners by global RCM-ranked positions |
| `-ilut_rcm 0\|1` | 0 | Independently RCM-order each local owner+halo ILUT matrix |
| `-ilut_row_scale 0\|1` | 0 | Apply local ILUT row equilibration |
| `-ilut_drop`, `-ilut_rate` | ILUT defaults | Control local factor dropping and fill |
| `lis_solver_set_precon_half()` | off | Opt into local packed FP16 factors |

All choices are opt-in; default global ILUT remains unchanged. Owner-core writes are nonoverlapping. Each local factor is built and applied using a **single internal thread**, while the outer RAS apply can distribute distinct subdomains across OpenMP threads. Setup remains serial because LIS solver/preconditioner creation mutates shared library state. Applying the *same* preconditioner concurrently from independent callers is unsupported.

### Limitations and safety

- Tested only on real-valued, serial, double CSR matrices. RAS is not an MPI domain decomposition. No adjoint `psolveh` is supplied; use a solver compatible with one-sided nonsymmetric preconditioning.
- RCM is an ordering heuristic, not a provably optimal geometric separator, and its two option flags mean different things.
- Factor-entry counts are **not peak RSS**. FGMRES true residuals and full setup+solve wall times must both be reported.
- No automatic diagonal shift or condition-number estimator has been validated. The separate diagnostic helper measures local inverse amplification on an **actual supplied RHS**; a single probe is not a proof of general factor stability.
- If local ILUT produces nonfinite output, setup/apply fails explicitly. A very large but finite response may still be numerically inappropriate; do not assume the defaults are safe for all contact/KKT cases.
- The Schur experiment remains in [fork Draft PR #66](https://github.com/3rav/lis/pull/66), not in this RAS-only proposal.

## Reproducible correctness tests

```sh
./test/testhalfilut --ras-3d-smoke
./test/testhalfilut --ras-rcm-smoke
./test/testhalfilut --ras-half-smoke
./test/testhalfilut --ras-row-scale-smoke
./test/testhalfilut --ras-local-rcm-smoke
```

The synthetic 512-row scalar 3D stencil checks finite output and serial/OpenMP consistency with reordered and row-scaled local factors; it is not real contact FEA. The routine `testhalfilut` additionally checks serial vs OpenMP local RAS apply at overlaps 0, 1 and 2 and an algebraic stabilized contact-KKT-like regression.

## Real matrix benchmark: negative result

Public `Boeing/bcsstk39` structural stiffness, 46,772 DOF, 2,089,294 imported entries, measured on GitHub Actions. It is **not** an indefinite Lagrange-multiplier contact-KKT snapshot.

An actual one-thread global ILUT baseline (rate=16, drop=1e-4) was compared to an eight-thread, eight-domain RAS (rate=24, drop=5e-5, overlap=1, local RCM on, row scaling on). Factor-entry storage was calibrated to **0.334666% relative difference**:

| Metric | Serial global ILUT | Parallel RAS |
|---|---:|---:|
| Stored-factor entry proxy | 18,980,706 | 18,917,184 |
| Iterations | **84** | 501 (500-iteration cap) |
| Original-system true relative residual | **9.78215856e-9** | 1.95363117e-8 |
| Setup / solve wall time | 9.428 / 3.421 s | 9.575 / 14.468 s |
| Total wall time | **12.848 s** | 24.043 s |
| Converged at 1e-8 | **yes** | **no** |

Source logs and archived artifact: https://github.com/3rav/lis/actions/runs/43403263560

**STOP for a claim of performance superiority** at this matched-fill setting. A green benchmark workflow means that the experiments completed and logs were preserved, **not** that both variants converged. Neither a setup-only speedup nor a smaller uncalibrated factor-entry count justifies calling this implementation a proven improvement. Further RAS-only review should focus on correctness and narrowly documented opt-in scope, not on asserting numerical benefit.

## Review criteria

- Full 12-matrix CI and ThreadSanitizer+Archer must pass on this **actual RAS-only branch**.
- The removed Schur option must not be registered or linked in the RAS-only code or tests.
- Preserve default ILUT behavior, OOM cleanup, local FP16 and row-scaled RCM test coverage.
- Keep this PR **Draft** until the branch CI and numerical-safety review are completed; even then, "Ready for review" is an invitation to review, not a guarantee of merge or improved FEA performance.
