# CalculiX-facing LIS: five-package validation guide

This branch is a combined validation snapshot, not an additional upstream proposal. Its source and regression files are derived from experimental integration head 5b6c0c4a32c04f24ab297cd968eff236b888bd56 plus review follow-up commits; it also preserves the existing three USER regression workflows from fork #169 (upstream #256–259). The only added file relative to the combined five packages is this guide. No workflow definition is changed relative to #169.

## Review and publication order

| Package | Fork PR | Dependency |
|---|---|---|
| Configuration, lifecycle, ABI | #187 | #169 baseline |
| MFAMG interpolation/composite | #188 | #187 |
| Linear/Schur/A-P sessions | #189 | #188 |
| Generalized eigensolver reporting | #190 | #169, upstream #257 |
| Native Krylov 64-bit/budgets | #191 | #169 |

#169 incorporates the four pending upstream PRs; when their actual merges are available, rebase the five packages onto the updated develop, verify their relative diffs, and submit sequentially. #179 and earlier experimental stacks are retained for traceability. The new packages are drafts, not proof of numerical validation.

## Build and regression gates

Use a fresh checkout and build directory for each configuration. Follow existing CI scripts for exact compiler/toolchain requirements. Example serial build from this branch:

```sh
autoreconf -fiv
./configure --enable-test --disable-fortran --disable-f90 --disable-mpi --disable-omp
make -j2
make check
```

Check the generated files against the checkout after regeneration. Repeat using the native project options for OpenMP, MPI (with MPI compiler and launcher), 64-bit LIS_INT, QUAD, complex and OOM testing where supported. Run the existing sanitizer and Archer recipes from ci-extended.yml. Do not treat skipped tests as passed. Compare actual job logs and test exit status, not merely workflow-level checks.

Relevant existing regressions: testsolvercapabilities, testsolveroptionreplace, testabi, testmfamg, testmfamgilutp (including --parallel-ras), testnullspacecandidates, testcondense_lis, testmatrixautosession, testmatrixoperatormpi, testgesolvepair, test_krylov_api_contract. Also run the USER callback/serial-contract/BiCGSTAB(l) tests introduced by the prerequisite PRs. Use make check to select supported configurations and follow the existing workflow for MPI tests; do not invoke serial-only tests as MPI support evidence.

## Application migration

* ABI: opt in to installed config using LIS_USE_INSTALLED_CONFIG before including lis.h. Some configurations still require the same explicit compiler flags as the LIS build (notably _LONG__LONG, _LONG__DOUBLE, _COMPLEX and USE_MPI, plus the OpenMP toolchain flag). Compare lis_abi_get_client_info and lis_abi_get_info through lis_abi_is_layout_compatible with LIS_ABI_LAYOUT_FLAGS. This checks sizes and selected flags, not every offset or calling convention. Rebuild consumers when the solver layout changes; never pass objects from mismatched headers/libraries.
* Profiles: lis_solver_set_option_replace_if_changed preserves reusable state on equivalent textual profiles and handles RAS/Schur fields. Runtime half-storage policy is separate: set it explicitly with lis_solver_set_precon_half.
* AMG: use lis_solver_prepare_amg_coarse_basis(A,solver,n,modes,groups,tol,&diagnostics) or lis_solver_set_amg_coarse_basis. These are geometric interpolation modes, not a declaration of exact or near nullspace. Application supplies correct DOF geometry and contact labels; accepted modes are owned copies. Retire caller-owned preconditioners before configuration changes. Use lis_solver_copy_amg_coarse_basis_vector to obtain a copy of a mode, not a raw owned pointer. Serial default precision is the initial preparation scope.
* A/P: use lis_matrix_auto_session_solve_with_precon_matrix for supported unscaled serial profiles. Publish generation changes for both inputs, keep borrowed inputs alive, and reset the session before freeing them. For the distinct-P MFAMG profiles supply an explicit AMG basis. Use FGMRES for MFAMG_ILUTP. Global ILUT/ILUTP uses one worker; explicit RAS uses supported domain limits. MPI's separate collective native CSR/Jacobi entry point is experimental; it requires congruent world communicators, identical configuration, and does not guarantee recovery from arbitrary internal OOM. No distributed USER or MFAMG support is claimed.
* Schur: create/update/solve/destroy lis_sc_solver_session; use REUSE_ALL only for identical matrix and metadata, NUMERIC for matching pattern, FULL after clearing dimension-bound configuration. Check both returned status and report.solution_committed. The full solution is accepted only after physical residual certification.
* Eigensolver: read LIS_GESOLVE_SESSION_REPORT with lis_gesolve_matrix_session_get_report after an attempt; inspect return_code, stage, solver_called, diagnostic availability and results_committed. A failed preflight can leave earlier results accessible, but those are not results of the failed attempt.
* Native Krylov: use lis_krylov_problem64 and lis_krylov_symmetric_budgeted64 for large local real-double vectors and callback budgets. Zero budgets mean no user limit; existing conservative total cap still applies. Handle BUDGET_EXHAUSTED separately from callback failure. The projected nev/ncv remain bounded int. Inertia counts can use int64_t; they do not prove full spectral extraction by themselves.

## CalculiX acceptance and performance

Use identical CCX 2.24 sources and compiler flags, model, load increments, thread count, physical residual tolerance, and hardware for baseline and candidate. Apply the full integrated adapter patch to pristine sources, or the upgrade patch only to the matching previous adapter. The adapter needs the combined LIS API set.

For each model compare physical residual, displacement/stress or eigenpair result, convergence status, Newton iterations, linear iterations, preconditioner rebuild count, elapsed time and peak resident memory. Test static, unconstrained rigid-body modes, multiple bodies, contact, frequency, buckling, dynamics, condensed DOFs, 64-bit and error paths. Use at least three repeat runs after a warm-up for timing; report median and range, and do not compare runs with different convergence tolerances.

Current evidence: source-tree equality was checked for the five packages. Selected adapter logic was tested locally with UBSan and mocks (12 bodies/72 candidate modes). This is not a full CCX build, native LIS numerical test, MPI run, sanitizer suite or benchmark. The local checkout could not access GitHub DNS. Final CI of all five packages and real CCX benchmarks remain required.

## Review follow-up (Claude and Gemini)

The consolidated packages now include explicit AMG child-basis semantics, partial-allocation cleanup, Schur NUMERIC warm-start retention and exact-nullspace rebind requirements, MPI negative-error reduction and null-communicator validation, conservative eigen iteration availability and additional Krylov64/shift-invert tests. These fixes have been committed but are not yet certified by CI. The earlier baseline residual policy, callback/OOM paths, installed configuration isolation, public report versioning, and CI path/skip audits remain open. No new workflow or PR has been introduced.
