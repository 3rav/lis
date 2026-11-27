# LIS 2.1.25 — BUCKLE fixes (stage 10)

## Scope and provenance

Patches were developed against the **uploaded LIS 2.1.25 source archive**, with a serial real-double OpenMP build. They are not a verified build of GitHub develop commit `91668f3197b58d9c964c1b4d484284d7c365dc8c`. Source inspection of that commit confirmed the relevant NEV condition and the missing typed linear-options propagation in transformed block GLI; rebase and compile/test against the actual develop tree before opening pull requests. No remote repository was modified.

## Apply order

1. `01_block_nev_status.patch`: transformed generalized Block Lanczos only. Keep residual certification for every computed vector, but derive the convergence status from the leading explicitly requested pairs (the transformed solver's decreasing-|nu| order). Report MAXITER when rank deflation returns fewer pairs than requested. Preserve legacy behavior for absent/zero NEV; do not change Ritz-value-only mode or the numerical algorithm. No ABI changes.
2. `02_block_gli_linear_options.patch`: apply `lis_esolver_set_linear_options()` configuration before the global command-line options in transformed Block GLI, consistent with scalar GLI. No ABI changes.
3. `03_gli_selective_shift_csr.patch`: depends on patch 02. Materialize the shifted operator into owned CSR for serial MFAMG, MFAMG_ILUTP and ILUT with pivoting enabled; retain the operator path for other preconditioners and existing distributed path. Do the same when recovering from a shifted-solve failure. Retain the owned matrix until the preconditioner/solver is released. The block path now calls the same preconditioner rebuild helper as the scalar path. No ABI changes.

From a clean 2.1.25 source tree:

```sh
patch -p1 < /path/to/01_block_nev_status.patch
patch -p1 < /path/to/02_block_gli_linear_options.patch
patch -p1 < /path/to/03_gli_selective_shift_csr.patch
```

## Regression evidence

`verified/summary.csv`: 29/29 expected checks passed (10 NEV before/after probes, 12 preconditioner/representation probes, 7 existing library regression executables with all three changes). Logs are in `verified/`. The test fixture uses diagonal synthetic pencils of sizes 80 and 12 with known eigenvalues. The NEV fixture also tests repeated eigenvalues and rank deflation returning fewer pairs than requested. It asserts the solver's reported status (not only the API return), the reported count for the rank-deficient case, and the leading eigenvalue/residual values. The CSR probe uses GNU linker `--wrap=lis_precon_create` ONLY to inject a synthetic near-nullspace vector and observe the matrix representation. The tests do not demonstrate a production near-nullspace propagation API, nor do they establish complete iterative convergence on the large CCX model. No MUMPS is used in these probes.

To reproduce on an unmodified, configured and compiled LIS 2.1.25 source tree (GCC, real double, OpenMP, non-MPI, static `liblis.a`):

```sh
python3 run_regressions.py /path/to/pristine/lis-2.1.25
# or specify the pristine static archive separately:
python3 run_regressions.py /path/to/pristine/lis-2.1.25 --lib /path/to/liblis.a
```

The script compiles alternate eigensolver objects and links them ahead of the static library; it does not edit the source tree. It sets `OMP_NUM_THREADS=1` and `LIS_MFAMG_SPARSE_TERMINAL=1`. A failed assertion terminates the run. The native MFAMG_ILUTP test uses FGMRES; BiCGSTAB with this near-nullspace configuration is rejected by the existing capability contract, so do not remove that validation solely to enable the benchmark.

## Open items before submission

* Build and run the full regression suite and MPI/other precision builds on actual develop; verify selective CSR recovery and failure cleanup under memory-allocation failures and row-backed user matrices. Patch 03 is experimental pending these checks.
* Add a supported ownership-aware eigensolver API for passing application AMG coarse basis / near-nullspace to the inner solver (also in nested refinement and retries). Develop has separate `lis_solver_set_amg_coarse_basis()` semantics; do not conflate this basis with exact nullspace deflation. The stage-10 adapter is test-only, not a proposed public API.
* Full-size iterative CCX convergence remains unresolved. A small regression with successful MFAMG does not prove performance or convergence at 152,592 DOF.
* For a public PR, use only these synthetic fixtures and source changes; do not attach private CCX matrices or private benchmark input data.
