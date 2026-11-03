# Native Krylov: proposed composable API and integration review

This is an **experimental integration branch**, not an upstream pull request.
It reconciles the prototypes in upstream PRs #174 (adaptive subspace),
#175 (reusable allocation), #176 (thresholded Jacobi), and #177
(multi-direction restart seeding).

## ABI first

The original `lis_krylov_problem`, `lis_krylov_options`, and
`lis_krylov_stats` definitions are unchanged. In particular, this branch
**does not append fields** to the two public option/statistics structures.
Previously compiled callers of `lis_krylov_symmetric()` continue to use
the same struct sizes and call signature.

New features use two **new**, explicitly versioned structs and a single
composable entry point. Version 1 requires `struct_size == sizeof(v1)`
and `version == LIS_KRYLOV_EXT_VERSION`; unrecognized layouts are rejected.
Future versions should define and support compatible layouts explicitly,
without blindly reading fields past the supplied size.

## Example: adaptive sizing + multiple restart seeds + workspace

~~~c
lis_krylov_ext_options ext;
lis_krylov_ext_stats extra;
lis_krylov_workspace *workspace = NULL;
lis_krylov_ext_options_default(&ext);
lis_krylov_ext_stats_default(&extra);
ext.adaptive_ncv = 1;
ext.ncv_max = 24;
ext.restart_directions = 3;

if (lis_krylov_workspace_create(&problem, &options, &ext,
                                &workspace) == LIS_KRYLOV_SUCCESS) {
    /* Prepare the appropriate callbacks and inverse for each sigma. */
    status = lis_krylov_symmetric_ex(
        &problem, &options, &ext, workspace,
        eigenvalues, eigenvectors, residuals, &stats, &extra);
    lis_krylov_workspace_destroy(workspace);
}
~~~

For a one-off solve, pass `NULL` as `workspace` to the extended
entry point. For the unchanged legacy operation, use
`lis_krylov_symmetric()`. All extended controls are opt-in.
A default-initialized extended configuration sets
`restart_directions=1` and `adaptive_ncv=0`.

### Workspace and active capacity

`options.ncv` denotes the **initial active subspace** (or the legacy
default when zero). With adaptation enabled, the workspace reserves
`ncv_max` (or the documented automatic maximum) **at creation**,
and projected-matrix stride uses that reserved capacity. The active
dimension can grow without allocating or invalidating its vectors.
Both the initialized `options` and the extended configuration must
describe exactly the reserved problem shape on each reused call.
Changes to dimension, transform, `nev`, initial `ncv`, hard-locking
setting or reserved capacity are rejected. Matvec callbacks, context,
spectral shift and non-shape numerical controls can be refreshed.

The workspace retains **storage, not eigensolutions or factorization**.
Each call resets the convergence state, locked-vector count and
extended counters. The caller must prepare the shifted inverse for each
new value of `problem.sigma`; the callback interface does not offer
a safe way to infer whether a cached factorization corresponds to the
current shift. Workspaces must not be used concurrently.

### Restart seeding and diagnostics

`restart_directions=1` follows the original restart seed sequence.
Larger values attempt additional deterministic, metric-orthogonal
vectors on restarts. These are **scalar** callback operations;
this is not block Lanczos or a multi-right-hand-side solver.
Additional directions can be rejected for near-dependence or low
metric rank. Diagnostics report `extra_seed_attempts` and
`extra_seed_accepted` separately, so a requested count is not
silently mistaken for an achieved count.

`lis_krylov_ext_stats` also records initial/final `ncv`,
growth events, accepted basis vectors, and projected eigendecompositions.
The original `lis_krylov_stats` remains ABI-identical.

### Dense Rayleigh–Ritz

Initial Jacobi sweeps skip pivots below 5% of the largest
off-diagonal entry **at the start of each sweep**. This is a heuristic;
later sweeps retain the original machine-epsilon cutoff and stopping
criterion. The changed rotation order can influence roundoff and
iteration counts. It is orthogonal to workspace lifetime and
restart policy.

## Numerical and performance acceptance criteria

The integration test `test_krylov_unified` checks:

- Independent physical residuals, metric normalization and
  B-orthogonality for repeated and clustered finite eigenvalues
- Equality of the original one-direction path and the extended
  one-direction path, including iteration and callback counts
- Adaptive subspace and multi-direction restart **in the same solve**
  using a reusable workspace
- Repeated solves at distinct shifts, capacity mismatch rejection,
  invalid version rejection and callback error propagation
- Extra seed attempts versus actual accepted directions

`test_krylov_workspace` additionally compares reusable versus
independent one-off solves, and `test_krylov_projected` tests
clustered dense projections. `test_krylov_adaptive_ex` explicitly
forces **actual** subspace growth from `ncv=6` to `ncv=24`
using a reusable workspace; its combined 3-direction restart run
converges all four requested eigenpairs with four growth events.
Existing native and LIS adapter regressions must continue to pass.

On this particular n=64 diagonal test, adaptive scalar restart
used 164 operator applications and adaptive 3-direction restart
used 202. This is **not** evidence that extra restart directions
accelerate the problem; it only validates their compatibility and
numerical correctness. Additional independent failure-to-success
examples would be necessary to justify exposing this control upstream.

Post-integration test run (GCC, Clang, OpenMP, ASan/UBSan and full
Autotools): https://github.com/3rav/lis/actions/runs/41689162906

Performance claims require equal-memory fixed `ncv=ncv_max`
baselines, separate end-to-end and Jacobi-kernel times, solver work
counts, and independent repeated runs on multiple matrix families.
The previous individual PR benchmarks are **not** proof of broad
speedups or of usefulness for every combination.

## Preliminary integrated measurements (single runner)

A standalone, deterministic synthetic chain benchmark compares fixed
initial/mid/maximum `ncv` with adaptive sizing, using the same
reserved maximum capacity for the fixed-maximum and adaptive cases.
The current prototype produced:

| Case | Configuration | Operator calls | Restarts | Process CPU |
|---|---|---:|---:|---:|
| n=80, nev=3 | Fixed ncv=24 | 20 | 0 | 0.000673 s |
| n=80, nev=3 | Adaptive 6→9, cap=24 | 26 | 4 | 0.000382 s |
| n=120, nev=4 | Fixed ncv=24 | 27 | 1 | 0.001510 s |
| n=120, nev=4 | Adaptive 8→12, cap=24 | 33 | 5 | 0.000860 s |

Adaptive sizing did **not** reduce operator calls compared with fixed
`ncv=ncv_max`. These CPU measurements are sub-millisecond to millisecond
and have not been replicated on production-size problems; they do
**not** establish a general improvement. The observed direction
shows only that projected-subspace work may be relevant even when
operator calls increase.

An alternating-order, five-trial comparison of 250 independent
shifted solves with an exact diagonal callback measured:

| Configuration | Median CPU time / solve | Min–max total CPU per trial |
|---|---:|---:|
| One-shot allocation | 572.796 µs | 0.143086–0.144698 s |
| Reused workspace | 571.924 µs | 0.142812–0.144153 s |

**No meaningful benefit from workspace reuse was demonstrated in this
case**: the difference is approximately 0.15%, and observed ranges
overlap. It would be premature to add the workspace API to upstream
based on this performance evidence alone.

Integrated benchmark and regression run:
https://github.com/3rav/lis/actions/runs/41688434138

## Integration and upstream disposition

- This branch intentionally composes all four features in one
  experimental implementation; it is **not** evidence that all four
  original PRs should be merged as they stand.
- Original PRs have pairwise Git merge conflicts (all six pairs).
  These include algorithm-code conflicts among #174, #175 and #177.
- Prioritize review of the internal #176 optimization, then review
  #175 with workload data; #174's original ABI change must be
  redesigned; #177 needs evidence of a convergence benefit.
- No new dependencies, application-specific code or permanent
  GitHub Actions workflows are intended.
