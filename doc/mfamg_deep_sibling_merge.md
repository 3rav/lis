# Opt-in MFAMG deep-level sibling grouping

MFAMG can form a long coarse hierarchy when many short branches attach
to a shared hub. This optional, per-solver grouping step reduces the
number of aggregates that choose the same strong neighboring anchor
before the multivector QR prolongator is constructed.

## Configuration

```text
-i gmres -p mfamg -mfamg_deep_sibling_merge 1
```

The option is `0` (default; unchanged MFAMG hierarchy) or `1`.
Invalid values are rejected. Configuration is stored per solver and
participates in transactional option-profile replacement and change
detection, with no environment-variable global state.

## Algorithm and safeguards

On serial multivector hierarchy levels with `level_id >= 4`, after the
existing connected aggregation:

- A candidate anchor edge must be strong relative to the **diagonals
  of both endpoints**:
  `|a_ij| >= LIS_MFAMG_STRENGTH_THETA * sqrt(|a_ii| |a_jj|)`.
  The existing threshold is 0.05; zero, invalid, and weak edges do not
  qualify.
- An aggregate selects its strongest qualifying anchor. Equal
  magnitudes use the **lowest aggregate ID**, independently of sparse
  row traversal order (for a fixed aggregate numbering).
- Siblings sharing the same anchor can be packed up to the bounded
  target `max(existing_target, min(local_n, 8 * mode_count))`.
  The multiplier is capped before multiplication to avoid integer
  overflow. The minimum depth `4` and multiplier `8` limit the cost
  of enlarged local multivector QR calculations; they are chosen
  heuristics, not universal optimal parameters.
- The normal local multivector QR/rank construction still follows.
  No near-nullspace mode is removed explicitly.

The transformation **does not run on matrices with more than one MPI
rank**. Those matrices use the same distributed aggregation path as
`-mfamg_deep_sibling_merge 0`, rather than attempting unvalidated
cross-rank or rank-local grouping.

A shared anchor does not imply direct connectivity between siblings.

**The anchor (hub) itself remains in its own group**, rather than
joining the siblings that point to it. The selection threshold uses
diagonal-normalized strength, but the winner among qualifying edges
is ranked by raw `|a_ij|`, not a normalized score. Therefore the
anchor choice can change under diagonal scaling. Two aggregates
pointing to each other are handled as separate anchor lists and are
not automatically combined into one sibling group. These are
intentional, documented heuristics, not invariance guarantees.

A shared anchor does not imply direct connectivity between siblings.
Consequently, some grouped aggregates may be disconnected in the
coarse graph. The strong-edge filter limits weak coupling but does
**not** prove an energy-stable interpolation for arbitrary operators.
The option remains disabled by default.

## Numerical regression coverage

`test/testmfamgdeepmerge.c` covers:

1. A 121x121 shifted star Laplacian, verifying the structural grouping
   count (3) and group-size bound.
2. A **1025x1025 SPD star-of-chains** operator (32 arms, 32 unknowns per
   arm, diagonal shift `1e-3 I`) with two smooth near-nullspace modes.
   Preconditioner construction verifies that depth `>=4` is reached,
   grouping actually fires, and the next-level dimension decreases.
3. Full `lis_solve()` with GMRES/MFAMG for option states
   `0 -> 1 -> 0` on **the same LIS_SOLVER**. Each solve checks success,
   finite positive iteration count and independently evaluates the
   physical relative residual `||b-Ax||_2/||b||_2 < 1e-8`.
4. The profile transition back to zero must restore the original
   hierarchy dimensions, including level count; internal per-instance
   counters verify that no merge runs when disabled.
5. An MPI **two-rank numerical solve** on a larger SPD distributed
   line exercises a hierarchy deeper than level 4 and verifies the
   serial-only guard: zero merge attempts, identical hierarchy
   dimensions and physical residuals for `0 -> 1 -> 0`.
   A separate validation-only negative control removes the outer
   `nprocs==1` check: the test must then fail upon reaching the
   helper's independent rank guard. The negative control is never
   included in the submitted production code.

The internal `lis_precon_mfamg_get_deep_merge_info` test accessor
reports attempts, successful events, reduced aggregate count and
deepest affected level. It is **not part of the public LIS API**.

### Example on Linux serial (one deterministic matrix)

| Setting | MFAMG levels | Unknowns at level 5 | Merged aggregates | GMRES iterations | Relative physical residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| Disabled | 13 | 30 | 0 | 37 | 7.63e-10 |
| Enabled | 6 | 6 | 26 (depth 4) | 37 | 8.55e-10 |
| Disabled again | 13 | 30 | 0 | 37 | 7.63e-10 |

This demonstrates a **smaller hierarchy**, not a measured reduction in
Krylov iterations or runtime. The option remains application-dependent;
there is no claim about contact matrices or universal performance.

## Larger single-node CPU-time comparison

A separate benchmark variant uses a **16385-unknown SPD star-of-chains**
matrix (32 arms of length 512), the same two-mode near-nullspace and
GMRES/MFAMG solver, with options `0 -> 1 -> 0`. It reports CPU time
from `clock()` for one explicit preconditioner construction and one
complete `lis_solve()` (which builds and uses the preconditioner).
The benchmark is on a dedicated fork branch, not in the default
`make check` runtime.

| Option | Levels | GMRES iter. | Build CPU (ms) | Solve CPU (ms) | Relative true residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| `0` | 15 | 51 | 3.477 | 55.454 | 9.94e-10 |
| `1` | 7 | 52 | 3.319 | 56.962 | 8.27e-10 |
| `0` restored | 15 | 51 | 3.374 | 56.352 | 9.94e-10 |

The enabled run reduced hierarchy depth but needed **one additional
Krylov iteration**. These single-run CPU times are close enough that
no reliable runtime improvement can be established. The benefit is
therefore **structural** on these examples; the feature is not
presented as a general solver-performance optimization.

A **seven-repetition measurement** on a separate CI runner (same
16385-unknown matrix and serial build) provides a less noisy comparison:

| Option | Build CPU, mean ± sample SD (ms) | Solve CPU, mean ± sample SD (ms) | GMRES iterations |
| --- | ---: | ---: | ---: |
| `0` | 2.486 ± 0.017 | 46.249 ± 0.658 | 51 |
| `1` | 2.422 ± 0.014 | 47.143 ± 0.586 | 52 |
| `0` restored | 2.499 ± 0.021 | 45.842 ± 0.211 | 51 |

The enabled configuration shows **no demonstrated end-to-end speedup**
and one more Krylov iteration, despite the shorter hierarchy. CPU time
depends on the CI runner, so these values are not production benchmarks.

Reproducible seven-run benchmark source and CI log:
https://github.com/3rav/lis/actions/runs/43170700768

Initial single-run comparison:
https://github.com/3rav/lis/actions/runs/43170161116

## Memory and platform validation

The integration test destroys all preconditioners, vectors, solver and
matrix after repeated option transitions. It also requires exact
restoration of the baseline iteration count and physical residual
in the deterministic `0 -> 1 -> 0` sequence. Run it with sanitizers
(ASan/LSan/UBSan) as well as the regular Linux serial/OpenMP/quad,
MPI two-rank and Windows UCRT64 `make check` configurations to check
memory safety on the exercised paths.

This change does not modify CalculiX, the default MFAMG configuration
or the already upstreamed terminal block-Schur solver.

## ABI note

Adding `mfamg_deep_sibling_merge` to the visible `LIS_SOLVER_STRUCT`
changes its binary layout. Consumers relying on the public structure
layout must be rebuilt against matching LIS headers and library;
this addition is **not ABI-compatible with an older binary build**.
It does not introduce a new public function.
