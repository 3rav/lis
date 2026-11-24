# LIS GitHub Actions: focused PR checks

The full historical Actions configuration could expand to dozens of independent
jobs for an unrelated pull request, often recompiling the same LIS sources
several times. Native regression **sources and Automake tests are not removed**.

## Automatic pull-request validation

- `ci.yml`: 8 baseline build/test jobs (tracked Autotools files, Linux serial,
  OpenMP and half, 64-bit `LIS_INT`, MPI, and two MinGW configurations).
- Specialized workflows: trigger only for relevant implementation, test or
  workflow paths, not for every PR.
- Repeated 32/64-bit × Release/ASan/UBSan specialized grids use two
  complementary PR profiles: **int32 Release** and **int64 ASan+UBSan**.
- Initial-residual validation also retains one QUAD profile on PRs.
- CSR-transpose and CA-MFAMG similarly select focused PR configurations.

The exact number of checks depends on changed paths. A small generic solver
change without matching specialized paths executes the **8 baseline** checks,
rather than the former stack of unrelated workflows.

## Broader validation

- `ci-extended.yml`: the original expensive FGCRO-DR leak proofs, OpenMP
  TSan/Archer, memory-allocation failure injections, MPI OOM and MFAMG scalar
  combinations; runs on `push` to `develop`, on PRs touching critical paths, or via `workflow_dispatch`.
- Specialized workflows run full mode/index matrices on matching pushes to
  `develop` and manual dispatch where available. Shared headers, build files
  and relevant infrastructure are included in the path filters.
- GitHub requires a workflow to exist on the default branch for manual
  dispatch. Until these workflows are present on `master`, do not rely on
  manual dispatch; use the automatic PR and integration triggers instead.

## Merge policy

Before merging an algorithm/ABI-wide change, require the relevant targeted
PR checks and extended CI where applicable. Full specialized matrices run on
matching integration pushes; if full pre-merge coverage is required and manual
dispatch is unavailable, arrange it on a validation branch with the maintainer.

Branch protection may explicitly require old job names. Update any required
checks when adopting this split; skipped/removed required jobs can otherwise
leave a pull request blocked.

This reduction affects the **CI orchestration only**, not the LIS numerical
algorithms, native test programs, matrix-free USER implementation or the
CalculiX adapter.
