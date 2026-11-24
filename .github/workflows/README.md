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
  combinations; runs on `push` to `develop` or via `workflow_dispatch`.
- Specialized workflow `workflow_dispatch` (and existing push events)
  retains its original full mode/index matrices.

## Merge policy

Before merging an algorithm/ABI-wide change, run the corresponding
specialized workflow at full scope. For sweeping `include/**`, memory
allocator or build-system changes, use extended CI as well.

Branch protection may explicitly require old job names. Update any required
checks when adopting this split; skipped/removed required jobs can otherwise
leave a pull request blocked.

This reduction affects the **CI orchestration only**, not the LIS numerical
algorithms, native test programs, matrix-free USER implementation or the
CalculiX adapter.
