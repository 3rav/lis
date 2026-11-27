# Reconstructed historical RCI foundation (qualification only)

Based on the supplied historical transcript (v0.42, v0.46 and v0.51).
Not recovered original source; no claim of its ABI or numerical equivalence.

- Separate size/version-checked opaque handle request with A/B, PRECON,
  shifted solve, zero and complex AXPBY; generation and token guards.
- Block dispatch, with deterministic per-column fallback and error handling.
- Non-MPI contour projection, with complex shifts and weights and caller-owned
  vector storage, using the same machine for standard and generalized cases.
  The caller implements the backend, including shifted solves. No result is
  published on error; scratch can be modified and is not rolled back.
- Synchronous nested real-double solver adapter uses existing four core methods;
  it does not implement complex shifted solves or a native LIS preconditioner.
- Tests: independent diagonal complex reference, block/column dispatch, reuse,
  error/generation/cancellation injection at each phase, nested four methods,
  version, size and aliasing checks.

No public headers installed, no configure/Autotools changes. The source change
adds an internal qualification dimension accessor, to avoid reading beyond an
incompatible child context. Keep this branch separate from the upstream PR until
the new API design and qualification have been reviewed.

Not yet FEAST outer iteration, shared classic solver kernels, GCRO-DR/FGCRO-DR,
stackless Krylov, split-complex bridge or extended precision LIS builds. Complex
coefficients/opaque vectors are not equivalent to full complex LIS support.
