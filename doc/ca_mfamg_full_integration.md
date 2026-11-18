# CA-MFAMG full integration (experimental)

This branch combines opt-in CA-MFAMG aggregation from fork PR #94,
the current upstream native MFAMG+ILUTP residual composition, and an
application-driven static condensation API for arbitrary local internal blocks.

## Generality

- Caller-supplied DOF partition labels and near-nullspace vectors are
  independent of element types, block size, number of bodies and numbering.
- Condensation accepts variable-sized local blocks, retains all physical
  couplings, reconstructs eliminated DOFs and verifies the original residual.
- Explicit FULL / NUMERIC / REUSE_ALL Schur reuse modes. Incompatible or
  failed updates do not replace a previously successful session.
- Native mfamg_ilutp passes CA metadata to its MFAMG child, retaining the
  original native ILUTP correction; no external numerical solver is used.
- Callers explicitly construct the Schur session and pass the reduced matrix
  to lis_solve. -p ca-mfamg alone does NOT condense DOFs.

## Limitations

The condensation bridge supports serial real-double unsplit CSR. Distributed
MPI or complex variants return UNSUPPORTED. The CA hierarchy and composite
retain their existing MFAMG/ILUTP constraints, including the applicable
OpenMP policy. No C3D8I, nine-DOF, two-body or individual snapshot tuning.

## Testing

Local C core, independent randomized differential (500 cases), allocation
fault injection and mock adapter tests passed. The combined implementation
passed GitHub Actions on fork commit a056ad8c759265a3e51634d3be2c562fe1cdd3a0:
5/5 native integration jobs (serial, OpenMP, 64-bit LIS_INT, ASan/UBSan,
MPI), 20/20 standard Linux/MSYS2 jobs, and 3/3 native ILUTP regression jobs.
The serial native make check reported 119 PASS, 0 FAIL and 0 SKIP.
GitHub Actions must be run again whenever the reviewed head changes.

**Important MPI scope:** the two full condensation profile scripts explicitly
SKIP in MPI builds because the Schur adapter supports serial real-double
unsplit CSR only. Other MFAMG regressions (including testmfamgfineagg)
run under MPI. Passing the MPI job does NOT establish distributed Schur
condensation or CA-MFAMG+ILUTP composition.

## ABI and API compatibility

The opt-in implementation appends fields to the public LIS_SOLVER_STRUCT,
changing its binary layout. Users of the new release must rebuild their
applications/plugins against the matching LIS headers and binary. This is
source-level backward-compatible with supported options but not a promise
of ABI compatibility with already compiled clients. The maintainer should
assess the library version/SONAME policy before release.
The condensation bridge exposes public headers lis_condense_core.h and
lis_condense_lis.h; its lis_sc_lis_labels API requires explicit LIS_INT
element counts for source and destination buffers.

## Build integration

Configure/build and regression checks pass after autoreconf -fiv.
Generated Makefile.in files for new condensation sources and tests must be
regenerated and reviewed before merging if downstream uses the checked-in
Autotools build files directly.

Private CalculiX matrix snapshots are excluded from this PR.
