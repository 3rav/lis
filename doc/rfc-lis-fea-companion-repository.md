# RFC: separate `lis-fea` companion repository for reproducible FEA examples

**Status:** proposal; creating this document does **not** create any repository.

## Motivation

LIS is a reusable numerical library. FEM/CalculiX models, Newton-contact matrices,
RHS vectors, DOF maps, near-nullspace examples, numerical breakdown reproductions,
and benchmark scripts should have their own lifecycle, rather than adding large
test inputs to the LIS source tree or broadening the scope of #234.

Suggested repository: **`anishida/lis-fea`** (name subject to maintainer approval).
This is a **companion project**, not a mandatory LIS dependency.

## Suggested scope

- Small, deterministic synthetic sparse systems (Matrix Market A/b, expected
  residual and failure-mode tests); CI can run these without large downloads.
- Minimal CalculiX example meshes and input files where redistribution is approved.
- Export and map scripts: matrix/RHS, active DOFs, element connectivity, solver
  options and near-nullspace vectors.
- Reference numerical configurations: standalone MFAMG, native ILUTP,
  and MFAMG+ILUTP with global or opt-in RAS correction.
- Reproduction recipes: exact LIS commit/variant, matrix hash, RHS hash,
  FGMRES restart, tolerance, thread count, true residual, setup/solve times,
  max RSS, and expected convergence/breakdown.
- Optional large real datasets only with explicit data-owner authorization.

## Initial layout

```text
lis-fea/
  README.md
  CONTRIBUTING.md
  LICENSE                 # chosen by maintainers
  .gitignore
  docs/
    data-policy.md
    benchmark-format.md
  examples/
    synthetic/
    calculix/
  benchmarks/
    manifests/
    results/
  prototypes/
    c3d8i-condensation/  # separate, model-informed research example
  scripts/
  .github/workflows/
```

## Example numerical regressions

Current **privately held** experiments associated with #234 include 15 frozen
CalculiX contact/Newton linear systems with ~60,265 unknowns each, right-hand
sides, DOF maps and mesh components. The results show both convergent cases
and a repeatable eliminated zero-pivot breakdown in native ILUTP. These are
valuable *potential* regression cases, not automatically redistributable inputs.

Propose initial public tests based on synthetic data and small, newly
constructed minimal reproducers. Include the exact failure expected where
appropriate (a clean `LIS_BREAKDOWN` is a valid regression outcome).

### Data and licensing gates

**Do not upload any private A/b snapshots, DOF maps, mesh geometry, mode
vectors, solutions or original archives under this RFC.** The full source
datasets require explicit owner consent, a license/provenance review and a
check for engineering/confidential details before publication.

The archived snapshots are hundreds of MB compressed, and individual expanded
Matrix Market files can exceed GitHub's normal per-file size limit. Do not put
them in regular Git history. If authorized, store immutable versioned datasets
in Releases or a suitable research data archive, publish SHA-256 checksums,
and keep a small machine-readable manifest in this repo. Git LFS is optional
and subject to quotas and bandwidth costs.

## Relationship to LIS PR #234

- #234 remains a general, experimental **opt-in MFAMG+ILUTP preconditioner**.
- A specialized C3D8I static-condensation adapter, if published later,
  is an optional example, **not** a dependency or universal assumption.
- Large numerical matrices should not be added to #234 or its CI.
- This proposed repository has a distinct maintenance, governance and
  contribution policy.

## Maintainer decision requested

Would you be willing to create `anishida/lis-fea` as a companion repository
(or nominate a preferred name), and accept an initial contribution consisting
only of documentation, synthetic examples, a manifest schema and reproducible
scripts? Real datasets would be a separate, explicit decision.

If preferred, the project can be incubated in a contributor-owned repository
and transferred or linked once its scope and maintenance have been agreed.

**Note:** GitHub pull requests cannot themselves create a second repository
under another owner's namespace. This is a request for design/governance
approval, not a request to merge the examples into LIS.

Related: #234.
