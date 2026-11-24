# Optional CalculiX 2.19 → LIS_MATRIX_USER validation adapter

This is **not a new LIS matrix format** and is not linked into the production
CalculiX binary. The two files `ccx_lis_user_adapter.[ch]` demonstrate the
existing `LIS_MATRIX_USER` split-triangular view backed by native CalculiX
storage. The CI harness checks out the *unmodified* `3rav/CalculiX` at
`84d5dd9ca1bdaf6ad6e92c5c0d50b8f1965975ba` for `CalculiX.h`.

## Storage contract

- `n` equations; `nnz` strictly lower entries
- `icol[j]`: count of strictly lower entries for column `j`
- `irow[k]`: **1-based**, strictly increasing row indices per column
- `ad[n]`: real64 diagonal, **borrowed**, not copied
- `au[0:nnz]`: real64 strict-lower numerical values, **borrowed**
- `au[nzs3:nzs3+nnz]`: optional opposite upper values for asymmetric matrices
- `au_capacity`: total allocated length of `au` in **doubles**; required for
  reliable `nzs3` validation, as it is not present in the legacy solver call
- Origin and sparsity are validated *before* LIS sees them. Positive `nzs3`
  does not imply that the source allocation is large enough
- Index metadata is built once in LIS-owned arrays; no CSR duplicate is
  created for eligible DIRECT preconditioners

## Lifecycle

`lis_initialize()` first, then:

1. `ccx_lis_user_open(..., "-i cg -p jacobi ...", &ctx)`, or GMRES for
   nonsymmetric/contact cases
2. `ccx_lis_user_solve(ctx, b, solution, &relative_residual, &policy, &action)`
3. On Newton value changes, call `ccx_lis_user_update(ctx, ad, au,
   au_capacity, nzs3)` **even if the numerical pointers did not change**.
   This advances LIS's numerical generation and invalidates stale data
4. On contact-like sparsity changes, call `ccx_lis_user_rebind(...)`.
   This creates a replacement topology transactionally and resets old solver
   state **before** destroying the old source. An invalid replacement leaves
   the old state usable
5. `ccx_lis_user_close(ctx)` before `lis_finalize()`

All calls are serial, nonconcurrent, default precision; the borrowed `ad/au`
arrays must remain valid until they are updated or the context is closed.
The interface does **not** take ownership of CCX input or output arrays.
The RHS remains unchanged and the solution is only published on success.
The validation adapter scans diagonal/lower/upper numerical values and rejects
NaN/Inf before binding or updating a view; these scans are O(n+nnz) and can
be made configurable in a separately reviewed production integration.
The adapter also prechecks dimensions against the configured LIS_INT range.

**Important borrowed-data caveat:** rejection is transactional for a *new*
candidate buffer or topology. If CalculiX mutates the currently borrowed
`ad`/`au` arrays *in place* and introduces invalid numbers before calling
`ccx_lis_user_update()`, the adapter cannot restore the previous physical
values without making an additional numerical copy. The caller must restore
those arrays before attempting another solve; a rejected update does not
magically roll back external storage.

## Build / verification

The `.github/workflows/calculix-user-path.yml` GitHub Actions workflow compiles
`testccxuserpath` (low-level storage ABI) and `testccxadapter`
(checked adapter). Both run in int32/int64 Release and ASan+UBSan variants,
with the same LIS and CalculiX integer ABI macros.

The adapter remains deliberately outside CCX's Makefile. The existing
`SPOOLES` selection is unchanged. This test does **not** assemble a physical
`.inp` model, validate contact convergence on real meshes, or guarantee
MPI USER collective error agreement.
