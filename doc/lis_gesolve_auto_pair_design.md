# LIS-C: safe AUTO matrix-pair session for generalized eigensolvers

**Status:** design only. No new numerical API is implemented by this proposal.
**Baseline:** `anishida/lis` `develop`, commit `b7c91ad66f37944c6adc27094f5f44879a105ae0` (2026-10-31).

## Motivation

An embedding application such as CalculiX supplies stiffness and mass matrices
as a pair. Existing LIS cached CSR materialization operates one source at a
time; native shared split-triangular topology applies to source views, not to
two independently materialized CSR outputs.

A generic helper must respect eigensolver state lifetimes. `lis_gesolve_once()`
stores A/B on `LIS_ESOLVER`, and an eigensolver retains work vectors, eigenpairs,
results and possibly nested linear-solver state between calls. A matrix pair
cache alone does not prove all those references are dead before CSR refresh.

## Proposed API (names and signature to be reviewed)

```c
typedef struct LIS_GESOLVE_MATRIX_SESSION_STRUCT *LIS_GESOLVE_MATRIX_SESSION;

LIS_INT lis_gesolve_matrix_session_create(
    LIS_GESOLVE_MATRIX_SESSION *session);
LIS_INT lis_gesolve_matrix_session_get_esolver(
    LIS_GESOLVE_MATRIX_SESSION session,
    LIS_ESOLVER *borrowed_esolver);
LIS_INT lis_gesolve_matrix_session_solve(
    LIS_GESOLVE_MATRIX_SESSION session,
    LIS_MATRIX A, LIS_MATRIX B, LIS_VECTOR x,
    LIS_SCALAR *eigenvalue, LIS_INT *policy_A, LIS_INT *policy_B);
LIS_INT lis_gesolve_matrix_session_reset(
    LIS_GESOLVE_MATRIX_SESSION session);
void lis_gesolve_matrix_session_destroy(
    LIS_GESOLVE_MATRIX_SESSION session);
```

The session owns the eigensolver and both caches, but borrows A/B and the
application's input/output vector. Expose no borrowed cache-owned CSR to callers
except under an explicit invalidation contract.

## Preconditions to implementation

1. **Map retained references.** Audit `lis_gesolve()`, inner solvers and
   transformed generalized algorithms to identify every A/B and derived-matrix
   reference retained after an eigen solve. Verify freeing only
   `esolver->work` is enough before using it; assume it is **not** until proven.
2. **Classify direct support.** For every relevant eigensolver type and
   spectral transformation, verify whether USER/OPERATOR matvec, matvech,
   diagonal and row access suffice. Do not copy the ordinary linear solver's
   direct-vs-CSR decision blindly.
3. **Track both sources.** Include recursively nested operator generations;
   reset incompatible state before cache mutation. Publication of a successful
   pair is atomic at the session boundary. Failed B preparation must not
   expose a prepared A as a successful pair.
4. **Own result lifetimes.** Document whether selected eigenpairs and stored
   eigenvectors remain readable after a subsequent solve. Ensure no use after
   free when rebuilding caches.
5. **Shared sparse structure.** Reuse a source shared topology where present.
   Sharing materialized CSR ptr/index arrays between K and B is a **separate**
   optimization requiring explicit reference-counted ownership, identical
   pattern verification and independently owned numerical values.
6. **Safety limits.** No concurrent session use or source mutation. No claim
   of multi-rank distributed cached CSR until implemented and tested.

## Proposed tests

- Generalized modal and buckling eigenproblems; compare selected eigenpairs
  and physical residuals with fresh, independent CSR baselines.
- K-only value refresh, B-only value refresh, both refreshed and both
  structure changes. Verify serial and OpenMP.
- USER and nested OPERATOR sources. Verify option/solver changes cause reset.
- Failure in K and failure in B row callbacks; neither scenario publishes a
  partially prepared pair; next solve recovers.
- Previously selected eigenvectors remain valid if their lifetime is
  guaranteed by the API; otherwise access after reset is prohibited and
  documented.
- Real/complex, 32-/64-bit LIS indices, ASan/UBSan; MPI one-rank only.
- Benchmark memory and time versus two independent materializations. No
  performance claim until measurements establish a benefit.

## Do not implement as a naive two-cache wrapper

A wrapper which simply calls `lis_matrix_materialize_csr_cached(A)` and
`lis_matrix_materialize_csr_cached(B)` before `lis_gesolve()` can invalidate
previously retained eigensolver references. Preserving the cache is not by
itself enough to make persistent generalized-eigensolver integration safe.
