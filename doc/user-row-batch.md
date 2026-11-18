# Optional batched row access for `LIS_MATRIX_USER`

This experimental, opt-in extension reduces per-row function-pointer dispatch in
USER matrix-vector operations and row-source consumers such as matrix
materialization and sparse preconditioner setup. It does **not** alter
ordinary CSR or the existing `LIS_MATRIX_USER_GET_ROW` contract.

## API

First initialize the USER matrix with `lis_matrix_set_user()` and install
the existing single-row callback:

```c
lis_matrix_set_user(A, user_data, NULL, NULL);
lis_matrix_set_user_row(A, get_one_row);
lis_matrix_set_user_rows_batch(A, get_multiple_rows);
lis_matrix_set_user_row_thread_safe(A, LIS_TRUE);
```

The batch callback signature is:

```c
LIS_INT get_multiple_rows(
    void *user_data, LIS_INT row_begin, LIS_INT row_end,
    LIS_INT *offsets, LIS_INT capacity,
    LIS_INT *cols, LIS_SCALAR *vals);
```

The row interval is half-open `[row_begin,row_end)`, using locally owned
row indices, while returned columns are zero-based **global** column IDs.
`offsets` contains `row_end-row_begin+1` entries, starts at zero, and
must be nondecreasing. The final offset is the number of nonzeros.

For the first, *count* call, `cols` and `vals` are `NULL` and
`capacity=0`. The callback must fill `offsets`. For the subsequent
*fill* call, `capacity` is at least the required number of nonzeros and the
callback must reproduce the same offsets while filling the arrays.

The matrix must remain unchanged between the two calls. LIS checks for
inconsistent offsets and out-of-range column indices. The application must
not retain any of the scratch pointers beyond the call.

## Performance and safety

- Each row-source instance caches up to 16 adjacent rows; the next request
  outside that range rebuilds the cache.
- USER matvec and matrix materialization create worker-private row-source
  state; cache buffers and offsets are never shared across workers.
- Only serial/local USER column numbering is supported by the row-source
  implementation (same as the preexisting USER direct-row path).
- A new batch callback resets the thread-safe flag. Call
  `lis_matrix_set_user_row_thread_safe(A,LIS_TRUE)` only if **both** installed
  callbacks and their application state are safe for concurrent calls.
- Explicit application-provided matvec still takes precedence over the
  generic callback-backed matvec.
- The single-row callback remains mandatory in this first version, for
  existing direct-row API compatibility and automatic fallback.
- Do not mutate matrix coefficients or callbacks while another thread is
  constructing or using a batch cache.

## Compatibility

This adds an optional function pointer at the end of the public
`LIS_MATRIX_STRUCT`. It is a source-compatible extension but **not a
promise of binary compatibility** for applications compiled with older LIS
headers. Rebuild applications and review SONAME/version policy before
release.

## Validation

`test/testuserbatch.c` compares scalar and batched matvec/matvech, checks
out-of-order row requests, validates malformed offset rejection, and counts
callback dispatches. The separate `user-batch-regression.yml` workflow
runs one and eight OpenMP threads with normal and ASan/UBSan builds.

The expected performance benefit depends on callback overhead and row
density. No speedup is claimed without representative application
benchmarks.
