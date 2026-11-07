/* Optional, solver-independent Matrix Market snapshot exporter.
 * No FEA model, mesh, DOF numbering or contact assumptions are made.
 * This advanced API is explicitly invoked by the embedding application.
 */
#ifndef LIS_SNAPSHOT_H
#define LIS_SNAPSHOT_H

#include "lis.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Export the exact current assembled linear operator and optional vectors.
 *
 * Creates <prefix>_A.mtx, optionally _b.mtx, _x0.mtx, _x.mtx and _B.mtx.
 * B is an n-by-basis_count Matrix Market array (column-major).
 *
 * Serial real-valued matrices only in this first version. Explicit CSR,
 * other convertible explicit formats, and row-backed USER/OPERATOR sources
 * are supported; unsupported/matrix-free sources without row access fail.
 * This does NOT solve or modify the operator, run FEA, or create a DOF map.
 *
 * All arguments are borrowed. A failed call returns nonzero; files from
 * earlier successful writes may remain, so discard the entire prefix set
 * on error. Use a unique prefix for each Newton/nonlinear snapshot.
 */
LIS_INT lis_snapshot_export_mm(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x0,
    LIS_VECTOR x,
    LIS_INT basis_count,
    LIS_VECTOR basis[],
    const char *prefix);

#ifdef __cplusplus
}
#endif

#endif
