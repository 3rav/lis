/* Optional CalculiX 2.19 -> generic LIS_MATRIX_USER bridge.
 * This is an out-of-tree validation fixture, not part of the ccx build.
 * Requires the unmodified CalculiX.h and an initialized LIS library.
 */
#ifndef CCX_LIS_USER_ADAPTER_H
#define CCX_LIS_USER_ADAPTER_H
#include <stddef.h>
#include <stdint.h>
#include "CalculiX.h"
#include "lis.h"
#include "lis_matrix_materialize.h"

typedef struct CCX_LIS_USER_CONTEXT CCX_LIS_USER_CONTEXT;

/* icol counts strictly-lower entries by column; irow uses origin 1.
 * ad/au are borrowed double arrays, and au_capacity is their actual
 * allocation length (number of doubles). Asymmetric upper values start
 * at au[nzs3], not necessarily au[nnz].
 *
 * The context retains the option string only through the LIS solver, and
 * must be closed before lis_finalize(). All calls are serial only.
 */
LIS_INT ccx_lis_user_open(
    ITG n, ITG nnz, const ITG icol[], const ITG irow[],
    const double ad[], const double au[], size_t au_capacity,
    ITG nzs3, int asymmetric, const char *options,
    CCX_LIS_USER_CONTEXT **out);

/* Notify LIS even when the borrowed values were updated in place. */
LIS_INT ccx_lis_user_update(
    CCX_LIS_USER_CONTEXT *ctx, const double ad[], const double au[],
    size_t au_capacity, ITG nzs3);

/* Transactional topology/borrowed-array replacement, fixed n. The
 * original matrix/session remains valid on invalid-input/allocation
 * failure. The session resets reusable preconditioners before releasing
 * the old matrix on a successful change of structure. */
LIS_INT ccx_lis_user_rebind(
    CCX_LIS_USER_CONTEXT *ctx, ITG nnz, const ITG icol[],
    const ITG irow[], const double ad[], const double au[],
    size_t au_capacity, ITG nzs3);

/* The supplied RHS is not modified. On solver failure solution is not
 * published. policy/action may be NULL; relative_residual may be NULL. */
LIS_INT ccx_lis_user_solve(
    CCX_LIS_USER_CONTEXT *ctx, const double rhs[], double solution[],
    double *relative_residual, LIS_INT *policy, LIS_INT *action);

void ccx_lis_user_close(CCX_LIS_USER_CONTEXT *ctx);
#endif
