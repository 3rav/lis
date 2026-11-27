/* Experimental full-lifecycle native solve over versioned RCI requests.
 * AI-assisted reconstruction under proposed LIS COPYING terms. */
#ifndef LIS_RCI_NATIVE_FULL_H
#define LIS_RCI_NATIVE_FULL_H
#include "lis_rci_native_session.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_full lis_rci_native_full;
/* Full setup, recurrence and cleanup use pinned native LIS source.
 * All operands and solver objects remain live and exclusive until terminal.
 * Current qualification: real double, no MPI; linear CSR or USER with
 * NONE/Jacobi precon (USER requires NONE), no scaling/storage conversion,
 * nullspace correction, or per-solver thread override. FEAST CSR A/B,
 * default mode and no automatic subspace retry. Unsupported configurations
 * are rejected before native setup. Keep the same host thread and inherited
 * OpenMP thread count for all calls. Native preconditioner construction
 * remains synchronous; operator and preconditioner *applications* in the
 * recurrence and final true-residual matvec are RCI requests.
 * Numerical failure and public API/dispatch return are reported separately.
 * Cancellation unwinds native cleanup and does not guarantee a valid iterate.
 */
lis_rci_native_full *lis_rci_native_full_create(void);
void lis_rci_native_full_destroy(lis_rci_native_full *);
int lis_rci_native_full_start_linear(lis_rci_native_full *,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER,
 lis_rci_native_session_request *,size_t);
int lis_rci_native_full_start_feast(lis_rci_native_full *,LIS_MATRIX,LIS_MATRIX,LIS_VECTOR,LIS_SCALAR *,LIS_ESOLVER,
 lis_rci_native_session_request *,size_t);
int lis_rci_native_full_advance(lis_rci_native_full *,uint64_t,LIS_INT,lis_rci_native_session_request *,size_t);
int lis_rci_native_full_cancel(lis_rci_native_full *,uint64_t);
int lis_rci_native_full_status(const lis_rci_native_full *,LIS_INT *api_code,LIS_INT *native_code);
#ifdef __cplusplus
}
#endif
#endif
