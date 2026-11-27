/* Experimental versioned native session API, optional installed.
 * AI-assisted reconstruction, proposed under LIS COPYING terms. */
#ifndef LIS_RCI_NATIVE_SESSION_H
#define LIS_RCI_NATIVE_SESSION_H
#include "lis_rci_native_protocol.h"
#include "lis.h"
#include "lis_krylov.h"
#ifdef __cplusplus
extern "C" {
#endif
#define LIS_RCI_NATIVE_SESSION_VERSION 1u
enum { LIS_RCI_NATIVE_VECTORS=1u, LIS_RCI_NATIVE_RAW_DOUBLE=2u,
 LIS_RCI_NATIVE_BLOCK=4u, LIS_RCI_NATIVE_SPLIT=8u,
 LIS_RCI_NATIVE_SHIFT_A_MINUS_SIGMA_B=16u };
enum { LIS_RCI_NATIVE_STALE=-20 };
typedef struct { uint64_t a_structure,a_values,b_structure,b_values; } lis_rci_native_generations;
typedef struct {
 uint32_t version,size,operation,flags;
 uint32_t index_bytes,scalar_bytes,real_bytes,reserved;
 uint64_t token,rows,columns;
 lis_rci_native_generations generations;
 const void *const *in;
 void *const *out;
 LIS_REAL shift_real,shift_imag;
 LIS_MATRIX A,B;
 LIS_SOLVER solver;
} lis_rci_native_session_request;
typedef struct lis_rci_native_session lis_rci_native_session;
/* All objects borrowed exclusively until terminal; descriptor arrays and
 * handles borrowed until next accepted reply/cancel/destroy. Inputs read-only.
 * Split FEAST requests have one logical RHS but two handle entries, real
 * followed by imaginary; solve (sigma B-A). Krylov has raw double vectors
 * and solve (A-sigma B), explicitly flagged. Native block requests preserve
 * block dispatch even for one column. Native scalar type is the LIS ABI.
 * Prepared native methods require the normal LIS setup/work/preconditioner.
 * This API does NOT prepare or finalize a full public lis_solve/lis_esolve.
 * Native generation counters checked automatically. Raw Krylov caller must
 * supply current generations on advance. Cannot safely detect destruction
 * or unsynchronized concurrent mutation of borrowed objects.
 * Request bytes must cover this version's complete structure; wrong size,
 * wrong token, or raw missing generations do not consume outstanding work.
 */
lis_rci_native_session *lis_rci_native_session_create(void);
void lis_rci_native_session_destroy(lis_rci_native_session *);
int lis_rci_native_session_start_linear(lis_rci_native_session *,LIS_SOLVER,LIS_INT method,
 lis_rci_native_session_request *,size_t);
int lis_rci_native_session_start_feast(lis_rci_native_session *,LIS_ESOLVER,
 lis_rci_native_session_request *,size_t);
int lis_rci_native_session_start_krylov(lis_rci_native_session *,const lis_krylov_problem64 *,
 const lis_krylov_options *,const lis_krylov_ext_options *,lis_krylov_workspace *,
 const lis_krylov_budget *,double *,double *,double *,lis_krylov_stats *,lis_krylov_ext_stats *,
 const lis_rci_native_generations *,lis_rci_native_session_request *,size_t);
int lis_rci_native_session_advance(lis_rci_native_session *,uint64_t,LIS_INT,
 const lis_rci_native_generations *,lis_rci_native_session_request *,size_t);
int lis_rci_native_session_cancel(lis_rci_native_session *,uint64_t);
int lis_rci_native_session_pending(const lis_rci_native_session *,lis_rci_native_session_request *,size_t);
int lis_rci_native_session_status(const lis_rci_native_session *,LIS_INT *native_code);
#ifdef __cplusplus
}
#endif
#endif
