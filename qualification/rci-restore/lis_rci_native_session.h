/* Candidate typed native Generic RCI session; not an installed ABI. */
#ifndef LIS_RCI_NATIVE_SESSION_H
#define LIS_RCI_NATIVE_SESSION_H
#include "lis_rci_restore.h"
#include "lis.h"
#include "lis_krylov.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct lis_rci_native_session lis_rci_native_session;
enum { LIS_NR_LINEAR=1,LIS_NR_GCRO=2,LIS_NR_FGCRO=3,LIS_NR_FEAST=4,LIS_NR_KRYLOV=5 };
enum { LIS_NR_VECTOR=1,LIS_NR_SPLIT=2,LIS_NR_RAW64=3 };
enum { LIS_NR_STALE=-9 };
/* Same operation codes as the opaque Generic RCI protocol. Vector arrays
 * are borrowed native handles; raw64 arrays are borrowed double buffers.
 * Shift sign explicitly identifies sigma B-A (+1) versus A-sigma B (-1).
 * Scalar/block dispatch is preserved even for a one-column block.
 * The prepared solver/eigensolver is exclusively borrowed; setup and
 * teardown are still caller responsibilities, not a public lis_solve API.
 * Context kind is fixed, tokens survive restarts, wrong tokens are retryable.
 * Output descriptors are unchanged on nonconsuming errors. */
typedef struct {
 uint64_t token;
 uint32_t operation,representation;
 LIS_INT count;
 int block,shift_sign;
 LIS_MATRIX A,B;
 LIS_SOLVER solver;
 LIS_VECTOR *in,*out;
 LIS_VECTOR in_real,in_imag,out_real,out_imag;
 LIS_REAL zr,zi;
 int64_t n64;
 const double *raw_in;
 double *raw_out;
 double sigma;
} lis_rci_native_session_request;
lis_rci_native_session *lis_rci_native_session_create(int kind);
void lis_rci_native_session_destroy(lis_rci_native_session *);
int lis_rci_native_session_start_linear(lis_rci_native_session *,LIS_SOLVER,LIS_INT method,lis_rci_native_session_request *);
int lis_rci_native_session_start_feast(lis_rci_native_session *,LIS_ESOLVER,lis_rci_native_session_request *);
int lis_rci_native_session_start_krylov(lis_rci_native_session *,const lis_krylov_problem64 *,const lis_krylov_options *,const lis_krylov_ext_options *,lis_krylov_workspace *,const lis_krylov_budget *,double *,double *,double *,lis_krylov_stats *,lis_krylov_ext_stats *,lis_rci_native_session_request *);
int lis_rci_native_session_advance(lis_rci_native_session *,uint64_t,LIS_INT,lis_rci_native_session_request *);
int lis_rci_native_session_cancel(lis_rci_native_session *,uint64_t);
int lis_rci_native_session_pending(const lis_rci_native_session *,lis_rci_native_session_request *);
int lis_rci_native_session_status(const lis_rci_native_session *,LIS_INT *native,uint64_t *a,uint64_t *p,lis_krylov_budget_stats *);
#ifdef __cplusplus
}
#endif
#endif
