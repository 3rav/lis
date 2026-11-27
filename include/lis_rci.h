#ifndef LIS_RCI_RECONSTRUCTED_H
#define LIS_RCI_RECONSTRUCTED_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Independent reconstructed API; not the original LIS Generic RCI ABI. */
typedef struct lis_rci_context lis_rci_context;
typedef enum { LIS_RCI_CG=1, LIS_RCI_BICGSTAB, LIS_RCI_GMRES, LIS_RCI_FGMRES } lis_rci_method;
typedef enum { LIS_RCI_DONE=0, LIS_RCI_MATVEC=1, LIS_RCI_PRECON=2, LIS_RCI_MAXITER=-1,
 LIS_RCI_BREAKDOWN=-2, LIS_RCI_BADARG=-3, LIS_RCI_CALLBACK=-4, LIS_RCI_NONFINITE=-5,
 LIS_RCI_STATE=-6 } lis_rci_status;
typedef struct { size_t n, maxiter, restart; double rtol, atol; lis_rci_method method; } lis_rci_options;
typedef struct { const double *in; double *out; size_t n; } lis_rci_request;
typedef struct { size_t iterations, matvecs, precons; double residual; int status, callback_code; } lis_rci_info;
lis_rci_context *lis_rci_create(const lis_rci_options *options);
void lis_rci_destroy(lis_rci_context *c);
/* b and x0 are copied. x0 may be NULL (zero guess). Cannot start while active. */
int lis_rci_start(lis_rci_context *c,const double *b,const double *x0,lis_rci_request *q);
/* Complete previous request; callback_code=0 on success. Request buffers valid until next call. */
int lis_rci_advance(lis_rci_context *c,int callback_code,lis_rci_request *q);
const double *lis_rci_solution(const lis_rci_context *c);
int lis_rci_get_info(const lis_rci_context *c,lis_rci_info *info);
#ifdef __cplusplus
}
#endif
#endif
