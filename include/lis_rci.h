/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#ifndef LIS_RCI_RECONSTRUCTED_H
#define LIS_RCI_RECONSTRUCTED_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Experimental reverse communication for real double vectors.
 * A context is not reentrant. Use distinct contexts and request buffers.
 * Operator/preconditioner must stay fixed and linear except that FGMRES
 * permits variable preconditioning; CG requires SPD operator and preconditioner.
 * GMRES and FGMRES share the right-preconditioned Arnoldi implementation.
 * The caller must complete each request exactly once, filling every output
 * entry, without modifying the input. Identity is a valid PRECON response.
 * DONE means the verified absolute residual satisfies max(atol,rtol*||b||).
 * info.residual is the last computed true absolute residual (or NaN);
 * on interruption it need not describe the most recent iterate.
 * Allocation failure or invalid options returns NULL from create.
 * No API/ABI stability promise is made for this experimental interface. */
typedef struct lis_rci_context lis_rci_context;
typedef enum { LIS_RCI_CG=1, LIS_RCI_BICGSTAB, LIS_RCI_GMRES, LIS_RCI_FGMRES } lis_rci_method;
typedef enum { LIS_RCI_DONE=0, LIS_RCI_MATVEC=1, LIS_RCI_PRECON=2, LIS_RCI_DOT=3, LIS_RCI_NORM=4, LIS_RCI_MAXITER=-1,
 LIS_RCI_BREAKDOWN=-2, LIS_RCI_BADARG=-3, LIS_RCI_CALLBACK=-4, LIS_RCI_NONFINITE=-5,
 LIS_RCI_STATE=-6, LIS_RCI_PROTOCOL=-7 } lis_rci_status;
typedef struct { size_t n, maxiter, restart; double rtol, atol; lis_rci_method method; } lis_rci_options;
typedef struct { const double *in; double *out; size_t n; } lis_rci_request;
typedef struct { size_t iterations, matvecs, precons; double residual; int status, callback_code; } lis_rci_info;
lis_rci_context *lis_rci_create(const lis_rci_options *options);
void lis_rci_destroy(lis_rci_context *c);
/* b and x0 are copied. x0 may be NULL (zero guess). Cannot start while active. */
int lis_rci_start(lis_rci_context *c,const double *b,const double *x0,lis_rci_request *q);
/* Complete previous request; callback_code=0 on success. Request buffers valid until next call. */
int lis_rci_advance(lis_rci_context *c,int callback_code,lis_rci_request *q);
/* Borrowed context-owned array; valid until destroy, contents may change on calls.
 * Only DONE guarantees convergence. Before start, contains zeros. */
const double *lis_rci_solution(const lis_rci_context *c);
int lis_rci_get_info(const lis_rci_context *c,lis_rci_info *info);
/* Extended RCI: all dot products and Euclidean norms are external requests.
 * MATVEC/PRECON use in/out; DOT uses in/in2; NORM uses in only.
 * Scalar requests have out == NULL; vector requests have in2 == NULL.
 * The scalar completion argument is ignored for vector requests.
 * Tokens are monotonic within one context, including across solves; never
 * compare tokens of unrelated contexts. Wrong tokens return PROTOCOL without
 * consuming a request. Do not write through expired output pointers.
 * Legacy and extended entry points must not be mixed within a solve.
 * All vectors are local full-length vectors; this is not a distributed solver. */
typedef struct {
 uint64_t token;
 int operation;
 size_t n;
 const double *in,*in2;
 double *out;
} lis_rci_request_ex;
int lis_rci_start_ex(lis_rci_context *c,const double *b,const double *x0,lis_rci_request_ex *q);
int lis_rci_advance_ex(lis_rci_context *c,uint64_t token,int callback_code,long double scalar,lis_rci_request_ex *q);
/* Optional transport-independent status consensus. One record per unique rank.
 * token and tag are opaque uint64_t values compared bitwise (not doubles).
 * Merge is componentwise unsigned MAX. Error 0 means success; the lowest
 * rank with a nonzero signed 32-bit error wins. Disagreements and the selected
 * error are independent; the caller defines precedence. No MPI calls. */
typedef struct {
 int token_mismatch,tag_mismatch;
 uint32_t error_rank; /* UINT32_MAX when no error; inspect error to disambiguate */
 int32_t error;
} lis_rci_consensus_result;
int lis_rci_consensus_pack(uint64_t token,uint64_t tag,uint32_t rank,int32_t error,uint64_t words[5]);
int lis_rci_consensus_merge(uint64_t dst[5],const uint64_t src[5]);
int lis_rci_consensus_unpack(const uint64_t words[5],lis_rci_consensus_result *result);
#ifdef __cplusplus
}
#endif
#endif
