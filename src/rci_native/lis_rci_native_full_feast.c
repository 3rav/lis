/* Native FEAST setup/finalization remains in its own source translation unit. */
#include "lis_rci_native_full_internal.h"
#include "full_feast_continuation.inc"
int full_native_feast_frame(lis_rci_native_full *S,LIS_MATRIX A,LIS_MATRIX B,LIS_VECTOR x,LIS_SCALAR *value,LIS_ESOLVER e)
{
 struct full_frames_feast *C=calloc(1,sizeof(*C));if(!C)return LIS_ERR_OUT_OF_MEMORY;
 C->eigen.A=A;C->eigen.B=B;C->eigen.x=x;C->eigen.evalue0=value;C->eigen.esolver=e;
 S->frame=C;S->stage=1;
 return LIS_SUCCESS;
}
int full_native_feast_step(lis_rci_native_full *S)
{
 /* Scope belongs to this C call only; never leave native TLS pointing
  * at a suspended frame. Explicit external callback errors use completion. */
 LIS_CALLBACK_SCOPE scope;lis_callback_scope_begin(&scope);
 int st=full_step_eigen(S);lis_callback_scope_end(&scope);return st;
}
