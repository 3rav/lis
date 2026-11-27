#ifndef LIS_RCI_NATIVE_FULL_INTERNAL_H
#define LIS_RCI_NATIVE_FULL_INTERNAL_H
#include "lis_rci_native_full.h"
enum { FULL_AGAIN=3,FULL_LINEAR=1,FULL_FEAST=2,FULL_CHILD=1,FULL_CERT=2 };
struct lis_rci_native_full {
 int mode,stage,pending,status,forced,threads;
 uint64_t token,child_token;
 LIS_INT completion,api_code,native_code,forced_code;
 LIS_MATRIX A,B;
 LIS_SOLVER solver;
 LIS_ESOLVER esolver;
 LIS_SCALAR *evalue;
 lis_rci_native_generations generations;
 lis_rci_native_session *child;
 lis_rci_native_session_request request;
 const void *cert_in[1];
 void *cert_out[1];
 void *frame;
};
int full_native_finish(lis_rci_native_full *,LIS_INT);
int full_native_begin_linear(lis_rci_native_full *,LIS_SOLVER,LIS_INT);
int full_native_begin_feast(lis_rci_native_full *,LIS_ESOLVER);
int full_native_cert(lis_rci_native_full *,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);
int full_native_feast_frame(lis_rci_native_full *,LIS_MATRIX,LIS_MATRIX,LIS_VECTOR,LIS_SCALAR *,LIS_ESOLVER);
int full_native_feast_step(lis_rci_native_full *);
#endif
