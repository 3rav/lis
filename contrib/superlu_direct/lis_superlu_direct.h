#ifndef LIS_SUPERLU_DIRECT_H
#define LIS_SUPERLU_DIRECT_H
#include "lis.h"
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Register optional SuperLU full LU as a user-defined LIS preconditioner.
 * Call once per lis_initialize/lis_finalize session, before -p superlu.
 * Requires serial real-double CSR and explicit libsuperlu linking.
 */
LIS_INT lis_superlu_direct_register(void);
#ifdef __cplusplus
}
#endif
#endif
