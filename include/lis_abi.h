/* ABI-stable runtime description of the linked LIS library.
 *
 * This header intentionally does NOT include lis.h.  A client may call
 * lis_abi_get_info() before making any calls that use LIS_INT/LIS_SCALAR,
 * even if its generated lis_config.h differs from the linked library.
 */
#ifndef LIS_ABI_H
#define LIS_ABI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define LIS_ABI_INFO_VERSION 1u
#define LIS_ABI_SCALAR_REAL 0u
#define LIS_ABI_SCALAR_COMPLEX 1u

#define LIS_ABI_FLAG_MPI          UINT32_C(0x00000001)
#define LIS_ABI_FLAG_OPENMP       UINT32_C(0x00000002)
#define LIS_ABI_FLAG_QUAD_SOLVER  UINT32_C(0x00000004)
#define LIS_ABI_FLAG_LONG_DOUBLE  UINT32_C(0x00000008)
#define LIS_ABI_FLAG_FLOAT128     UINT32_C(0x00000010)

typedef struct LIS_ABI_INFO_STRUCT
{
    uint32_t version;
    uint32_t struct_size;
    uint32_t lis_int_bytes;
    uint32_t lis_scalar_bytes;
    uint32_t lis_real_bytes;
    uint32_t scalar_kind;
    uint32_t flags;
    uint32_t reserved;
} LIS_ABI_INFO;

/*
 * Returns 0 on success, nonzero on invalid args or an undersized buffer.
 * No LIS-defined typedef is present in this signature or structure.
 * The caller must pass at least sizeof(LIS_ABI_INFO) bytes.
 * On invalid arguments the destination is not modified.
 *
 * This is a runtime compatibility query, NOT an ABI translation layer.
 * Embedders must reject incompatible linked libraries before using any
 * interface passing LIS_INT/LIS_SCALAR or LIS-owned structures.
 */
int lis_abi_get_info(LIS_ABI_INFO *info, size_t info_size);

#ifdef __cplusplus
}
#endif
#endif
