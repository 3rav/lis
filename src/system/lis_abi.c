/* Fixed-layout runtime build ABI query for embedding applications. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif
#include <string.h>
#include "lislib.h"
#include "lis_abi.h"

int lis_abi_get_info(LIS_ABI_INFO *info, size_t info_size)
{
    LIS_ABI_INFO result;

    if(info==NULL || info_size<sizeof(result))
        return 1;

    memset(&result,0,sizeof(result));
    result.version=LIS_ABI_INFO_VERSION;
    result.struct_size=(uint32_t)sizeof(result);
    result.lis_int_bytes=(uint32_t)sizeof(LIS_INT);
    result.lis_scalar_bytes=(uint32_t)sizeof(LIS_SCALAR);
    result.lis_real_bytes=(uint32_t)sizeof(LIS_REAL);
#ifdef _COMPLEX
    result.scalar_kind=LIS_ABI_SCALAR_COMPLEX;
#else
    result.scalar_kind=LIS_ABI_SCALAR_REAL;
#endif
#ifdef USE_MPI
    result.flags|=LIS_ABI_FLAG_MPI;
#endif
#ifdef _OPENMP
    result.flags|=LIS_ABI_FLAG_OPENMP;
#endif
#ifdef USE_QUAD_PRECISION
    result.flags|=LIS_ABI_FLAG_QUAD_SOLVER;
#endif
#ifdef _LONG__DOUBLE
    result.flags|=LIS_ABI_FLAG_LONG_DOUBLE;
#endif
#if defined(_COMPLEX) && defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
    result.flags|=LIS_ABI_FLAG_FLOAT128;
#endif
    memcpy(info,&result,sizeof(result));
    return 0;
}
