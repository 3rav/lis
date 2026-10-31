#include "lis_abi.h"
#include "lis.h"
#include <stdio.h>
#include <string.h>

int main(void)
{
    LIS_ABI_INFO info,old;
    int fail=0;
    memset(&info,0x5a,sizeof(info));
    old=info;
    if(lis_abi_get_info(NULL,sizeof(info))==0 ||
       lis_abi_get_info(&info,sizeof(info)-1)==0 ||
       memcmp(&old,&info,sizeof(info))!=0)
        fail=1;
    if(lis_abi_get_info(&info,sizeof(info))!=0 ||
       info.version!=LIS_ABI_INFO_VERSION ||
       info.struct_size!=sizeof(info) ||
       info.lis_int_bytes!=sizeof(LIS_INT) ||
       info.lis_scalar_bytes!=sizeof(LIS_SCALAR) ||
       info.lis_real_bytes!=sizeof(LIS_REAL) ||
       info.reserved!=0)
        fail=1;
#ifdef _COMPLEX
    if(info.scalar_kind!=LIS_ABI_SCALAR_COMPLEX) fail=1;
#else
    if(info.scalar_kind!=LIS_ABI_SCALAR_REAL) fail=1;
#endif
#ifdef USE_MPI
    if(!(info.flags & LIS_ABI_FLAG_MPI)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_MPI) fail=1;
#endif
#ifdef _OPENMP
    if(!(info.flags & LIS_ABI_FLAG_OPENMP)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_OPENMP) fail=1;
#endif
#ifdef USE_QUAD_PRECISION
    if(!(info.flags & LIS_ABI_FLAG_QUAD_SOLVER)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_QUAD_SOLVER) fail=1;
#endif
#ifdef _LONG__DOUBLE
    if(!(info.flags & LIS_ABI_FLAG_LONG_DOUBLE)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_LONG_DOUBLE) fail=1;
#endif
#if defined(_COMPLEX) && defined(_LONG__DOUBLE) && defined(HAVE_QUADMATH_H)
    if(!(info.flags & LIS_ABI_FLAG_FLOAT128)) fail=1;
#else
    if(info.flags & LIS_ABI_FLAG_FLOAT128) fail=1;
#endif
    printf("LIS_ABI_INFO=%s\n",fail?"FAIL":"PASS");
    return fail;
}
