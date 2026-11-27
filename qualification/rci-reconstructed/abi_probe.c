#include <lis.h>
#include <lis_abi.h>
#include <stdio.h>
#include <string.h>
int main(int argc,char **argv){
 LIS_ABI_INFO a={0};
 if(argc!=3||lis_abi_get_info(&a,sizeof a)!=LIS_ABI_SUCCESS)return 1;
 unsigned expected=strcmp(argv[1],"int64")==0?8:4;
 unsigned omp=strcmp(argv[2],"openmp")==0;
 printf("ABI int=%u scalar=%u real=%u flags=%u header_int=%zu expected=%u\n",a.lis_int_bytes,a.lis_scalar_bytes,a.lis_real_bytes,a.flags,sizeof(LIS_INT),expected);
 if(sizeof(LIS_INT)!=expected||a.lis_int_bytes!=expected||a.lis_scalar_bytes!=sizeof(LIS_SCALAR)||a.lis_real_bytes!=sizeof(LIS_REAL)||a.scalar_kind!=LIS_ABI_SCALAR_REAL||sizeof(LIS_SCALAR)!=8||sizeof(LIS_REAL)!=8||a.flags&LIS_ABI_FLAG_MPI||!!(a.flags&LIS_ABI_FLAG_OPENMP)!=omp)return 2;
 puts("ABI PASS");return 0;
}
