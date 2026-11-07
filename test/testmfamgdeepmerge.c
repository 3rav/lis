#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <stdio.h>
#include "lislib.h"
#define N 121
extern LIS_INT lis_precon_mfamg_merge_siblings_for_test(
    LIS_MATRIX A, LIS_INT target, LIS_INT *aggregate, LIS_INT *nc);
int main(int argc, char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_INT id[N],sizes[N],nc=N,err=0,i;
    int failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1) { lis_finalize(); return 0; }
    }
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,N);
    if(err) goto done;
    for(i=0;i<N;i++) {
        id[i]=i; sizes[i]=0;
        err=lis_matrix_set_value(
            LIS_INS_VALUE,i,i,(i==0 ? 121.0 : 2.0),A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,0,i,-1.0,A);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,i,0,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_precon_mfamg_merge_siblings_for_test(A,96,id,&nc);
    if(err) goto done;
    if(nc!=3) {
        fprintf(stderr,"MFAMG DEEP MERGE expected 3 aggs, got %d\n",(int)nc);
        goto done;
    }
    for(i=0;i<N;i++) {
        if(id[i]<0 || id[i]>=nc) goto done;
        if(++sizes[id[i]]>96) goto done;
    }
    printf("MFAMG DEEP MERGE star test PASSED n=%d aggs=%d\n",N,(int)nc);
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG DEEP MERGE LIS error=%d\n",(int)err);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
