#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

#if !defined(USE_MPI) && !defined(_COMPLEX)
typedef struct { LIS_INT n,calls,fail; LIS_SCALAR diag; } CTX;
static LIS_INT apply(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)ptr;
    LIS_INT i;
    ++c->calls;
    if(c->fail && c->calls==c->fail)
    {
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<c->n;i++) y[i]=c->diag*x[i];
    return LIS_SUCCESS;
}
static LIS_INT create(CTX *c,LIS_MATRIX *A)
{
    LIS_INT err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err)return err;
    err=lis_matrix_set_size(*A,0,c->n);
    if(err)return err;
    return lis_matrix_set_user(*A,c,apply,apply);
}
static int run(int buckling,LIS_INT afail,LIS_INT bfail)
{
    CTX ac={4,0,afail,4.0},bc={4,0,bfail,2.0};
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR value=0;
    LIS_REAL residual=0;
    LIS_INT err,status,i,iter;
    int failed=1;
    if(create(&ac,&A))goto done;
    if(create(&bc,&B))goto done;
    if(lis_vector_duplicate(A,&x) || lis_esolver_create(&e))goto done;
    if(lis_esolver_set_option(buckling ?
        "-e gii -etransform buckling -shift 1 -etol 1e-10 -emaxiter 50 -eprint none" :
        "-e gii -etransform shift-invert -etol 1e-10 -emaxiter 50 -eprint none",e))goto done;
    if(lis_vector_set_all(1.0,x))goto done;
    err=lis_gesolve(A,B,x,&value,e);
    if(lis_esolver_get_status(e,&status) || lis_esolver_get_iter(e,&iter) ||
       lis_esolver_get_residualnorm(e,&residual))goto done;
    if(afail || bfail)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG || iter!=0 ||
           !isnan((double)residual) ||
           (afail && ac.calls!=afail) || (bfail && bc.calls<bfail))
        {
            fprintf(stderr,"fault gen=%d A=%d B=%d api=%d status=%d iter=%d calls=%d/%d resid=%g\n",
                buckling,(int)afail,(int)bfail,(int)err,(int)status,(int)iter,
                (int)ac.calls,(int)bc.calls,(double)residual);
            goto done;
        }
    }
    else if(err!=LIS_SUCCESS || status!=LIS_SUCCESS ||
            !(fabs((double)value-2.0)<=1e-9) ||
            !(residual<=1e-10))goto done;
    for(i=0;i<x->n;i++)if(!isfinite((double)x->value[i]))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"USER eigen inverse failure gen=%d A=%d B=%d\n",
                      buckling,(int)afail,(int)bfail);
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    return failed;
}
#endif
int main(int argc,char **argv)
{
#if defined(USE_MPI) || defined(_COMPLEX)
    (void)argc;(void)argv;
    puts("USER eigen inverse SKIP: real serial regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(0,0,1);
    failed|=run(1,1,0);
    failed|=run(0,0,0);
    failed|=run(1,0,0);
    lis_finalize();
    printf("USER_EIGEN_INVERSE %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
