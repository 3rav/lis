#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct { LIS_INT n,calls,fail_on; } USER_FAULT;
static LIS_INT apply(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    USER_FAULT *ctx=(USER_FAULT *)raw;
    LIS_INT i;
    ++ctx->calls;
    if(ctx->fail_on>0 && ctx->calls==ctx->fail_on)
        return LIS_ERR_ILL_ARG;
    for(i=0;i<ctx->n;i++)
        y[i]=(LIS_SCALAR)(4.0+0.2*(double)i)*x[i];
    return LIS_SUCCESS;
}

static int run(LIS_INT fail_on)
{
    USER_FAULT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,b=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,err=LIS_SUCCESS,status=LIS_SUCCESS;
    int failed=1;
    ctx.n=8;ctx.calls=0;ctx.fail_on=fail_on;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,apply,apply);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.1;
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(
        "-i tfqmr -p none -scale none -maxiter 200 -tol 1e-9 -print none",
        solver);
    if(err)goto done;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status))goto done;
    if(fail_on)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG || ctx.calls<fail_on)
        {
            fprintf(stderr,"TFQMR fault %d api=%d status=%d calls=%d\n",
                   (int)fail_on,(int)err,(int)status,(int)ctx.calls);
            goto done;
        }
    }
    else
    {
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)
        {
            fprintf(stderr,"TFQMR healthy api=%d status=%d calls=%d\n",
                    (int)err,(int)status,(int)ctx.calls);
            goto done;
        }
        for(i=0;i<ctx.n;i++)
        {
            double exact=(1.0+0.125*(double)i)/(4.0+0.2*(double)i);
            if(fabs((double)x->value[i]-exact)>1e-7)goto done;
        }
    }
    failed=0;
done:
    if(failed)fprintf(stderr,"USER_TFQMR FAIL inject=%d api=%d\n",
                      (int)fail_on,(int)err);
    if(solver)lis_solver_destroy(solver);
    if(b)lis_vector_destroy(b);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(2);
    failed|=run(3);
    failed|=run(4);
    failed|=run(5);
    failed|=run(0);
    lis_finalize();
    printf("USER_TFQMR %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
}
