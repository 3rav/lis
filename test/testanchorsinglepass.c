/* Native, serial real-double prototype regression. */
#include "lislib.h"
#include "lis_solver_options.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static LIS_INT matrix(int kind,LIS_MATRIX *A)
{
    LIS_INT i,rc;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,A)) || (rc=lis_matrix_set_size(*A,0,80)))return rc;
    for(i=0;i<80;i++){
        double d=kind==0?1:kind==1?(double)(i+1):pow(10.0,-8.0+8.0*i/79.0);
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,*A)))return rc;
    }
    if((rc=lis_matrix_set_type(*A,LIS_MATRIX_CSR)))return rc;
    return lis_matrix_assemble(*A);
}
int main(int argc,char **argv)
{
    LIS_MATRIX P=NULL,A=NULL,B=NULL;LIS_VECTOR b=NULL,x=NULL,y=NULL;
    lis_anchor_pilot *ctx=NULL;lis_anchor_pilot_report r1,r2;
    int fail=1;
    LIS_INT q;
    if(lis_initialize(&argc,&argv))return 1;
    if(matrix(0,&P)||matrix(1,&A)||matrix(2,&B)||lis_vector_duplicate(P,&b)||
       lis_vector_duplicate(P,&x)||lis_vector_duplicate(P,&y)||lis_vector_set_all(1,b))goto done;
    if(lis_anchor_pilot_create(P,&ctx)||lis_anchor_pilot_solve(ctx,A,b,x,&r1))goto done;
    lis_anchor_pilot_destroy(ctx);ctx=NULL;
    if(lis_anchor_pilot_create_single_pass(P,&ctx)||lis_anchor_pilot_solve(ctx,A,b,y,&r2))goto done;
    printf("CONTINUATION old=(%lld+%lld) new=(%lld+%lld) residual=%.15e / %.15e\n",
        (long long)r1.pilot_iterations,(long long)r1.continuation_iterations,
        (long long)r2.pilot_iterations,(long long)r2.continuation_iterations,r1.relative_residual,r2.relative_residual);
    if(!r1.accepted||!r1.continued||r1.fallback||!r2.accepted||!r2.continued||r2.fallback||r2.factorizations||
       r2.pilot_iterations!=20||!(r2.relative_residual<=1e-8))goto done;
    if(lis_anchor_pilot_solve(ctx,B,b,x,&r2)||!r2.accepted||!r2.fallback||r2.continued||r2.pilot_iterations!=20||r2.factorizations!=1)goto done;
    if(lis_anchor_pilot_solve(ctx,P,b,x,&r2)||!r2.accepted||r2.fallback||r2.factorizations||r2.pilot_iterations!=1)goto done;
    lis_anchor_pilot_destroy(ctx);ctx=NULL;
    if(lis_anchor_pilot_create_restart_matched(P,&ctx)||lis_anchor_pilot_solve(ctx,A,b,x,&r1))goto done;
    printf("MATCHED_DRIVER pilot=%lld continuation=%lld residual=%.15e\n",
        (long long)r1.pilot_iterations,(long long)r1.continuation_iterations,r1.relative_residual);
    if(!r1.accepted||!r1.continued||r1.fallback||r1.factorizations||r1.pilot_iterations!=20||
       r1.continuation_iterations<1||r1.continuation_iterations>80||!(r1.relative_residual<=1e-8))goto done;
    /* The next RHS must begin with the original tolerance, not the relaxed one. */
    if(lis_anchor_pilot_solve(ctx,A,b,y,&r2)||!r2.accepted||r2.fallback||r2.factorizations||
       r1.pilot_iterations!=r2.pilot_iterations||r1.continuation_iterations!=r2.continuation_iterations||
       r1.relative_residual!=r2.relative_residual||memcmp(x->value,y->value,(size_t)x->n*sizeof(LIS_SCALAR)))goto done;
    if(lis_anchor_pilot_solve(ctx,B,b,x,&r2)||!r2.accepted||!r2.fallback||r2.continued||r2.pilot_iterations!=20||r2.factorizations!=1)goto done;
    if(lis_anchor_pilot_solve(ctx,P,b,x,&r2)||!r2.accepted||r2.fallback||r2.factorizations||r2.pilot_iterations!=1)goto done;
    /* Zero RHS must accept without dividing by zero, no fallback. */
    if(lis_vector_set_all(0,b)||lis_vector_set_all(7,x))goto done;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&r2)||!r2.accepted||r2.fallback||r2.relative_residual!=0)goto done;
    for(q=0;q<x->n;q++)if(x->value[q]!=0)goto done;
    /* Nonfinite RHS and matrix values are rejected before factor/solve;
     * output remains untouched and the context remains reusable. */
    if(lis_vector_set_all(1,b)||lis_vector_set_all(7,x))goto done;
    b->value[0]=NAN;
    if(lis_anchor_pilot_solve(ctx,A,b,x,&r2)!=LIS_ERR_ILL_ARG||r2.accepted||r2.factorizations)goto done;
    b->value[0]=INFINITY;
    if(lis_anchor_fresh_solve(A,b,x,&r2)!=LIS_ERR_ILL_ARG||r2.accepted||r2.factorizations)goto done;
    b->value[0]=1;
    {
        double old=A->value[0];
        A->value[0]=NAN;
        if(lis_anchor_pilot_solve(ctx,A,b,x,&r2)!=LIS_ERR_ILL_ARG||r2.accepted||r2.factorizations)goto done;
        A->value[0]=old;
        if(lis_matrix_mark_values_changed(A))goto done;
    }
#ifdef _OPENMP
    {
        LIS_INT rc=LIS_SUCCESS;
        #pragma omp parallel num_threads(2)
        {
            #pragma omp single
            {
                omp_set_num_threads(1);
                rc=lis_anchor_pilot_solve(ctx,A,b,x,&r2);
            }
        }
        if(rc!=LIS_ERR_ILL_ARG||r2.accepted||r2.factorizations)goto done;
    }
#endif
    for(q=0;q<x->n;q++)if(x->value[q]!=7)goto done;
    if(lis_anchor_pilot_solve(ctx,P,b,x,&r2)||!r2.accepted||r2.fallback)goto done;
    fail=0;
done:
    lis_anchor_pilot_destroy(ctx);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    if(P)lis_matrix_destroy(P);
    lis_finalize();
    printf("SINGLE_PASS_REGRESSION=%s\n",fail?"FAIL":"PASS");
    return fail;
}
