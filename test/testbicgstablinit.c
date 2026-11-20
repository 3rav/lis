#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

/* Regression: BiCGSTAB(l) once added the nonzero initial x twice:
 * for A=4I, b[0]=1, x0[0]=0.15 it returned x[0]=0.40 instead of 0.25.
 *
 * The legacy method may report LIS_BREAKDOWN for an exact Krylov
 * termination, so certify Ax=b independently rather than relabel
 * that status as success.
 */
typedef struct { LIS_INT n; } OPERATOR;
static LIS_INT ax(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    OPERATOR *o=(OPERATOR *)ptr;
    LIS_INT i;
    for(i=0;i<o->n;i++) y[i]=(LIS_SCALAR)4.0*x[i];
    return LIS_SUCCESS;
}

static int run(LIS_INT ell,int init_zeros)
{
    OPERATOR op;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=LIS_SUCCESS,status=LIS_SUCCESS,i;
    int failed=1;
    op.n=8;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,op.n);if(err)goto done;
    err=lis_matrix_set_user(A,&op,ax,ax);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<op.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.15;
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(
        "-i bicgstabl -p none -scale none -maxiter 120 -tol 1e-9 -print none",
        solver);
    if(err)goto done;
    solver->options[LIS_OPTIONS_ELL]=ell;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=init_zeros?LIS_TRUE:LIS_FALSE;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status))goto done;
    /* The solver's true-residual post-certification must agree with
     * the numerical solution once x0 is added exactly once. */
    if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)
    {
        fprintf(stderr,"BiCGSTAB(l) api=%d status=%d ell=%d init_zeros=%d\n",
                (int)err,(int)status,(int)ell,init_zeros);
        goto done;
    }
    for(i=0;i<op.n;i++)
    {
        double exact=(1.0+0.125*(double)i)/4.0;
        if(fabs((double)x->value[i]-exact)>1e-8)
        {
            fprintf(stderr,"BiCGSTAB(l) x[%d]=%.17g exact=%.17g ell=%d init_zeros=%d status=%d\n",
                    (int)i,(double)x->value[i],exact,
                    (int)ell,init_zeros,(int)status);
            goto done;
        }
    }
    failed=0;
done:
    if(failed)fprintf(stderr,"BICGSTABL_NONZERO_INIT FAILED ell=%d zero=%d\n",
                      (int)ell,init_zeros);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(1,0);
    failed|=run(2,0);
    failed|=run(1,1);
    failed|=run(2,1);
    lis_finalize();
    printf("BICGSTABL_NONZERO_INIT %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
}
