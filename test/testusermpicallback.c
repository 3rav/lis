#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <stdio.h>
#include <math.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

/*
 * Distributed USER with APPLICATION-owned matvec and matvech callbacks.
 * This diagonal operator needs no halo data. This is NOT a test of
 * LIS-native distributed row callbacks or USER->CSR materialization.
 */
typedef struct {
    LIS_INT n;
    LIS_INT first_global;
    LIS_INT forward_calls;
    LIS_INT adjoint_calls;
} CONTEXT;

static LIS_SCALAR diagonal(LIS_INT global_i)
{
    return (LIS_SCALAR)(4.0+0.25*(double)global_i);
}
static LIS_INT ax(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CONTEXT *ctx=(CONTEXT *)ptr;
    LIS_INT i;
    ctx->forward_calls++;
    for(i=0;i<ctx->n;i++) y[i]=diagonal(ctx->first_global+i)*x[i];
    return LIS_SUCCESS;
}
static LIS_INT ahx(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CONTEXT *ctx=(CONTEXT *)ptr;
    LIS_INT i;
    ctx->adjoint_calls++;
    for(i=0;i<ctx->n;i++) y[i]=diagonal(ctx->first_global+i)*x[i];
    return LIS_SUCCESS;
}
static int solve_method(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,
                        LIS_SOLVER solver,const char *options)
{
    LIS_INT err,status,i;
    err=lis_solver_set_option((char *)options,solver);
    if(err)return 1;
    lis_vector_set_all((LIS_SCALAR)0.0,x);
    err=lis_solve(A,b,x,solver);
    if(err)return 1;
    err=lis_solver_get_status(solver,&status);
    if(err || status!=LIS_SUCCESS)return 1;
    for(i=0;i<A->n;i++)
    {
        LIS_INT g=A->is+i;
        LIS_SCALAR exact=(LIS_SCALAR)(1.0+0.125*(double)g);
        if(fabs((double)(x->value[i]-exact))>1e-7)return 1;
    }
    return 0;
}
int main(int argc,char **argv)
{
    CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL,b=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=LIS_SUCCESS,i,n=0,gn=0,is=0,ie=0;
    int failed=0,global_failed=0,rank=0,nranks=1;
#ifdef USE_MPI
    int mpi_err;
#endif
    err=lis_initialize(&argc,&argv);
    if(err)return 1;
#ifdef USE_MPI
    MPI_Comm_rank(LIS_COMM_WORLD,&rank);
    MPI_Comm_size(LIS_COMM_WORLD,&nranks);
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err){failed=1;goto done;}
    err=lis_matrix_set_size(A,0,16);if(err){failed=1;goto done;}
    err=lis_matrix_get_size(A,&n,&gn);if(err){failed=1;goto done;}
    err=lis_matrix_get_range(A,&is,&ie);if(err){failed=1;goto done;}
    if(gn!=16 || ie-is!=n){failed=1;goto done;}
    ctx.n=n;ctx.first_global=is;ctx.forward_calls=0;ctx.adjoint_calls=0;
    err=lis_matrix_set_user(A,&ctx,ax,ahx);if(err){failed=1;goto done;}
    err=lis_vector_duplicate(A,&x);if(err){failed=1;goto done;}
    err=lis_vector_duplicate(A,&y);if(err){failed=1;goto done;}
    err=lis_vector_duplicate(A,&b);if(err){failed=1;goto done;}
    for(i=0;i<n;i++)
    {
        LIS_INT g=is+i;
        x->value[i]=(LIS_SCALAR)(1.0+0.125*(double)g);
        b->value[i]=diagonal(g)*x->value[i];
    }
    err=lis_matvec(A,x,y);
    if(err){failed=1;goto done;}
    for(i=0;i<n;i++)
       if(fabs((double)(y->value[i]-b->value[i]))>1e-11)
       {failed=1;goto done;}
    err=lis_matvech(A,x,y);
    if(err){failed=1;goto done;}
    for(i=0;i<n;i++)
       if(fabs((double)(y->value[i]-b->value[i]))>1e-11)
       {failed=1;goto done;}
    err=lis_solver_create(&solver);if(err){failed=1;goto done;}
    if(solve_method(A,b,x,solver,
        "-i cg -p none -scale none -maxiter 200 -tol 1e-10 -print none"))
    {failed=1;goto done;}
    if(solve_method(A,b,x,solver,
        "-i bicg -p none -scale none -maxiter 200 -tol 1e-10 -print none"))
    {failed=1;goto done;}
    if(ctx.forward_calls<2 || ctx.adjoint_calls<2)
    {failed=1;goto done;}
done:
    if(solver)lis_solver_destroy(solver);
    if(b)lis_vector_destroy(b);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
#ifdef USE_MPI
    mpi_err=MPI_Allreduce(&failed,&global_failed,1,MPI_INT,MPI_MAX,LIS_COMM_WORLD);
    if(mpi_err!=MPI_SUCCESS)global_failed=1;
#else
    global_failed=failed;
#endif
    lis_finalize();
    if(rank==0)
       printf("USER_MPI_CALLBACK %s ranks=%d\n",
              global_failed?"FAILED":"PASSED",nranks);
    return global_failed?1:0;
}
