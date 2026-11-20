#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <stdio.h>
#include <math.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

/* Cyclic off-rank couplings with duplicate remote columns. Dynamic row
 * changes require a fresh exchange map at every call. */
typedef struct { LIS_INT global_n,is,local_n; int rank,fail_rank,fail; } ROWS;
static LIS_INT get_row(void *data,LIS_INT row,LIS_INT *nnz,
                        LIS_INT *cols,LIS_SCALAR *vals)
{
    ROWS *s=(ROWS *)data;
    LIS_INT g;
    if(row<0 || row>=s->local_n || !nnz)return LIS_ERR_ILL_ARG;
    if(s->fail==1 && s->rank==s->fail_rank && row==0)
        return LIS_ERR_ILL_ARG;
    if(!cols && !vals)
    {
        *nnz=(s->fail==4 ? 1 : (s->fail==5 ? 0 : 5));
        return LIS_SUCCESS;
    }
    if(s->fail==5){*nnz=0;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<(s->fail==4 ? 1 : 5))
        return LIS_ERR_ILL_ARG;
    g=s->is+row;
    cols[0]=g;
    vals[0]=(LIS_SCALAR)4.0;
    if(s->fail==4){*nnz=1;return LIS_SUCCESS;}
    cols[1]=(g+1)%s->global_n;
    vals[1]=(LIS_SCALAR)0.5;
    if(s->fail==2 && s->rank==s->fail_rank && row==0)
        cols[1]=s->global_n;
    cols[2]=(g+s->global_n-1)%s->global_n;
    vals[2]=(LIS_SCALAR)0.25;
    /* Two distinct entries with the SAME off-rank global column.
     * Their contributions must both survive Ax and A^H x. */
    cols[3]=(g+s->global_n/2)%s->global_n;
    vals[3]=(LIS_SCALAR)0.2;
    cols[4]=cols[3];
    vals[4]=(LIS_SCALAR)0.125;
    *nnz=(s->fail==3 && s->rank==s->fail_rank && row==0) ? 4 : 5;
    return LIS_SUCCESS;
}
static LIS_SCALAR true_x(LIS_INT g)
{
    return (LIS_SCALAR)(1.0+0.1*(double)g);
}
static LIS_SCALAR ax_expected(LIS_INT g,LIS_INT n)
{
    return (LIS_SCALAR)4.0*true_x(g)
           +(LIS_SCALAR)0.5*true_x((g+1)%n)
           +(LIS_SCALAR)0.25*true_x((g+n-1)%n)
           +(LIS_SCALAR)0.325*true_x((g+n/2)%n);
}
static LIS_SCALAR ah_expected(LIS_INT g,LIS_INT n)
{
    return (LIS_SCALAR)4.0*true_x(g)
           +(LIS_SCALAR)0.25*true_x((g+1)%n)
           +(LIS_SCALAR)0.5*true_x((g+n-1)%n)
           +(LIS_SCALAR)0.325*true_x((g+n/2)%n);
}
int main(int argc,char **argv)
{
    ROWS rows;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,b=NULL,y=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,cap=0,err=LIS_SUCCESS,status=LIS_SUCCESS;
    int failed=0,global_failed=0,rank=0,nranks=1;
    if(lis_initialize(&argc,&argv))return 1;
#ifdef USE_MPI
    MPI_Comm_rank(LIS_COMM_WORLD,&rank);
    MPI_Comm_size(LIS_COMM_WORLD,&nranks);
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err){failed=1;goto done;}
    err=lis_matrix_set_size(A,0,12);if(err){failed=1;goto done;}
    rows.global_n=A->gn;rows.local_n=A->n;rows.is=A->is;
    rows.rank=rank;rows.fail_rank=(nranks>1?1:0);rows.fail=0;
    err=lis_matrix_set_user(A,&rows,NULL,NULL);if(err){failed=1;goto done;}
    err=lis_matrix_set_user_row(A,get_row);if(err){failed=1;goto done;}
    err=lis_matrix_get_capabilities(A,&cap);if(err){failed=1;goto done;}
    if(!(cap&LIS_MATRIX_CAP_MATVEC) ||
       !(cap&LIS_MATRIX_CAP_MATVECH)){failed=1;goto done;}
    err=lis_vector_duplicate(A,&x);if(err){failed=1;goto done;}
    err=lis_vector_duplicate(A,&y);if(err){failed=1;goto done;}
    err=lis_vector_duplicate(A,&b);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)x->value[i]=true_x(A->is+i);
    err=lis_matvec(A,x,y);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)(y->value[i]-ax_expected(A->is+i,A->gn)))>1e-11)
        {failed=1;goto done;}
    lis_vector_copy(y,b);
    err=lis_matvech(A,x,y);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)(y->value[i]-ah_expected(A->is+i,A->gn)))>1e-11)
        {failed=1;goto done;}

    err=lis_solver_create(&solver);if(err){failed=1;goto done;}
    err=lis_solver_set_option(
      "-i gmres -p none -scale none -restart 8 -maxiter 100 -tol 1e-10 -print none",
      solver);if(err){failed=1;goto done;}
    lis_vector_set_all((LIS_SCALAR)0.0,x);
    err=lis_solve(A,b,x,solver);
    if(err){failed=1;goto done;}
    err=lis_solver_get_status(solver,&status);
    if(err || status!=LIS_SUCCESS){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)(x->value[i]-true_x(A->is+i)))>1e-7)
        {failed=1;goto done;}

    /* Rows can change to local-only and empty between matvec calls:
     * the new exchange must not reuse a stale request pattern. */
    rows.fail=4;
    for(i=0;i<A->n;i++)x->value[i]=true_x(A->is+i);
    err=lis_matvec(A,x,y);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)(y->value[i]-4.0*true_x(A->is+i)))>1e-11)
        {failed=1;goto done;}
    err=lis_matvech(A,x,y);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)(y->value[i]-4.0*true_x(A->is+i)))>1e-11)
        {failed=1;goto done;}

    rows.fail=5;
    err=lis_matvec(A,x,y);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)y->value[i])>1e-11){failed=1;goto done;}
    err=lis_matvech(A,x,y);if(err){failed=1;goto done;}
    for(i=0;i<A->n;i++)
        if(fabs((double)y->value[i])>1e-11){failed=1;goto done;}

    rows.fail=0;
    if(nranks>1)
    {
        /* In-place Ax and A^H x require the local source to be saved
         * before overwriting the destination. */
        for(i=0;i<A->n;i++)y->value[i]=true_x(A->is+i);
        err=lis_matvec(A,y,y);if(err){failed=1;goto done;}
        for(i=0;i<A->n;i++)
            if(fabs((double)(y->value[i]-ax_expected(A->is+i,A->gn)))>1e-11)
            {failed=1;goto done;}
        for(i=0;i<A->n;i++)y->value[i]=true_x(A->is+i);
        err=lis_matvech(A,y,y);if(err){failed=1;goto done;}
        for(i=0;i<A->n;i++)
            if(fabs((double)(y->value[i]-ah_expected(A->is+i,A->gn)))>1e-11)
            {failed=1;goto done;}
    }
    /* An error on just one rank must be visible to *all* ranks. */
    rows.fail=1;
    err=lis_matvec(A,x,y);
    if(err!=LIS_ERR_ILL_ARG){failed=1;goto done;}
    err=lis_matvech(A,x,y);
    if(err!=LIS_ERR_ILL_ARG){failed=1;goto done;}
    /* Validate global column bounds and query/fill consistency on
     * one rank. The error must be broadcast collectively by both Ax
     * and A^H even though other ranks have valid rows. */
    rows.fail=2;
    err=lis_matvec(A,x,y);
    if(err!=LIS_ERR_ILL_ARG){failed=1;goto done;}
    err=lis_matvech(A,x,y);
    if(err!=LIS_ERR_ILL_ARG){failed=1;goto done;}
    rows.fail=3;
    err=lis_matvec(A,x,y);
    if(err!=LIS_ERR_ILL_ARG){failed=1;goto done;}
    err=lis_matvech(A,x,y);
    if(err!=LIS_ERR_ILL_ARG){failed=1;goto done;}
done:
    if(solver)lis_solver_destroy(solver);
    if(b)lis_vector_destroy(b);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
#ifdef USE_MPI
    MPI_Allreduce(&failed,&global_failed,1,MPI_INT,MPI_MAX,LIS_COMM_WORLD);
#else
    global_failed=failed;
#endif
    lis_finalize();
    if(rank==0)printf("USER_MPI_SPARSE %s ranks=%d\n",
        global_failed?"FAILED":"PASSED",nranks);
    return global_failed?1:0;
}
