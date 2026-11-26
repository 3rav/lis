/* Internal FGMRES stop policy regression; no reference driver dependency. */
#include "lislib.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(EXPECT_INT64)
typedef char check_lis_int64[(sizeof(LIS_INT)>=8)?1:-1];
#endif
#if defined(EXPECT_MPI) && !defined(USE_MPI)
#error MPI configuration not enabled
#endif
#if defined(EXPECT_COMPLEX) && !defined(USE_COMPLEX) && !defined(_COMPLEX)
#error complex configuration not enabled
#endif
#if defined(EXPECT_QUAD) && !defined(USE_QUAD_PRECISION)
#error quad configuration not enabled
#endif

typedef struct { LIS_INT at,last,calls; int invalid; } history;
static LIS_INT stop(LIS_INT iter,LIS_REAL res,void *data)
{
    history *h=data;
    if(iter!=h->last+1 || !isfinite((double)res)){
        h->invalid=1;
        return LIS_TRUE;
    }
    h->last=iter;
    ++h->calls;
    return iter==h->at;
}

#if defined(USE_MPI) || defined(_COMPLEX) || defined(USE_COMPLEX) || defined(_LONG__DOUBLE) || defined(USE_QUAD_PRECISION)
int main(void)
{
    if(lis_solve_kernel_fgmres_stop(NULL,NULL,NULL,NULL,NULL,stop,NULL)!=LIS_ERR_NOT_IMPLEMENTED)
        return 1;
    puts("FGMRES_STOP_UNSUPPORTED=PASS");
    return 0;
}
#else
static LIS_INT matrix(int identity,LIS_MATRIX *A)
{
    LIS_INT i,rc;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,A)) ||
       (rc=lis_matrix_set_size(*A,0,80))) return rc;
    for(i=0;i<80;i++)
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,identity?1.0:(double)(i+1),*A)))return rc;
    if((rc=lis_matrix_set_type(*A,LIS_MATRIX_CSR)))return rc;
    return lis_matrix_assemble(*A);
}
static LIS_INT prepare(LIS_MATRIX P,LIS_SOLVER *s,LIS_PRECON *pc,LIS_INT maxiter)
{
    LIS_INT rc;
    char options[]="-i fgmres -p ilut -scale none -print none -restart 8 -maxiter 100 -tol 1e-14 -ilutp 1 -ilutp_full 1 -ilutp_full_limit 60000000 -initx_zeros true";
    if((rc=lis_solver_create(s)) ||
       (rc=lis_solver_set_option(options,*s)) ||
       (rc=lis_solver_set_num_threads(*s,1)) ||
       (rc=lis_solver_set_maxiter(*s,maxiter)))return rc;
    return lis_solver_sequence_prepare_precon(P,*s,pc);
}
static int run(LIS_MATRIX P,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_INT maxiter,history *h,LIS_INT *status)
{
    LIS_SOLVER s=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT rc=prepare(P,&s,&pc,maxiter);
    if(!rc)rc=h?lis_solve_kernel_fgmres_stop(A,b,x,s,pc,stop,h):lis_solve_kernel(A,b,x,s,pc);
    if(!rc)rc=lis_solver_get_status(s,status);
    if(s)lis_solver_destroy(s);
    return rc!=LIS_SUCCESS;
}
static int compare(LIS_MATRIX P,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR y,LIS_INT limit,LIS_INT at,LIS_INT expected,LIS_INT calls)
{
    LIS_INT a,c;
    history h={at,0,0,0};
    if(run(P,A,b,x,limit,NULL,&a)||run(P,A,b,y,100,&h,&c))return 1;
    printf("FGMRES_STOP_CASE at=%lld limit=%lld calls=%lld status=%lld/%lld\n",(long long)at,(long long)limit,(long long)h.calls,(long long)a,(long long)c);
    return a!=c || (expected>=0 && a!=expected) || h.invalid ||
        (calls>=0 && h.calls!=calls) ||
        memcmp(x->value,y->value,(size_t)x->n*sizeof(LIS_SCALAR))!=0;
}
static int thread_guard(LIS_MATRIX P,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x)
{
#ifdef _OPENMP
    LIS_SOLVER s=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT rc,parallel_rc=LIS_SUCCESS;
    history h={1,0,0,0};
    rc=prepare(P,&s,&pc,100);
    if(rc){if(s)lis_solver_destroy(s);return 1;}
    omp_set_num_threads(2);
    rc=lis_solve_kernel_fgmres_stop(A,b,x,s,pc,stop,&h);
    omp_set_num_threads(1);
    #pragma omp parallel num_threads(2)
    {
        #pragma omp single
        {
            /* Force a one-thread nested team setting, but remain in parallel. */
            omp_set_num_threads(1);
            parallel_rc=lis_solve_kernel_fgmres_stop(A,b,x,s,pc,stop,&h);
        }
    }
    lis_solver_destroy(s);
    if(rc!=LIS_ERR_ILL_ARG||parallel_rc!=LIS_ERR_ILL_ARG||h.calls)return 1;
#else
    (void)P;(void)A;(void)b;(void)x;
#endif
    return 0;
}
int main(int argc,char **argv)
{
    LIS_MATRIX P=NULL,A=NULL;
    LIS_VECTOR b=NULL,x=NULL,y=NULL;
    int fail=1;
    LIS_INT at;
    if(lis_initialize(&argc,&argv))return 1;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif
    if(matrix(1,&P)||matrix(0,&A)||lis_vector_duplicate(P,&b)||
       lis_vector_duplicate(P,&x)||lis_vector_duplicate(P,&y)||lis_vector_set_all(1,b))goto done;
    /* m-1, m, m+1: compare the actual approximation with an ordinary
     * FGMRES solve capped at the same iteration, including a restart. */
    for(at=7;at<=9;at++)
        if(compare(P,A,b,x,y,at,at,LIS_MAXITER,at))goto done;
    {
        LIS_INT ordinary,hooked;
        history h={9,0,0,0};
        if(run(P,A,b,x,9,NULL,&ordinary)||run(P,A,b,y,9,&h,&hooked)||
           ordinary!=LIS_MAXITER||hooked!=ordinary||h.calls!=9||h.invalid||
           memcmp(x->value,y->value,(size_t)x->n*sizeof(LIS_SCALAR)))goto done;
    }
    if(compare(P,A,b,x,y,100,0,-1,-1))goto done;
    /* Convergence and requested stop in the same iteration. */
    if(compare(P,P,b,x,y,100,1,LIS_SUCCESS,1))goto done;
    if(thread_guard(P,A,b,x))goto done;
    if(lis_vector_set_all(0,b)||compare(P,A,b,x,y,100,0,LIS_SUCCESS,0))goto done;
    fail=0;
done:
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    if(P)lis_matrix_destroy(P);
    lis_finalize();
    printf("FGMRES_STOP_REGRESSION=%s\n",fail?"FAIL":"PASS");
    return fail;
}
#endif
