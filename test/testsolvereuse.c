/*
 * Public fixed-matrix reused-preconditioner API regression.
 *
 * Covers CG+Jacobi and BiCGSTAB+ILU(0), multiple right-hand sides,
 * explicit numeric refresh, changed options, and injected work OOM.
 */
#include <stdio.h>
#include <math.h>
#include "lis.h"
#include "lis_solver.h"

extern LIS_SOLVER_MALLOC_WORK lis_solver_malloc_work[];

static LIS_INT fail_work(LIS_SOLVER solver)
{
    (void)solver;
    return LIS_OUT_OF_MEMORY;
}
#define REQ(call) do { \
    LIS_INT result=(call); \
    if(result!=LIS_SUCCESS) { \
        fprintf(stderr,"%s failed rc=%lld line=%d\n", \
                #call,(long long)result,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)
#define VERIFY(cond) do { \
    if(!(cond)) { \
        fprintf(stderr,"Check failed %s at line %d\n",#cond,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

static int match_solution(LIS_VECTOR x,LIS_VECTOR exact)
{
    LIS_REAL norm;
    if(lis_vector_axpy(-1.0,exact,x)!=LIS_SUCCESS) return 0;
    if(lis_vector_nrm2(x,&norm)!=LIS_SUCCESS) return 0;
    return norm<1.0e-8;
}

static int run_method(LIS_MATRIX A,LIS_VECTOR exact,LIS_VECTOR b,
                      LIS_VECTOR x,LIS_VECTOR normal,int ilu)
{
    LIS_SOLVER solver=NULL,fresh=NULL;
    LIS_SOLVE_REUSE ctx=NULL;
    LIS_SOLVER_MALLOC_WORK saved=NULL;
    LIS_INT i,round,err,iter_before;
    LIS_INT mode=ilu?LIS_SOLVER_BICGSTAB:LIS_SOLVER_CG;
    const char *options=ilu
        ? "-i bicgstab -p ilu -scale none -tol 1.0e-11 -maxiter 200 -print none"
        : "-i cg -p jacobi -scale none -tol 1.0e-11 -maxiter 200 -print none";
    int failed=0,hooked=0;

    REQ(lis_solver_create(&solver));
    REQ(lis_solver_create(&fresh));
    REQ(lis_solver_set_option((char *)options,solver));
    REQ(lis_solver_set_option((char *)options,fresh));

    /* This API explicitly rejects unsafe scaling configurations. */
    REQ(lis_solver_set_option("-scale jacobi",solver));
    ctx=(LIS_SOLVE_REUSE)1; /* ensure the failed call clears output */
    err=lis_solve_reuse_create(A,solver,&ctx);
    VERIFY(err==LIS_ERR_NOT_IMPLEMENTED && ctx==NULL);
    REQ(lis_solver_set_option("-scale none",solver));

    REQ(lis_solve_reuse_create(A,solver,&ctx));
    VERIFY(ctx!=NULL);

    for(round=0;round<6;round++) {
        for(i=0;i<A->n;i++) {
            LIS_SCALAR v=(LIS_SCALAR)(1.0+0.125*((i+round)%7));
            REQ(lis_vector_set_value(LIS_INS_VALUE,i,v,exact));
        }
        REQ(lis_matvec(A,exact,b));
        REQ(lis_vector_set_all(0.0,x));
        REQ(lis_vector_set_all(0.0,normal));
        REQ(lis_solve_reuse_execute(ctx,b,x));
        VERIFY(solver->retcode==LIS_SUCCESS);
        iter_before=solver->iter;
        VERIFY(match_solution(x,exact));
        REQ(lis_solve(A,b,normal,fresh));
        VERIFY(fresh->retcode==LIS_SUCCESS && fresh->iter==iter_before);
        VERIFY(match_solution(normal,exact));
    }

    /* Reject changed options rather than applying stale cached factors. */
    REQ(lis_solver_set_option("-tol 1.0e-7",solver));
    err=lis_solve_reuse_execute(ctx,b,x);
    VERIFY(err==LIS_ERR_ILL_ARG);
    REQ(lis_solver_set_option("-tol 1.0e-11",solver));

    /* OOM in the work allocator must not destroy the cached preconditioner. */
    saved=lis_solver_malloc_work[mode];
    lis_solver_malloc_work[mode]=fail_work;
    hooked=1;
    err=lis_solve_reuse_execute(ctx,b,x);
    lis_solver_malloc_work[mode]=saved;
    hooked=0;
    VERIFY(err==LIS_OUT_OF_MEMORY);
    REQ(lis_vector_set_all(0.0,x));
    REQ(lis_solve_reuse_execute(ctx,b,x));
    VERIFY(match_solution(x,exact));

    /* In-place numeric update must be paired with an explicit refresh. */
    for(i=0;i<A->n;i++) {
        LIS_INT k;
        for(k=A->ptr[i];k<A->ptr[i+1];k++) {
            if(A->index[k]==i) A->value[k]=(LIS_SCALAR)5.0;
        }
    }
    REQ(lis_solve_reuse_refresh(ctx));
    REQ(lis_matvec(A,exact,b));
    REQ(lis_vector_set_all(0.0,x));
    REQ(lis_vector_set_all(0.0,normal));
    REQ(lis_solve_reuse_execute(ctx,b,x));
    VERIFY(match_solution(x,exact));
    REQ(lis_solve(A,b,normal,fresh));
    VERIFY(match_solution(normal,exact));

    /* Restore matrix for the next preconditioner case. */
    for(i=0;i<A->n;i++) {
        LIS_INT k;
        for(k=A->ptr[i];k<A->ptr[i+1];k++) {
            if(A->index[k]==i) A->value[k]=(LIS_SCALAR)4.0;
        }
    }

cleanup:
    if(hooked) lis_solver_malloc_work[mode]=saved;
    if(ctx && ctx!=(LIS_SOLVE_REUSE)1) lis_solve_reuse_destroy(ctx);
    if(fresh) lis_solver_destroy(fresh);
    if(solver) lis_solver_destroy(solver);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR exact=NULL,b=NULL,x=NULL,normal=NULL;
    LIS_INT err,i;
    int failed=0;
    const LIS_INT n=96;
    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    REQ(lis_matrix_create(LIS_COMM_WORLD,&A));
    REQ(lis_matrix_set_size(A,0,n));
    REQ(lis_matrix_malloc(A,3,NULL));
    for(i=0;i<n;i++) {
        REQ(lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));
        if(i>0) REQ(lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A));
        if(i+1<n) REQ(lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A));
    }
    REQ(lis_matrix_assemble(A));
    REQ(lis_vector_duplicate(A,&exact));
    REQ(lis_vector_duplicate(A,&b));
    REQ(lis_vector_duplicate(A,&x));
    REQ(lis_vector_duplicate(A,&normal));

    if(run_method(A,exact,b,x,normal,0) ||
       run_method(A,exact,b,x,normal,1)) {
        failed=1;
        goto cleanup;
    }

    puts("FIXED_MATRIX_PRECON_REUSE_PASSED");

cleanup:
    if(normal) lis_vector_destroy(normal);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
