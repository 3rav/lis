/*
 * Repeated solves with the same LIS_SOLVER must release the preceding
 * residual-history allocation without changing the numerical result.
 */
#include <stdio.h>
#include "lis.h"

#define CHECK(expr) do { \
    LIS_INT rc=(expr); \
    if(rc!=LIS_SUCCESS) { \
        fprintf(stderr,"%s failed: %ld at line %d\n",#expr,(long)rc,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)
#define ASSERT_TRUE(expr) do { \
    if(!(expr)) { \
        fprintf(stderr,"Assertion failed: %s at line %d\n",#expr,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR exact=NULL,b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_REAL error_norm;
    LIS_REAL *previous_history=NULL;
    LIS_INT err,i,round;
    int failed=0;
    const LIS_INT n=12;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;
    CHECK(lis_matrix_create(LIS_COMM_WORLD,&A));
    CHECK(lis_matrix_set_size(A,0,n));

    /* Symmetric, strictly diagonally dominant tridiagonal matrix. */
    for(i=0;i<n;i++) {
        CHECK(lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));
        if(i>0)
            CHECK(lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A));
        if(i+1<n)
            CHECK(lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A));
    }
    CHECK(lis_matrix_assemble(A));
    CHECK(lis_vector_duplicate(A,&exact));
    CHECK(lis_vector_duplicate(A,&b));
    CHECK(lis_vector_duplicate(A,&x));
    CHECK(lis_solver_create(&solver));
    CHECK(lis_solver_set_option(
        "-i cg -p jacobi -scale none -tol 1.0e-12 -maxiter 200 -print none",
        solver));

    for(round=0;round<8;round++) {
        for(i=0;i<n;i++) {
            LIS_SCALAR value=(LIS_SCALAR)(1.0+(double)(i+round)*0.125);
            CHECK(lis_vector_set_value(LIS_INS_VALUE,i,value,exact));
        }
        CHECK(lis_matvec(A,exact,b));
        CHECK(lis_vector_set_all(0.0,x));

        previous_history=solver->rhistory;
        CHECK(lis_solve(A,b,x,solver));
        ASSERT_TRUE(solver->retcode==LIS_SUCCESS);
        ASSERT_TRUE(solver->rhistory!=NULL);
        ASSERT_TRUE(lis_is_malloc(solver->rhistory)==LIS_TRUE);
        if(previous_history!=NULL) {
            ASSERT_TRUE(previous_history!=solver->rhistory);
            ASSERT_TRUE(lis_is_malloc(previous_history)==LIS_FALSE);
        }
        CHECK(lis_vector_axpy(-1.0,exact,x));
        CHECK(lis_vector_nrm2(x,&error_norm));
        ASSERT_TRUE(error_norm<1.0e-8);
    }
    puts("SOLVER_RESIDUAL_HISTORY_REUSE_PASSED");

cleanup:
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
