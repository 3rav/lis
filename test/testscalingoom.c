/*
 * Fault injection for an early lis_solve_kernel scaling allocation failure.
 *
 * Run with GNU ld --wrap=lis_malloc (test/Makefile.am).
 * The first vector duplicate in the kernel creates xx; after history
 * allocation, the second vector duplicate (solver->d) is forced to fail.
 *
 * Regression criterion: xx and rhistory must no longer be live, while the
 * caller-owned preconditioner is still live and a subsequent solve succeeds.
 */
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include "lis.h"
#include "lis_precon.h"

extern void *__real_lis_malloc(size_t size, char *tag);

static int fault_active=0;
static int saw_history=0;
static int injected=0;
static int first_vectors=0;
static void *observed_xx=NULL;
static void *observed_history=NULL;

void *__wrap_lis_malloc(size_t size, char *tag)
{
    void *p;
    if(fault_active && tag && strcmp(tag,"lis_vector_duplicateex::vout")==0) {
        if(saw_history) {
            ++injected;
            return NULL;
        }
        ++first_vectors;
        p=__real_lis_malloc(size,tag);
        if(observed_xx==NULL) observed_xx=p;
        return p;
    }
    p=__real_lis_malloc(size,tag);
    if(fault_active && tag && strcmp(tag,"lis_solve::rhistory")==0) {
        observed_history=p;
        saw_history=1;
    }
    return p;
}

#define REQUIRE(call) do { \
    LIS_INT rc=(call); \
    if(rc!=LIS_SUCCESS) { \
        fprintf(stderr,"%s returned %lld at line %d\n", \
                #call,(long long)rc,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

#define ASSERT_TRUE(cond) do { \
    if(!(cond)) { \
        fprintf(stderr,"Assertion failed: %s at line %d\n",#cond,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON borrowed=NULL;
    LIS_INT err,i;
    LIS_SCALAR result;
    int failed=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    REQUIRE(lis_matrix_create(LIS_COMM_WORLD,&A));
    REQUIRE(lis_matrix_set_size(A,0,4));
    for(i=0;i<4;i++) {
        REQUIRE(lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));
    }
    REQUIRE(lis_matrix_assemble(A));
    REQUIRE(lis_vector_duplicate(A,&b));
    REQUIRE(lis_vector_duplicate(A,&x));
    REQUIRE(lis_vector_set_all(4.0,b));
    REQUIRE(lis_vector_set_all(0.0,x));
    REQUIRE(lis_solver_create(&solver));

    /* CG+Jacobi scaling will ask for solver->d after xx and rhistory. */
    REQUIRE(lis_solver_set_option(
        "-i cg -p none -scale jacobi -tol 1.0e-12 -print none",
        solver));
    solver->A=A;
    REQUIRE(lis_precon_create(solver,&borrowed));
    ASSERT_TRUE(lis_is_malloc(borrowed)==LIS_TRUE);
    ASSERT_TRUE(solver->rhistory==NULL);

    fault_active=1;
    err=lis_solve_kernel(A,b,x,solver,borrowed);
    fault_active=0;

    ASSERT_TRUE(err==LIS_OUT_OF_MEMORY);
    ASSERT_TRUE(solver->retcode==LIS_OUT_OF_MEMORY);
    ASSERT_TRUE(first_vectors==1 && saw_history && injected==1);
    ASSERT_TRUE(observed_xx!=NULL && observed_history!=NULL);
    ASSERT_TRUE(lis_is_malloc(observed_xx)==LIS_FALSE);
    ASSERT_TRUE(lis_is_malloc(observed_history)==LIS_FALSE);
    ASSERT_TRUE(lis_is_malloc(borrowed)==LIS_TRUE);
    ASSERT_TRUE(solver->rhistory==NULL);
    ASSERT_TRUE(solver->work==NULL);
    REQUIRE(lis_precon_destroy(borrowed));
    borrowed=NULL;

    /* Recover on the same solver object after the injected failure. */
    REQUIRE(lis_solver_set_option(
        "-i cg -p none -scale none -tol 1.0e-12 -print none",
        solver));
    REQUIRE(lis_solve(A,b,x,solver));
    ASSERT_TRUE(solver->retcode==LIS_SUCCESS);
    REQUIRE(lis_vector_get_value(x,0,&result));
    ASSERT_TRUE(result==(LIS_SCALAR)1.0);
    puts("SOLVER_SCALING_ALLOCATION_FAILURE_PASSED");

cleanup:
    fault_active=0;
    if(borrowed && lis_is_malloc(borrowed)) lis_precon_destroy(borrowed);
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
