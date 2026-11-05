/*
 * Deterministic ownership regression: inject a failure exactly at the
 * lis_solver_malloc_work[CG] dispatch, after preconditioner creation.
 * The direct kernel must leave borrowed precon allocated; lis_solve must
 * free its locally owned precon even when the kernel fails.
 */
#include <stdio.h>
#include "lis.h"
#include "lis_precon.h"
#include "lis_solver.h"

/* The fork defines this internal mutable dispatch table in lis_solver.c. */
extern LIS_SOLVER_MALLOC_WORK lis_solver_malloc_work[];

static LIS_PRECON seen_precon=NULL;
static LIS_INT observe_create(LIS_SOLVER solver,LIS_PRECON precon)
{
    (void)solver;
    seen_precon=precon;
    return LIS_SUCCESS;
}
static LIS_INT identity_solve(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    (void)solver;
    return lis_vector_copy(b,x);
}
static LIS_INT fail_work(LIS_SOLVER solver)
{
    (void)solver;
    return LIS_OUT_OF_MEMORY;
}
#define REQUIRE_OK(call) do { \
    LIS_INT result=(call); \
    if(result!=LIS_SUCCESS) { \
        fprintf(stderr,"%s failed, rc=%lld at line %d\n", \
                #call,(long long)result,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)
#define REQUIRE_TRUE(cond) do { \
    if(!(cond)) { \
        fprintf(stderr,"Assertion failed: %s at line %d\n", \
                #cond,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON owned=NULL;
    LIS_SCALAR value;
    LIS_SOLVER_MALLOC_WORK saved_work=NULL;
    LIS_INT err,i;
    int failed=0,hooked=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;
    REQUIRE_OK(lis_matrix_create(LIS_COMM_WORLD,&A));
    REQUIRE_OK(lis_matrix_set_size(A,0,3));
    for(i=0;i<3;i++) {
        REQUIRE_OK(lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));
    }
    REQUIRE_OK(lis_matrix_assemble(A));
    REQUIRE_OK(lis_vector_duplicate(A,&b));
    REQUIRE_OK(lis_vector_duplicate(A,&x));
    REQUIRE_OK(lis_vector_set_all(4.0,b));
    REQUIRE_OK(lis_solver_create(&solver));
    REQUIRE_OK(lis_precon_register("ownershiptest",
                                   observe_create,identity_solve,identity_solve));
    REQUIRE_OK(lis_solver_set_option(
        "-i cg -p ownershiptest -scale none -print none -tol 1.0e-12",
        solver));
    solver->A=A;
    REQUIRE_OK(lis_precon_create(solver,&owned));
    REQUIRE_TRUE(owned!=NULL && seen_precon==owned);
    REQUIRE_TRUE(lis_is_malloc(owned)==LIS_TRUE);

    saved_work=lis_solver_malloc_work[LIS_SOLVER_CG];
    lis_solver_malloc_work[LIS_SOLVER_CG]=fail_work;
    hooked=1;
    err=lis_solve_kernel(A,b,x,solver,owned);
    lis_solver_malloc_work[LIS_SOLVER_CG]=saved_work;
    hooked=0;
    REQUIRE_TRUE(err==LIS_OUT_OF_MEMORY);
    REQUIRE_TRUE(lis_is_malloc(owned)==LIS_TRUE);
    REQUIRE_TRUE(solver->retcode==LIS_OUT_OF_MEMORY);
    REQUIRE_OK(lis_precon_destroy(owned));
    owned=NULL;
    solver->precon=NULL;

    /* lis_solve owns the preconditioner it creates and must free it. */
    seen_precon=NULL;
    lis_solver_malloc_work[LIS_SOLVER_CG]=fail_work;
    hooked=1;
    err=lis_solve(A,b,x,solver);
    lis_solver_malloc_work[LIS_SOLVER_CG]=saved_work;
    hooked=0;
    REQUIRE_TRUE(err==LIS_OUT_OF_MEMORY);
    REQUIRE_TRUE(solver->retcode==LIS_OUT_OF_MEMORY);
    REQUIRE_TRUE(seen_precon!=NULL);
    REQUIRE_TRUE(lis_is_malloc(seen_precon)==LIS_FALSE);
    REQUIRE_TRUE(solver->precon==NULL);

    /* Restore the dispatch and verify normal execution still works. */
    REQUIRE_OK(lis_vector_set_all(0.0,x));
    REQUIRE_OK(lis_solve(A,b,x,solver));
    REQUIRE_TRUE(solver->retcode==LIS_SUCCESS);
    REQUIRE_TRUE(solver->precon==NULL);
    REQUIRE_OK(lis_vector_get_value(x,0,&value));
    REQUIRE_TRUE(value==(LIS_SCALAR)1.0);

    puts("PRECON_KERNEL_OWNERSHIP_PASSED");

cleanup:
    if(hooked) lis_solver_malloc_work[LIS_SOLVER_CG]=saved_work;
    if(owned && lis_is_malloc(owned)) lis_precon_destroy(owned);
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    lis_precon_register_free();
    lis_finalize();
    return failed;
}
