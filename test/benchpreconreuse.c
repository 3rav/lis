/*
 * Controlled experiment: recreate a Jacobi preconditioner for every RHS
 * versus reuse one immutable preconditioner with lis_solve_kernel.
 *
 * Valid only for the same unchanged matrix, scale=none, default precision.
 * This is a benchmark, NOT a general-purpose preconditioner cache API.
 * Usage: benchpreconreuse [n=2000] [repetitions=12] [jacobi|ilu]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lis.h"
#include "lis_precon.h"

static int check_solution(LIS_VECTOR x,LIS_VECTOR expected)
{
    LIS_REAL error;
    LIS_INT err;
    err=lis_vector_axpy(-1.0,expected,x);
    if(err) return 1;
    err=lis_vector_nrm2(x,&error);
    if(err) return 1;
    if((double)error>1.0e-7) {
        fprintf(stderr,"solution error: %.12e\n",(double)error);
        return 1;
    }
    return 0;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR exact=NULL,b=NULL,x_fresh=NULL,x_reuse=NULL;
    LIS_SOLVER fresh=NULL,reuse=NULL;
    LIS_PRECON precon=NULL;
    LIS_INT n=2000,repetitions=12,i,round,err;
    LIS_INT iter_fresh=0,iter_reuse=0;
    const char *method=argc>3?argv[3]:"jacobi";
    const char *options=NULL;
    double fresh_s=0.0,reuse_s=0.0,create_s=0.0,t;
    int failed=0,precon_ownership_uncertain=0;

    if(argc>1) n=(LIS_INT)atol(argv[1]);
    if(argc>2) repetitions=(LIS_INT)atol(argv[2]);
    if(n<4 || n>20000 || repetitions<1 || repetitions>100 || argc>4) {
        fprintf(stderr,"Usage: %s [n 4..20000] [repetitions 1..100] [jacobi|ilu]\n",argv[0]);
        return 2;
    }
    if(strcmp(method,"jacobi")==0) {
        options="-i cg -p jacobi -scale none -tol 1.0e-11 -maxiter 200 -print none";
    } else if(strcmp(method,"ilu")==0) {
        options="-i bicgstab -p ilu -scale none -tol 1.0e-11 -maxiter 200 -print none";
    } else {
        fprintf(stderr,"Unsupported preconditioner: %s\n",method);
        return 2;
    }

    err=lis_initialize(&argc,&argv);
    if(err) return 1;
#define REQUIRE(call) do { \
    err=(call); \
    if(err) { \
        fprintf(stderr,"%s failed with code %lld at line %d\n", \
                #call,(long long)err,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

    REQUIRE(lis_matrix_create(LIS_COMM_WORLD,&A));
    REQUIRE(lis_matrix_set_size(A,0,n));
    REQUIRE(lis_matrix_malloc(A,3,NULL));
    for(i=0;i<n;i++) {
        REQUIRE(lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));
        if(i>0) {
            REQUIRE(lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A));
        }
        if(i+1<n) {
            REQUIRE(lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A));
        }
    }
    REQUIRE(lis_matrix_assemble(A));
    REQUIRE(lis_vector_duplicate(A,&exact));
    REQUIRE(lis_vector_duplicate(A,&b));
    REQUIRE(lis_vector_duplicate(A,&x_fresh));
    REQUIRE(lis_vector_duplicate(A,&x_reuse));
    REQUIRE(lis_solver_create(&fresh));
    REQUIRE(lis_solver_create(&reuse));

    /* This is the crucial restriction on the benchmark. */
    REQUIRE(lis_solver_set_option((char *)options,fresh));
    REQUIRE(lis_solver_set_option((char *)options,reuse));
    if(fresh->options[LIS_OPTIONS_SCALE]!=LIS_SCALE_NONE ||
       reuse->options[LIS_OPTIONS_SCALE]!=LIS_SCALE_NONE) {
        fprintf(stderr,"Refusing benchmark with in-place scaling enabled\n");
        failed=1;
        goto cleanup;
    }
    reuse->A=A;
    t=lis_wtime();
    REQUIRE(lis_precon_create(reuse,&precon));
    create_s=lis_wtime()-t;

    for(round=0;round<repetitions;round++) {
        for(i=0;i<n;i++) {
            LIS_SCALAR v=(LIS_SCALAR)(1.0+0.125*((i+round)%7));
            REQUIRE(lis_vector_set_value(LIS_INS_VALUE,i,v,exact));
        }
        REQUIRE(lis_matvec(A,exact,b));
        REQUIRE(lis_vector_set_all(0.0,x_fresh));
        REQUIRE(lis_vector_set_all(0.0,x_reuse));

        t=lis_wtime();
        REQUIRE(lis_solve(A,b,x_fresh,fresh));
        fresh_s+=lis_wtime()-t;
        if(fresh->retcode!=LIS_SUCCESS) {
            fprintf(stderr,"fresh solver did not converge at iteration %lld\n",
                    (long long)round);
            failed=1;
            goto cleanup;
        }
        iter_fresh=fresh->iter;
        if(check_solution(x_fresh,exact)) {
            fprintf(stderr,"fresh solution error, round=%lld\n",(long long)round);
            failed=1;
            goto cleanup;
        }

        t=lis_wtime();
        err=lis_solve_kernel(A,b,x_reuse,reuse,precon);
        reuse_s+=lis_wtime()-t;
        if(err) {
            /*
             * Some legacy kernel error paths destroy this borrowed precon.
             * Do not destroy it a second time on the benchmark error path.
             */
            precon_ownership_uncertain=1;
            fprintf(stderr,"kernel failed, rc=%lld\n",(long long)err);
            failed=1;
            goto cleanup;
        }
        if(reuse->retcode!=LIS_SUCCESS) {
            fprintf(stderr,"reuse solver did not converge, rc=%lld\n",
                    (long long)reuse->retcode);
            failed=1;
            goto cleanup;
        }
        iter_reuse=reuse->iter;
        if(check_solution(x_reuse,exact)) {
            fprintf(stderr,"reuse solution error, round=%lld\n",(long long)round);
            failed=1;
            goto cleanup;
        }
        if(iter_fresh!=iter_reuse) {
            fprintf(stderr,"iteration count mismatch: fresh=%lld reuse=%lld\n",
                    (long long)iter_fresh,(long long)iter_reuse);
            failed=1;
            goto cleanup;
        }
        printf("round=%lld iters=%lld verified=1\n",
               (long long)round,(long long)iter_fresh);
    }
    printf("precon=%s n=%lld repetitions=%lld fresh_total_s=%.9f "
           "precon_setup_s=%.9f reuse_kernel_s=%.9f "
           "reuse_total_s=%.9f\n",
           method,(long long)n,(long long)repetitions,
           fresh_s,create_s,reuse_s,create_s+reuse_s);
    puts("BENCH_PRECON_REUSE_PASSED");

cleanup:
    /* Solver borrows the reused precon, so break this association first. */
    if(reuse) reuse->precon=NULL;
    if(precon && !precon_ownership_uncertain) lis_precon_destroy(precon);
    if(reuse) lis_solver_destroy(reuse);
    if(fresh) lis_solver_destroy(fresh);
    if(x_reuse) lis_vector_destroy(x_reuse);
    if(x_fresh) lis_vector_destroy(x_fresh);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
