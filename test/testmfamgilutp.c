/*
 * Regression for opt-in serial MFAMG + pivoted ILUTP residual composition.
 * Covers actual lis_solve, repeated application, and missing near-nullspace.
 * This synthetic matrix is not an end-to-end contact simulation.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lis.h"
#include "lis_precon.h"

static LIS_INT
make_matrix(LIS_MATRIX *A,LIS_INT n)
{
    LIS_INT err,i,is,ie;
    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,n);
    if(err) return err;
    err=lis_matrix_get_range(*A,&is,&ie);
    if(err) return err;
    for(i=is;i<ie;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,*A);
        if(err) return err;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-0.5,*A);
            if(err) return err;
        }
        if(i<n-1)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,1.0,*A);
            if(err) return err;
        }
    }
    err=lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if(err) return err;
    return lis_matrix_assemble(*A);
}

static int
solve_profile(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR exact,LIS_VECTOR x,
              int composite)
{
    LIS_SOLVER solver=NULL;
    LIS_VECTOR modes[1];
    LIS_INT err,i;
    double max_error=0.0;
    int failed=0;
    modes[0]=exact;
    err=lis_solver_create(&solver);
    if(err) return 1;
    err=lis_solver_set_option(
        composite ?
        "-i fgmres -p mfamg_ilutp -scale none "
        "-ilutp_threshold 0.5 -ilut_equil 1 "
        "-ilut_drop 0.005 -ilut_rate 8 "
        "-mfamg_block_size 2 -restart 32 -maxiter 240 "
        "-tol 1e-9 -print none" :
        "-i fgmres -p mfamg -scale none "
        "-mfamg_block_size 2 -restart 32 -maxiter 240 "
        "-tol 1e-9 -print none",
        solver);
    if(err) { failed=1; goto cleanup; }
    err=lis_solver_set_near_nullspace(solver,1,modes);
    if(err) { failed=1; goto cleanup; }

    for(i=0;i<2;i++)
    {
        LIS_INT j;
        lis_vector_set_all(0.0,x);
        err=lis_solve(A,b,x,solver);
        if(err || solver->retcode!=LIS_SUCCESS)
        {
            fprintf(stderr,"%s solve #%d failed err=%d status=%d\n",
                composite?"MFAMG_ILUTP":"MFAMG",(int)i,
                (int)err,(int)solver->retcode);
            failed=1;
            goto cleanup;
        }
        max_error=0.0;
        for(j=0;j<x->n;j++)
        {
            double e=fabs((double)(x->value[j]-exact->value[j]));
            if(!isfinite(e)) { failed=1; goto cleanup; }
            if(e>max_error) max_error=e;
        }
        if(max_error>1e-6)
        {
            fprintf(stderr,"%s max exact solution error=%e\n",
                composite?"MFAMG_ILUTP":"MFAMG",max_error);
            failed=1;
            goto cleanup;
        }
    }
cleanup:
    if(solver!=NULL) lis_solver_destroy(solver);
    return failed;
}

static int
reject_missing_modes(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x)
{
    LIS_SOLVER solver=NULL;
    LIS_INT err;
    int failed=0;
    if(lis_solver_create(&solver)) return 1;
    err=lis_solver_set_option(
        "-i fgmres -p mfamg_ilutp -scale none "
        "-maxiter 20 -tol 1e-8 -print none",solver);
    if(err) { failed=1; goto cleanup; }
    lis_vector_set_all(0.0,x);
    err=lis_solve(A,b,x,solver);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"MFAMG_ILUTP unexpectedly accepted missing modes\n");
        failed=1;
    }
cleanup:
    lis_solver_destroy(solver);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,exact=NULL,x=NULL;
    LIS_INT err;
    int failed=0;
#ifdef _OPENMP
    int old_threads=omp_get_max_threads();
#endif
    err=lis_initialize(&argc,&argv);
    if(err) return 1;
#ifdef _OPENMP
    /* Global native ILUTP has a serial pivot-sequence restriction. */
    omp_set_num_threads(1);
#endif

    err=make_matrix(&A,128);
    if(err) { failed=1; goto cleanup; }
    if(A->nprocs!=1)
    {
        printf("MFAMG_ILUTP SKIP multi-rank test\n");
        goto cleanup;
    }
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&exact) ||
       lis_vector_duplicate(A,&x))
    {
        failed=1;
        goto cleanup;
    }
    lis_vector_set_all(1.0,exact);
    err=lis_matvec(A,exact,b);
    if(err) { failed=1; goto cleanup; }

    failed|=solve_profile(A,b,exact,x,0);
    if(!failed) failed|=solve_profile(A,b,exact,x,1);
    if(!failed) failed|=reject_missing_modes(A,b,x);
    if(!failed) printf("MFAMG_ILUTP serial regression PASSED\n");

cleanup:
    if(x) lis_vector_destroy(x);
    if(exact) lis_vector_destroy(exact);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
#ifdef _OPENMP
    omp_set_num_threads(old_threads);
#endif
    lis_finalize();
    return failed?1:0;
}
