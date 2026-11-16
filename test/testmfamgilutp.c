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
#include <string.h>
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
        "-mfamg_block_size 1 -restart 16 -maxiter 160 "
        "-tol 1e-9 -print none" :
        "-i fgmres -p mfamg -scale none "
        "-mfamg_block_size 1 -restart 16 -maxiter 160 "
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


/* Verify that all newly-owned state is initially NULL. */
static int test_zero_init(void)
{
    struct LIS_PRECON_STRUCT p;
    memset(&p,0xa5,sizeof(p));
    if(lis_precon_init(&p) || p.user_data!=NULL)
    {
        fprintf(stderr,"composite pointer not cleared by lis_precon_init\n");
        return 1;
    }
    return 0;
}

/* Check setup rejection while leaving the original matrix unmodified. */
static int reject_setup(LIS_MATRIX A,LIS_VECTOR exact,const char *label,
                        int adds,int bad_scale,int hybrid,int bad_threads,int extra)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR modes[1];
    LIS_INT err;
    int failed=0;
#ifdef _OPENMP
    int old_threads=omp_get_max_threads();
    if(bad_threads) omp_set_num_threads(2);
#else
    if(bad_threads) return 0;
#endif
    modes[0]=exact;
    err=lis_solver_create(&solver);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_option(hybrid ?
        "-i fgmres -p hybrid -hybrid_i gmres -hybrid_p mfamg_ilutp "
        "-scale none -print none" :
        "-i fgmres -p mfamg_ilutp -scale none -print none",solver);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_near_nullspace(solver,1,modes);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_matrix(A,solver);
    if(err) { failed=1; goto done; }
    if(adds) solver->options[LIS_OPTIONS_ADDS]=LIS_TRUE;
    if(bad_scale) solver->options[LIS_OPTIONS_SCALE]=LIS_SCALE_JACOBI;
    if(extra==1) solver->ilutp_threshold=0.0;
    if(extra==2) solver->ilut_ras_overlap=1;
    if(extra==3) solver->ilut_schur=1;
    err=lis_precon_create(solver,&precon);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"%s unexpectedly accepted\n",label);
        failed=1;
        lis_precon_destroy(precon);
    }
    /* A failed lis_precon_create has already destroyed *precon. */
    precon=NULL;
done:
    if(solver) lis_solver_destroy(solver);
#ifdef _OPENMP
    if(bad_threads) omp_set_num_threads(old_threads);
#endif
    return failed;
}

/* Six setup/apply/teardown cycles, each with an in-place (b==x) apply. */
static int test_owned_lifecycle(LIS_MATRIX A,LIS_VECTOR b,
                                LIS_VECTOR exact,LIS_VECTOR x)
{
    LIS_VECTOR reference=NULL,modes[1];
    int failed=0,k;
    LIS_INT err,j;
    modes[0]=exact;
    if(lis_vector_duplicate(A,&reference)) return 1;
    for(k=0;k<6;k++)
    {
        LIS_SOLVER solver=NULL;
        LIS_PRECON precon=NULL;
        err=lis_solver_create(&solver);
        if(err) { failed=1; goto finish_iteration; }
        err=lis_solver_set_option(
            "-i fgmres -p mfamg_ilutp -scale none "
            "-ilutp_threshold 0.5 -ilut_equil 1 -print none",solver);
        if(err) { failed=1; goto finish_iteration; }
        err=lis_solver_set_near_nullspace(solver,1,modes);
        if(err) { failed=1; goto finish_iteration; }
        err=lis_solver_set_matrix(A,solver);
        if(err) { failed=1; goto finish_iteration; }
        err=lis_precon_create(solver,&precon);
        if(err)
        {
            precon=NULL; /* factory already freed failed preconditioner */
            failed=1;
            goto finish_iteration;
        }
        solver->precon=precon;
        err=lis_psolve(solver,b,reference);
        if(err) { failed=1; goto finish_iteration; }
        err=lis_vector_copy(b,x);
        if(err) { failed=1; goto finish_iteration; }
        err=lis_psolve(solver,x,x);
        if(err) { failed=1; goto finish_iteration; }
        for(j=0;j<x->n;j++)
        {
            double d=fabs((double)(x->value[j]-reference->value[j]));
            if(!isfinite(d) || d>1e-10)
            {
                fprintf(stderr,"in-place apply differs from ordinary apply\n");
                failed=1;
                break;
            }
        }
finish_iteration:
        if(precon)
        {
            lis_precon_destroy(precon);
            if(solver) solver->precon=NULL;
        }
        if(solver) lis_solver_destroy(solver);
        if(failed) break;
    }
    lis_vector_destroy(reference);
    return failed;
}


/*
 * A contact/Newton application may mutate matrix values in place. Its
 * recorded numerical generation must invalidate both prebuilt factors.
 * Marking a mutation (without actually altering values) is sufficient.
 */
static int test_stale_matrix_rejection(LIS_MATRIX A,LIS_VECTOR b,
                                       LIS_VECTOR exact,LIS_VECTOR x)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR modes[1];
    LIS_INT err;
    int failed=0;
    modes[0]=exact;
    err=lis_solver_create(&solver);
    if(err) return 1;
    err=lis_solver_set_option(
        "-i fgmres -p mfamg_ilutp -scale none "
        "-ilutp_threshold 0.5 -print none",solver);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_near_nullspace(solver,1,modes);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_matrix(A,solver);
    if(err) { failed=1; goto done; }
    err=lis_precon_create(solver,&precon);
    if(err)
    {
        precon=NULL;
        failed=1;
        goto done;
    }
    solver->precon=precon;
    err=lis_psolve(solver,b,x);
    if(err) { failed=1; goto done; }

    err=lis_matrix_mark_values_changed(A);
    if(err) { failed=1; goto done; }
    err=lis_psolve(solver,b,x);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"stale matrix generation was not rejected\n");
        failed=1;
    }
done:
    if(precon)
    {
        lis_precon_destroy(precon);
        solver->precon=NULL;
    }
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

    err=make_matrix(&A,32);
    if(err) { failed=1; goto cleanup; }
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

    failed|=test_zero_init();
    if(!failed && A->nprocs!=1)
    {
        failed|=reject_setup(A,exact,"multi-rank matrix",0,0,0,0,0);
        if(!failed) printf("MFAMG_ILUTP multi-rank rejection PASSED\n");
        goto cleanup;
    }
#ifdef _COMPLEX
    if(!failed) failed|=reject_setup(A,exact,"complex scalar build",0,0,0,0,0);
    if(!failed) printf("MFAMG_ILUTP complex scalar rejection PASSED\n");
    goto cleanup;
#endif
    if(!failed) failed|=solve_profile(A,b,exact,x,1);
    if(!failed) failed|=solve_profile(A,b,exact,x,0);
    if(!failed) failed|=test_owned_lifecycle(A,b,exact,x);
    if(!failed) failed|=reject_missing_modes(A,b,x);
    if(!failed) failed|=reject_setup(A,exact,"Additive Schwarz",1,0,0,0,0);
    if(!failed) failed|=reject_setup(A,exact,"non-none scaling",0,1,0,0,0);
    if(!failed) failed|=reject_setup(A,exact,"Hybrid nesting",0,0,1,0,0);
    if(!failed) failed|=reject_setup(A,exact,"legacy ILUT threshold zero",0,0,0,0,1);
    if(!failed) failed|=reject_setup(A,exact,"ILUT RAS",0,0,0,0,2);
    if(!failed) failed|=reject_setup(A,exact,"ILUT Schur",0,0,0,0,3);
#ifdef _OPENMP
    if(!failed) failed|=reject_setup(A,exact,"multiple OpenMP threads",0,0,0,1,0);
#endif
    if(!failed) failed|=test_stale_matrix_rejection(A,b,exact,x);
    if(!failed) printf("MFAMG_ILUTP serial and rejection regressions PASSED\n");

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
