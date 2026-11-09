/*
 * Standalone (non-MFAMG) regression of within-solve FGCRO-DR recycling.
 * Uses the same matrix and right-hand side with the feature OFF and ON.
 * Verifies a real restart, recycle projection, true residual, and AZ=C.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"
#include "lis_solver.h"
#include "lislib.h"

#ifdef _COMPLEX
int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    puts("SKIP: FGCRO-DR positive-k real-only regression");
    return 77;
}
#else
#define N 120

static int set_intra(int enabled)
{
#ifdef _WIN32
    if(_putenv(enabled
        ? "LIS_FGCRODR_INTRA_SOLVE_RECYCLE=1"
        : "LIS_FGCRODR_INTRA_SOLVE_RECYCLE=0")!=0)
        return 1;
    if(_putenv("LIS_FGMRES_PROJECTED_RELATIVE=1")!=0)
        return 1;
#else
    if(setenv("LIS_FGCRODR_INTRA_SOLVE_RECYCLE",
              enabled ? "1" : "0",1)!=0)
        return 1;
    if(setenv("LIS_FGMRES_PROJECTED_RELATIVE","1",1)!=0)
        return 1;
#endif
    return 0;
}

static int solve_one(LIS_MATRIX A, LIS_VECTOR b, int intra)
{
    LIS_VECTOR x=NULL, ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT rc,dim=0,action=0,iter=0,status=0;
    LIS_UNSIGNED_INT builds=0,updates=0,uses=0,rejects=0;
    LIS_REAL orth=0.0,relation=0.0,bnorm=0.0,rnorm=0.0;
    int failed=1;
    if(set_intra(intra)) goto done;
    if((rc=lis_vector_duplicate(A,&x))) goto done;
    if((rc=lis_vector_duplicate(A,&ax))) goto done;
    if((rc=lis_solver_create(&solver))) goto done;
    if((rc=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 2 -maxiter 300 "
        "-tol 1e-5 -fgcrodr_k 1 -print none",solver))) goto done;
    if((rc=lis_vector_set_all(0.0,x))) goto done;
    rc=lis_solve(A,b,x,solver);
    status=solver->retcode;
    iter=solver->iter;
    if(rc!=LIS_SUCCESS || status!=LIS_SUCCESS || iter<=2)
        goto done;
    if((rc=lis_solver_get_fgcrodr_recycle(
        solver,&dim,&action,NULL,NULL))) goto done;
    if((rc=lis_fgcrodr_get_stats(
        solver,&builds,&updates,&uses,&rejects))) goto done;
    if((rc=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((rc=lis_matvec(A,x,ax))) goto done;
    if((rc=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((rc=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((rc=lis_vector_nrm2(b,&bnorm))) goto done;
    if(bnorm==0.0 || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8 ||
       builds<1 || rejects!=0 || dim!=1 ||
       (intra && (uses!=1 || updates<1)) ||
       (!intra && uses!=0))
        goto done;
    printf("FGCRODR_INTRA_%s_OK iter=%d true_relative=%.12e "
           "k=%d builds=%llu updates=%llu uses=%llu rejects=%llu "
           "orth=%.6e relation=%.6e\n",
           intra?"ON":"OFF",(int)iter,(double)(rnorm/bnorm),
           (int)dim,(unsigned long long)builds,
           (unsigned long long)updates,(unsigned long long)uses,
           (unsigned long long)rejects,(double)orth,(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_%s_FAIL status=%d iter=%d "
                "dim=%d updates=%llu uses=%llu rejects=%llu "
                "orth=%.6e relation=%.6e rnorm=%.6e\n",
                intra?"ON":"OFF",(int)status,(int)iter,(int)dim,
                (unsigned long long)updates,(unsigned long long)uses,
                (unsigned long long)rejects,(double)orth,
                (double)relation,(double)rnorm);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    return failed;
}


/* A synthetic changed-operator sequence.  Reusing the solver across
 * distinct matrix objects must not project the new residual against
 * a stale A1*Z=C relation for A2.
 */
static int check_changed_operator(LIS_MATRIX A1, LIS_VECTOR b)
{
    LIS_MATRIX A2=NULL,stored=NULL;
    LIS_VECTOR x=NULL,ax=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,valid=0,failed=1;
    LIS_UNSIGNED_INT sg=0,vg=0;
    LIS_REAL orth=0.0,relation=0.0,rnorm=0.0,bnorm=0.0;

    if(set_intra(1)) goto done;
    if((err=lis_matrix_create(LIS_COMM_WORLD,&A2))) goto done;
    if((err=lis_matrix_set_size(A2,0,N))) goto done;
    for(i=0;i<N;i++)
    {
        if((err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.25,A2)))
            goto done;
        if(i>0 &&
           (err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A2)))
            goto done;
        if(i<N-1 &&
           (err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.5,A2)))
            goto done;
    }
    if((err=lis_matrix_set_type(A2,LIS_MATRIX_CSR))) goto done;
    if((err=lis_matrix_assemble(A2))) goto done;
    if((err=lis_vector_duplicate(A1,&x))) goto done;
    if((err=lis_vector_duplicate(A1,&ax))) goto done;
    if((err=lis_solver_create(&solver))) goto done;
    if((err=lis_solver_set_option(
        "-i fgcrodr -p jacobi -restart 2 -maxiter 300 "
        "-tol 1e-5 -fgcrodr_k 1 -print none",solver))) goto done;

    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A1,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&sg,&vg))) goto done;
    if(!valid || stored!=A1) goto done;

    /* New assembled matrix object, same layout, modified diagonal. */
    if((err=lis_vector_set_all(0.0,x))) goto done;
    err=lis_solve_sequence(A2,b,x,solver);
    if(err!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS) goto done;
    if((err=lis_fgcrodr_get_lifecycle_state(
        solver,&valid,&stored,&sg,&vg))) goto done;
    if(!valid || stored!=A2 ||
       sg!=A2->structure_generation ||
       vg!=A2->value_generation)
        goto done;
    if((err=lis_fgcrodr_validate_recycle(
        solver,&orth,&relation))) goto done;
    if((err=lis_matvec(A2,x,ax))) goto done;
    if((err=lis_vector_xpay(b,-1.0,ax))) goto done;
    if((err=lis_vector_nrm2(ax,&rnorm))) goto done;
    if((err=lis_vector_nrm2(b,&bnorm))) goto done;
    if(bnorm==0.0 || !isfinite(rnorm) || rnorm/bnorm>1e-5 ||
       !isfinite(orth) || !isfinite(relation) ||
       orth>1e-8 || relation>1e-8)
        goto done;
    printf("FGCRODR_INTRA_CHANGED_OPERATOR_OK "
           "A1_to_A2_status=%d residual=%.12e relation=%.6e\n",
           (int)solver->retcode,(double)(rnorm/bnorm),(double)relation);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_CHANGED_OPERATOR_FAIL "
                "api=%d status=%d valid=%d stored_matches_A2=%d "
                "relation=%.6e\n",
                (int)err,solver?(int)solver->retcode:-1,
                (int)valid,(int)(stored==A2),(double)relation);
    if(solver) lis_solver_destroy(solver);
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(A2) lis_matrix_destroy(A2);
    return failed;
}

int main(int argc, char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR u=NULL,b=NULL;
    LIS_INT rc=0,i,failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,&A))) goto done;
    if((rc=lis_matrix_set_size(A,0,N))) goto done;
    for(i=0;i<N;i++)
    {
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A)))
            goto done;
        if(i>0 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A)))
            goto done;
        if(i<N-1 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.5,A)))
            goto done;
    }
    if((rc=lis_matrix_set_type(A,LIS_MATRIX_CSR))) goto done;
    if((rc=lis_matrix_assemble(A))) goto done;
    if((rc=lis_vector_duplicate(A,&u))) goto done;
    if((rc=lis_vector_duplicate(A,&b))) goto done;
    for(i=0;i<N;i++)
        u->value[i]=(LIS_SCALAR)(((i*37)%23)-11.0);
    if((rc=lis_matvec(A,u,b))) goto done;
    if(solve_one(A,b,0) || solve_one(A,b,1)) goto done;
    if(check_changed_operator(A,b)) goto done;
    failed=0;
done:
    if(failed)
        fprintf(stderr,"FGCRODR_INTRA_REGRESSION_FAIL rc=%d\n",(int)rc);
    if(b) lis_vector_destroy(b);
    if(u) lis_vector_destroy(u);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
#endif
