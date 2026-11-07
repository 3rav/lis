/* Standalone, serial deterministic MFAMG v3 V-cycle regression. */
#ifndef _WIN32
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#endif
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "lislib.h"

#define N 24
extern LIS_INT lis_precon_mfamg_get_vcycle_experiment_for_test(
    LIS_PRECON precon, LIS_INT *steps, LIS_INT *diag_calls);

static int set_steps(int steps)
{
#ifdef _WIN32
    if(steps==1) return _putenv("LIS_MFAMG_VCYCLE_SMOOTH_STEPS=");
    if(steps==2) return _putenv("LIS_MFAMG_VCYCLE_SMOOTH_STEPS=2");
    return _putenv("LIS_MFAMG_VCYCLE_SMOOTH_STEPS=3");
#else
    if(steps==1) return unsetenv("LIS_MFAMG_VCYCLE_SMOOTH_STEPS");
    return setenv("LIS_MFAMG_VCYCLE_SMOOTH_STEPS",
                  steps==2 ? "2" : "3",1);
#endif
}
static int set_diag(int enabled)
{
#ifdef _WIN32
    return _putenv(enabled ? "LIS_MFAMG_VCYCLE_DIAG=1"
                           : "LIS_MFAMG_VCYCLE_DIAG=");
#else
    return enabled ? setenv("LIS_MFAMG_VCYCLE_DIAG","1",1)
                   : unsetenv("LIS_MFAMG_VCYCLE_DIAG");
#endif
}

static int run_case(LIS_MATRIX A,LIS_VECTOR basis,
                    LIS_VECTOR b,LIS_VECTOR x,int expected)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR near[1],ax=NULL;
    LIS_INT err=0,got=0,called=0,i;
    LIS_REAL r2=0.0,bi,ri;
    int failed=1;
    near[0]=basis;
    if(set_steps(expected)!=0 || set_diag(1)!=0) return 1;
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i gmres -p mfamg -scale none -print none",solver);
    if(err) goto done;
    err=lis_solver_set_near_nullspace(solver,1,near);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_precon_create(solver,&precon);
    if(err) {precon=NULL;goto done;}
    solver->precon=precon;
    err=lis_precon_mfamg_get_vcycle_experiment_for_test(
        precon,&got,&called);
    if(err || got!=expected || called!=0) goto done;
    err=lis_psolve(solver,b,x);
    if(err) goto done;
    err=lis_precon_mfamg_get_vcycle_experiment_for_test(
        precon,&got,&called);
    if(err || called!=1) goto done;
    err=lis_psolve(solver,b,x);
    if(err) goto done;
    err=lis_precon_mfamg_get_vcycle_experiment_for_test(
        precon,&got,&called);
    if(err || called!=1) goto done;
    err=lis_vector_duplicate(A,&ax);
    if(err) goto done;
    err=lis_matvec(A,x,ax);
    if(err) goto done;
    for(i=0;i<N;i++) {
#ifdef _COMPLEX
        bi=(LIS_REAL)cabs(b->value[i]);
        ri=(LIS_REAL)cabs(b->value[i]-ax->value[i]);
#else
        bi=(LIS_REAL)fabs(b->value[i]);
        ri=(LIS_REAL)fabs(b->value[i]-ax->value[i]);
#endif
        r2+=ri*ri;
    }
    if(!isfinite((double)r2) || r2<0.0) goto done;
    printf("MFAMG VCYCLE steps=%d diag_calls=%d defect=%.9e PASS\n",
        expected,(int)called,(double)sqrt(r2));
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG VCYCLE test rc=%d steps=%d\n",
                    (int)err,expected);
    if(ax) lis_vector_destroy(ax);
    if(solver) solver->precon=NULL;
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR z=NULL,b=NULL,x=NULL;
    LIS_INT i,err=0;
    int failed=1,steps;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1) {lis_finalize();return 0;}
    }
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,N);
    if(err) goto done;
    for(i=0;i<N;i++) {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if(err) goto done;
        }
        if(i+1<N) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_vector_duplicate(A,&z);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    for(i=0;i<N;i++) {
        z->value[i]=1.0;
        b->value[i]=1.0+((LIS_REAL)(i%5))/5.0;
    }
    for(steps=1;steps<=3;steps++)
        if(run_case(A,z,b,x,steps)!=0) goto done;
    failed=0;
done:
    if(err) fprintf(stderr,"MFAMG VCYCLE setup rc=%d\n",(int)err);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(z) lis_vector_destroy(z);
    if(A) lis_matrix_destroy(A);
    set_steps(1);set_diag(0);
    lis_finalize();
    return failed;
}
