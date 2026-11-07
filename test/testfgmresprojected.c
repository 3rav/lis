/* FGMRES v5: projected residual must use the solver's relative norm. */
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
#include <stdio.h>
#include <stdlib.h>
#include "lislib.h"

static int set_fix(int enable)
{
#ifdef _WIN32
    return _putenv(enable ? "LIS_FGMRES_PROJECTED_RELATIVE=1" :
                            "LIS_FGMRES_PROJECTED_RELATIVE=");
#else
    return enable ? setenv("LIS_FGMRES_PROJECTED_RELATIVE","1",1)
                  : unsetenv("LIS_FGMRES_PROJECTED_RELATIVE");
#endif
}

static int run_scaled_case(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,
                           LIS_VECTOR r, LIS_REAL magnitude)
{
    LIS_SOLVER solver=NULL;
    LIS_INT err;
    LIS_REAL rn=0.0,bn=0.0;
    int failed=1;
    b->value[0]=magnitude;
    b->value[1]=magnitude;
    err=lis_vector_set_all(0.0,x);
    if(err) return 1;
    err=lis_solver_create(&solver);
    if(err) return 1;
    err=lis_solver_set_option(
        "-i fgmres -p none -scale none -print none "
        "-tol 0.3 -maxiter 1 -restart 2",solver);
    if(err) goto done;
    err=lis_solve(A,b,x,solver);
    if(err) goto done;
    if(solver->retcode!=LIS_SUCCESS || solver->iter!=1) goto done;
    err=lis_matvec(A,x,r);
    if(err) goto done;
    err=lis_vector_xpay(b,-1.0,r);
    if(err) goto done;
    err=lis_vector_nrm2(r,&rn);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || bn==0.0) goto done;
    if(!(rn/bn<0.3)) goto done;
    printf("FGMRES projected relative scale=%.3e iter=%d true=%.9e PASS\\n",
           (double)magnitude,(int)solver->iter,(double)(rn/bn));
    failed=0;
done:
    if(failed)
        fprintf(stderr,
            "FGMRES projected relative FAILED: scale=%.3e rc=%d status=%d iter=%d\\n",
            (double)magnitude,(int)err,
            solver ? (int)solver->retcode : -1,
            solver ? (int)solver->iter : -1);
    if(solver) lis_solver_destroy(solver);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,r=NULL;
    LIS_INT err=0;
    int failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1) {lis_finalize();return 0;}
    }
#endif
    if(set_fix(1)!=0) goto done;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,2);
    if(err) goto done;
    err=lis_matrix_set_value(LIS_INS_VALUE,0,0,2.0,A);
    if(err) goto done;
    err=lis_matrix_set_value(LIS_INS_VALUE,1,1,3.0,A);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&r);
    if(err) goto done;
    if(run_scaled_case(A,b,x,r,10000.0)!=0) goto done;
    if(run_scaled_case(A,b,x,r,0.01)!=0) goto done;
    failed=0;
done:
    if(err) fprintf(stderr,"FGMRES projected setup rc=%d\\n",(int)err);
    if(r) lis_vector_destroy(r);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    set_fix(0);
    lis_finalize();
    return failed;
}
