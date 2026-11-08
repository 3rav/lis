/* Regression: opt-in larger fine aggregates for multivector MFAMG.
 * Default target=2 remains unchanged.  With the opt-in enabled two
 * independent algebraic smooth vectors must allow useful coarsening.
 */
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
#include <stdio.h>
#include <stdlib.h>
#include "lislib.h"

#define TEST_N 240
#define TEST_LEVEL_CAP 32
#define TEST_VARIABLE "LIS_MFAMG_EXPERIMENT_FINE_AGG"

extern LIS_INT lis_precon_mfamg_get_hierarchy_info(
    LIS_PRECON precon, LIS_INT max_levels, LIS_INT *nlevels,
    LIS_INT *coarse_n, LIS_INT *level_n, LIS_INT *owns_A);

static int toggle_experiment(int enabled)
{
#ifdef _WIN32
    return _putenv(enabled ? TEST_VARIABLE "=1" : TEST_VARIABLE "=");
#else
    return enabled ? setenv(TEST_VARIABLE,"1",1) : unsetenv(TEST_VARIABLE);
#endif
}

static int set_factor(int factor)
{
#ifdef _WIN32
    if(factor==2) return _putenv("LIS_MFAMG_FINE_AGG_FACTOR=2");
    if(factor==3) return _putenv("LIS_MFAMG_FINE_AGG_FACTOR=3");
    if(factor==4) return _putenv("LIS_MFAMG_FINE_AGG_FACTOR=4");
    return _putenv("LIS_MFAMG_FINE_AGG_FACTOR=");
#else
    if(factor==0) return unsetenv("LIS_MFAMG_FINE_AGG_FACTOR");
    return setenv("LIS_MFAMG_FINE_AGG_FACTOR",
        factor==2 ? "2" : (factor==3 ? "3" : "4"),1);
#endif
}

static int run_test(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_VECTOR basis[2]={NULL,NULL};
    LIS_INT levels[TEST_LEVEL_CAP], owned[TEST_LEVEL_CAP];
    LIS_INT nlevels=0, coarse_n=0, err=0;
    LIS_INT i, factor, previous_coarse=TEST_N;
    int failed=1;

    if(toggle_experiment(1)!=0) {
        fprintf(stderr,"MFAMG FINE AGG: cannot set environment\n");
        return 1;
    }
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,TEST_N);
    if(err) goto done;
    for(i=0;i<TEST_N;i++) {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A);
        if(err) goto done;
        if(i>0) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if(err) goto done;
        }
        if(i+1<TEST_N) {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_vector_duplicate(A,&basis[0]);
    if(err) goto done;
    err=lis_vector_duplicate(A,&basis[1]);
    if(err) goto done;
    for(i=0;i<TEST_N;i++) {
        basis[0]->value[i] = 1.0;
        basis[1]->value[i] = ((double)i-(TEST_N-1)/2.0)/(double)TEST_N;
    }
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option("-i gmres -p mfamg -scale none -print none",solver);
    if(err) goto done;
    err=lis_solver_set_near_nullspace(solver,2,basis);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;

    /*
     * On the connected 1-D matrix with two independent smooth modes,
     * increasing the first-level target from 4 to 6 to 8 must reduce
     * the first coarse space without losing valid hierarchy levels.
     * This exercises the real preconditioner creation path, not merely
     * the environment parser.
     */
    for(factor=2;factor<=4;factor++)
    {
        if(set_factor((int)factor)!=0) goto done;
        err=lis_precon_create(solver,&precon);
        if(err) goto done;
        err=lis_precon_mfamg_get_hierarchy_info(
            precon,TEST_LEVEL_CAP,&nlevels,&coarse_n,levels,owned);
        if(err) goto done;
        if(nlevels<2 || levels[0]!=TEST_N ||
           levels[1]>=TEST_N || levels[1]>TEST_N/2 ||
           levels[1]>=previous_coarse || coarse_n<=0) {
            fprintf(stderr,
                "MFAMG FINE AGG factor=%d: bad hierarchy nlevels=%d "
                "first=%d previous=%d terminal=%d\n",
                (int)factor,(int)nlevels,(int)levels[1],
                (int)previous_coarse,(int)coarse_n);
            goto done;
        }
        printf("MFAMG FINE AGG: factor=%d n=%d coarse_level1=%d "
               "nlevels=%d terminal=%d PASS\n",
               (int)factor,TEST_N,(int)levels[1],
               (int)nlevels,(int)coarse_n);
        previous_coarse=levels[1];
        lis_precon_destroy(precon);
        precon=NULL;
    }
    failed=0;

done:
    if(err) fprintf(stderr,"MFAMG FINE AGG: LIS error=%d\n",(int)err);
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    if(basis[1]) lis_vector_destroy(basis[1]);
    if(basis[0]) lis_vector_destroy(basis[0]);
    if(A) lis_matrix_destroy(A);
    if(toggle_experiment(0)!=0 || set_factor(0)!=0) failed=1;
    return failed;
}

int main(int argc,char **argv)
{
    int failed;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef USE_MPI
    {
        int nprocs=1;
        MPI_Comm_size(LIS_COMM_WORLD,&nprocs);
        if(nprocs!=1) {
            lis_finalize();
            return 0; /* serial-only matrix experiment */
        }
    }
#endif
    failed=run_test();
    lis_finalize();
    return failed;
}
