#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Nonsymmetric, nonsingular matrix with a zero (0,0) diagonal:
   column pivoting is required for an unregularized exact ILUT factor. */
#define N 4
static const LIS_SCALAR a[N][N] = {
    {0.0, 4.0, 1.0, 0.0},
    {3.0, 0.0,-1.0, 0.0},
    {0.0, 2.0, 0.0, 5.0},
    {1.0, 0.0, 2.0, 1.0}
};
static const LIS_SCALAR truth[N] = {1.25,-0.5,2.0,0.75};

static int run_one(LIS_INT rcm, LIS_INT scale)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL, bt=NULL, x=NULL;
    LIS_INT err,i,j;
    LIS_REAL diff;
    int failed=1,swaps=0;
    char options[512];

    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto cleanup;
    err=lis_matrix_set_size(A,0,N);if(err)goto cleanup;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            if(a[i][j]!=0.0)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,j,a[i][j],A);
                if(err)goto cleanup;
            }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto cleanup;
    err=lis_matrix_assemble(A);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&bt);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        b->value[i]=bt->value[i]=0.0;
        for(j=0;j<N;j++)
        {
            b->value[i]+=a[i][j]*truth[j];
            bt->value[i]+=a[j][i]*truth[j];
        }
    }

    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    snprintf(options,sizeof(options),
        "-i fgmres -p ilut -ilutp 1 -ilutp_threshold 1 "
        "-ilut_rcm %d -ilut_row_scale %d "
        "-ilut_drop 0 -ilut_rate 100 "
        "-ilu_pivot_tol 0 -scale none -omp_num_threads 1",
        (int)rcm,(int)scale);
    err=lis_solver_set_option(options,solver);if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);if(err){pc=NULL;goto cleanup;}
    if(!pc->ilutp_colperm || !pc->ilut_perm || !pc->D)
        goto cleanup;
    if(!rcm && !scale)
        for(i=0;i<N;i++)
            if(pc->ilutp_colperm[i]!=i) swaps++;
    if(!rcm && !scale && !swaps)
        goto cleanup;

    solver->precon=pc;
    err=lis_psolve_ilut_csr(solver,b,x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        diff=fabs(x->value[i]-truth[i]);
        if(!isfinite((double)diff) || diff>1e-9)
        {
            fprintf(stderr,
                "ILUTP forward mismatch rcm=%d scale=%d i=%d diff=%g\n",
                (int)rcm,(int)scale,(int)i,(double)diff);
            goto cleanup;
        }
    }

    err=lis_psolveh_ilut_csr(solver,bt,x);if(err)goto cleanup;
    for(i=0;i<N;i++)
    {
        diff=fabs(x->value[i]-truth[i]);
        if(!isfinite((double)diff) || diff>1e-9)
        {
            fprintf(stderr,
                "ILUTP adjoint mismatch rcm=%d scale=%d i=%d diff=%g\n",
                (int)rcm,(int)scale,(int)i,(double)diff);
            goto cleanup;
        }
    }
    failed=0;

cleanup:
    if(failed)
        fprintf(stderr,"ILUTP test failed (rcm=%d scale=%d err=%d)\n",
            (int)rcm,(int)scale,(int)err);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(bt)lis_vector_destroy(bt);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int failed=0;
    LIS_INT err;
    LIS_SOLVER solver=NULL;
    err=lis_initialize(&argc,&argv);if(err)return 1;
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    failed|=run_one(0,0);
    failed|=run_one(1,0);
    failed|=run_one(0,1);
    failed|=run_one(1,1);
    err=lis_solver_create(&solver);
    if(err)failed=1;
    else
    {
        if(lis_solver_set_option("-p ilut -ilutp 1 -ilutp_threshold 0.5",
                                  solver)!=LIS_SUCCESS ||
           !solver->ilutp_enabled ||
           fabs(solver->ilutp_threshold-0.5)>1e-12 ||
           lis_solver_set_option("-ilutp_threshold 2.0",solver)==LIS_SUCCESS)
            failed=1;
        lis_solver_destroy(solver);
    }
    lis_finalize();
    if(!failed)printf("testilutp: PASS\n");
    return failed;
}
