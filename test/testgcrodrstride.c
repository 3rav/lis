/*
 * GCRO-DR Ritz cadence regression: default and explicit stride=1 must
 * follow the same numerical path; stride=4 must reduce harmonic updates,
 * while preserving true residuals and cross-solve recycling.
 * Also tests per-solver profile replacement and MPI rank mismatch safety.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include "lis_solver.h"
#include "lis_solver_options.h"
#include "lislib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _COMPLEX
int main(void) { puts("SKIP: GCRO-DR real-only regression"); return 77; }
#else
#define SIDE 24
#define N (SIDE*SIDE)

static int build_matrix(LIS_MATRIX *out)
{
    LIS_MATRIX a=NULL;
    LIS_INT rc,i,is,ie;
    rc=lis_matrix_create(LIS_COMM_WORLD,&a);
    if(rc) return 1;
    if((rc=lis_matrix_set_size(a,0,N)))goto fail;
    if((rc=lis_matrix_get_range(a,&is,&ie)))goto fail;
    for(i=is;i<ie;i++)
    {
        int row=(int)(i/SIDE),col=(int)(i%SIDE);
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.06+0.002*(i%11),a)))
            goto fail;
        if(row>0 && (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-SIDE,-1.0,a)))
            goto fail;
        if(row+1<SIDE && (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+SIDE,-1.0,a)))
            goto fail;
        if(col>0 && (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.16,a)))
            goto fail;
        if(col+1<SIDE && (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.84,a)))
            goto fail;
    }
    if((rc=lis_matrix_set_type(a,LIS_MATRIX_CSR)))goto fail;
    if((rc=lis_matrix_assemble(a)))goto fail;
    *out=a;
    return 0;
fail:
    if(a)lis_matrix_destroy(a);
    return 1;
}

static int check_options(void)
{
    LIS_SOLVER a=NULL,b=NULL;
    LIS_INT changed=0;
    if(lis_solver_create(&a))return 1;
    if(lis_solver_create(&b))goto fail;
    if(lis_gcrodr_get_ritz_stride(a)!=1 ||
       lis_solver_set_option("-gcrodr_ritz_stride 0",a)!=LIS_ERR_ILL_ARG ||
       lis_solver_set_option("-gcrodr_ritz_stride -1",a)!=LIS_ERR_ILL_ARG ||
       lis_solver_set_option("-gcrodr_ritz_stride 1.5",a)!=LIS_ERR_ILL_ARG ||
       lis_solver_set_option("-gcrodr_ritz_stride 9999999999999999999999",a)!=LIS_ERR_ILL_ARG ||
       lis_solver_set_option("-gcrodr_ritz_stride 4",a)!=LIS_SUCCESS ||
       lis_gcrodr_get_ritz_stride(a)!=4)
        goto fail;
    if(lis_fgcrodr_profile_equal(a,b)!=LIS_FALSE)goto fail;
    if(lis_fgcrodr_profile_prepare(b,a)!=LIS_SUCCESS)goto fail;
    lis_fgcrodr_profile_copy(b,a);
    if(lis_gcrodr_get_ritz_stride(b)!=4 ||
       lis_fgcrodr_profile_equal(a,b)!=LIS_TRUE)
        goto fail;
    if(lis_solver_set_option_replace_if_changed(
           "-i gcrodr -p jacobi -gcrodr_ritz_stride 4",b,&changed)!=LIS_SUCCESS ||
       changed!=LIS_TRUE || lis_gcrodr_get_ritz_stride(b)!=4)
        goto fail;
    changed=LIS_TRUE;
    if(lis_solver_set_option_replace_if_changed(
           "-i gcrodr -p jacobi -gcrodr_ritz_stride 4",b,&changed)!=LIS_SUCCESS ||
       changed!=LIS_FALSE)
        goto fail;
    if(lis_solver_set_option_replace_if_changed(
           "-i gcrodr -p jacobi -gcrodr_ritz_stride 1",b,&changed)!=LIS_SUCCESS ||
       changed!=LIS_TRUE || lis_gcrodr_get_ritz_stride(b)!=1)
        goto fail;
    lis_solver_destroy(b);
    lis_solver_destroy(a);
    puts("GCRODR_RITZ_OPTIONS_OK");
    return 0;
fail:
    if(b)lis_solver_destroy(b);
    if(a)lis_solver_destroy(a);
    fputs("GCRODR_RITZ_OPTIONS_FAIL\n",stderr);
    return 1;
}

static int run_case(LIS_MATRIX A,int requested_stride,
                    LIS_INT *iterations,LIS_INT *updates,LIS_INT *uses)
{
    LIS_SOLVER solver=NULL;
    LIS_VECTOR u=NULL,b=NULL,x=NULL,ax=NULL;
    LIS_INT rc=0,i,step;
    LIS_REAL norm=0.,bnorm=0.;
    char opts[300];
    int fail=1;
    *iterations=*updates=*uses=0;
    if((rc=lis_solver_create(&solver)))goto done;
    snprintf(opts,sizeof(opts),
        "-i gcrodr -p jacobi -scale none -print none -restart 8 "
        "-maxiter 1000 -tol 1e-10 -gcrodr_keep 4 "
        "-gcrodr_keep_quality 0.8 -gcrodr_reject_quality 0.2%s",
        requested_stride==0?"":" -gcrodr_ritz_stride 1");
    if((rc=lis_solver_set_option(opts,solver)))goto done;
    if(requested_stride==4 &&
       (rc=lis_solver_set_option("-gcrodr_ritz_stride 4",solver)))goto done;
    if((rc=lis_vector_duplicate(A,&u)))goto done;
    if((rc=lis_vector_duplicate(A,&b)))goto done;
    if((rc=lis_vector_duplicate(A,&x)))goto done;
    if((rc=lis_vector_duplicate(A,&ax)))goto done;
    for(step=0;step<3;step++)
    {
        for(i=0;i<u->n;i++)
        {
            LIS_REAL j=(LIS_REAL)(u->is+i);
            u->value[i]=(LIS_SCALAR)(1.0+0.05*sin(0.023*j+0.2*step)
                                         +0.06*cos(0.041*j-0.1*step));
        }
        if((rc=lis_matvec(A,u,b)))goto done;
        if((rc=lis_vector_set_all(0.0,x)))goto done;
        rc=lis_solve_sequence(A,b,x,solver);
        if(rc!=LIS_SUCCESS || solver->retcode!=LIS_SUCCESS ||
           solver->iter<=0)goto done;
        *iterations += solver->iter;
        if((rc=lis_matvec(A,x,ax)))goto done;
        if((rc=lis_vector_xpay(b,-1.0,ax)))goto done;
        if((rc=lis_vector_nrm2(ax,&norm)))goto done;
        if((rc=lis_vector_nrm2(b,&bnorm)))goto done;
        if(!isfinite(norm) || bnorm<=0.0 || norm/bnorm>1e-8)
            goto done;
    }
    *updates=solver->recycle_updates;
    *uses=solver->recycle_cross_uses;
    if(solver->recycle_dim<=0 || !solver->recycle_valid ||
       *uses<1 || *updates<0)
        goto done;
    printf("GCRODR_RITZ_CASE_OK stride=%d iterations=%lld "
           "updates=%lld cross_uses=%lld residual=%.10e\n",
           requested_stride==0?1:requested_stride,
           (long long)*iterations,(long long)*updates,
           (long long)*uses,(double)(norm/bnorm));
    fail=0;
done:
    if(fail)
        fprintf(stderr,"GCRODR_RITZ_CASE_FAIL stride=%d rc=%lld "
                "status=%lld iterations=%lld updates=%lld uses=%lld\n",
                requested_stride,(long long)rc,
                solver?(long long)solver->retcode:-1LL,
                (long long)*iterations,(long long)*updates,
                (long long)*uses);
    if(ax)lis_vector_destroy(ax);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(u)lis_vector_destroy(u);
    if(solver)lis_solver_destroy(solver);
    return fail;
}

static int check_mpi_mismatch(LIS_MATRIX A)
{
#ifdef USE_MPI
    int rank=0,size=1,fail=1;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR x=NULL,b=NULL;
    LIS_INT rc=0;
    MPI_Comm_rank(A->comm,&rank);
    MPI_Comm_size(A->comm,&size);
    if(size<2)return 0;
    if((rc=lis_solver_create(&solver)))goto done;
    if((rc=lis_solver_set_option(
        "-i gcrodr -p jacobi -restart 8 -maxiter 10 -tol 1e-8 "
        "-gcrodr_keep 4 -print none",solver)))goto done;
    if((rc=lis_solver_set_option(
        rank==0 ? "-gcrodr_ritz_stride 1" : "-gcrodr_ritz_stride 4",
        solver)))goto done;
    if((rc=lis_vector_duplicate(A,&b)))goto done;
    if((rc=lis_vector_duplicate(A,&x)))goto done;
    if((rc=lis_vector_set_all(1.0,b)))goto done;
    if((rc=lis_vector_set_all(0.0,x)))goto done;
    rc=lis_solve_sequence(A,b,x,solver);
    if((rc!=LIS_SUCCESS && rc!=LIS_ERR_ILL_ARG) ||
       solver->retcode!=LIS_ERR_ILL_ARG)
        goto done;
    fail=0;
    if(rank==0)puts("GCRODR_MPI_STRIDE_MISMATCH_REJECTED_OK");
done:
    if(fail)fprintf(stderr,"GCRODR_MPI_STRIDE_MISMATCH_FAIL rank=%d rc=%lld status=%lld\n",
        rank,(long long)rc,solver?(long long)solver->retcode:-1LL);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(solver)lis_solver_destroy(solver);
    return fail;
#else
    (void)A;
    return 0;
#endif
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_INT it[3],updates[3],uses[3];
    int fail=1;
    if(lis_initialize(&argc,&argv))return 1;
    if(build_matrix(&A))goto done;
    if(check_options())goto done;
    if(run_case(A,0,&it[0],&updates[0],&uses[0]))goto done;
    if(run_case(A,1,&it[1],&updates[1],&uses[1]))goto done;
    if(run_case(A,4,&it[2],&updates[2],&uses[2]))goto done;
    if(it[0]!=it[1] || updates[0]!=updates[1] ||
       uses[0]!=uses[1] || updates[2]>=updates[0])
    {
        fprintf(stderr,"GCRODR_RITZ_COUNTER_FAIL base=%lld explicit=%lld stride4=%lld\n",
                (long long)updates[0],(long long)updates[1],(long long)updates[2]);
        goto done;
    }
    if(check_mpi_mismatch(A))goto done;
    puts("GCRODR_RITZ_STRIDE_PASS");
    fail=0;
done:
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return fail;
}
#endif
