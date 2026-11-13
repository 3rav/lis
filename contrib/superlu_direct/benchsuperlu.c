/*
 * Reproducible native LIS+SuperLU solve of an original CCX Matrix Market A,b.
 * Complete LU is used once as a user-registered LIS preconditioner.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "lislib.h"
#include "lis_superlu_direct.h"

static int read_rhs(const char *name,LIS_VECTOR v,LIS_INT n)
{
    FILE *fp=fopen(name,"r");
    char buf[2048],banner[64],obj[64],fmt[64],dtype[64],shape[64];
    long rows,cols,k;
    double val;
    if(!fp)return 1;
    if(!fgets(buf,sizeof(buf),fp) ||
       sscanf(buf,"%63s %63s %63s %63s %63s",banner,obj,fmt,dtype,shape)!=5 ||
       strcmp(banner,"%%MatrixMarket") || strcmp(obj,"matrix") ||
       strcmp(fmt,"array") || strcmp(dtype,"real") || strcmp(shape,"general"))
    { fclose(fp);return 1; }
    do { if(!fgets(buf,sizeof(buf),fp)){fclose(fp);return 1;} }
    while(buf[0]=='%');
    if(sscanf(buf,"%ld %ld",&rows,&cols)!=2 ||
       rows!=(long)n || cols!=1) {fclose(fp);return 1;}
    for(k=0;k<rows;k++)
    {
        if(!fgets(buf,sizeof(buf),fp) ||
           sscanf(buf,"%lf",&val)!=1 || !isfinite(val))
        {fclose(fp);return 1;}
        v->value[k]=(LIS_SCALAR)val;
    }
    fclose(fp);return 0;
}

int main(int argc,char **argv)
{
    LIS_INT err=1,iter=0,rc=0;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,residual=NULL;
    LIS_SOLVER solver=NULL;
    LIS_REAL bn=0,rn=0;
    double sec=0,rel=0;
    clock_t t;
    if(argc!=3)
    {
        fprintf(stderr,"Usage: %s A.mtx b.mtx\n",argv[0]);
        return 2;
    }
    if(lis_initialize(&argc,&argv))return 1;
    err=lis_superlu_direct_register();if(err)goto cleanup;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto cleanup;
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto cleanup;
    err=lis_input_matrix(A,argv[1]);if(err)goto cleanup;
    if(A->nprocs!=1 || A->n!=A->gn)goto cleanup;
    err=lis_vector_duplicate(A,&b);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&x);if(err)goto cleanup;
    err=lis_vector_duplicate(A,&residual);if(err)goto cleanup;
    if(read_rhs(argv[2],b,A->n)){err=LIS_ERR_FILE_IO;goto cleanup;}
    err=lis_vector_nrm2(b,&bn);if(err||!(bn>0))goto cleanup;
    err=lis_solver_create(&solver);if(err)goto cleanup;
    err=lis_solver_set_option("-i fgmres -p superlu -scale none "
        "-tol 1e-10 -maxiter 20 -restart 10 -omp_num_threads 1",solver);
    if(err)goto cleanup;
    err=lis_vector_set_all((LIS_SCALAR)0,x);if(err)goto cleanup;
    t=clock();
    err=lis_solve(A,b,x,solver);
    sec=(double)(clock()-t)/CLOCKS_PER_SEC;
    rc=solver->retcode;
    iter=solver->iter;
    if(!err)
    {
        err=lis_matvec(A,x,residual);if(err)goto cleanup;
        err=lis_vector_axpy((LIS_SCALAR)-1,b,residual);if(err)goto cleanup;
        err=lis_vector_nrm2(residual,&rn);if(err)goto cleanup;
        rel=(double)(rn/bn);
    }
    printf("SUPERLU_LIS_NATIVE n=%lld nnz=%lld iter=%lld status=%lld "
           "true_relres=%.16g total_cpu_s=%.6f\n",
           (long long)A->n,(long long)A->nnz,(long long)iter,
           (long long)rc,rel,sec);
    if(rc!=LIS_SUCCESS || rel>1e-9)err=LIS_BREAKDOWN;
cleanup:
    if(err)fprintf(stderr,"SUPERLU_LIS_NATIVE error=%lld\n",(long long)err);
    if(solver)lis_solver_destroy(solver);
    if(residual)lis_vector_destroy(residual);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return err?1:0;
}
