/* Reproduce transformed Block GLI factory failure without linker wrappers. */
#include "lis.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do { if(!(c)){fprintf(stderr,"FAIL line=%d: %s\n",__LINE__,#c);goto done;} } while(0)
#define OK(c) CHECK((c)==LIS_SUCCESS)
static LIS_INT diag(LIS_MATRIX *A,int metric)
{
    LIS_INT i,rc;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,A)) ||
       (rc=lis_matrix_set_size(*A,0,12)))return rc;
    for(i=0;i<12;i++){
        LIS_SCALAR d=metric?(i==0?1.0:(i==1?1.0/1.01:1.0/(3.0+i*.15))):1.0;
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,*A)))return rc;
    }
    if((rc=lis_matrix_set_type(*A,LIS_MATRIX_CSR)))return rc;
    return lis_matrix_assemble(*A);
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR lam;
    LIS_INT rc,status,j;
    int failure=1,ilutp;
    char opts[256];
    const char *linear;
    if(argc!=3)return 2;
    ilutp=!strcmp(argv[1],"ilutp");
    linear=ilutp?"-i fgmres -p ilut -ilutp 1 -tol 1e-13 -maxiter 500":
        (!strcmp(argv[1],"mfamg")?"-i bicgstab -p mfamg -tol 1e-13 -maxiter 500":
        "-i fgmres -p mfamg_ilutp -tol 1e-13 -maxiter 500");
    snprintf(opts,sizeof opts,"-e gli -etransform %s -shift 0.98 -ss 12 -nev 2 -etol 1e-8 -emaxiter 500 -eprint none",argv[2]);
    OK(lis_initialize(&argc,&argv));
    OK(diag(&A,0));OK(diag(&B,1));
    OK(lis_vector_duplicate(A,&x));
    OK(lis_esolver_create(&e));
    OK(lis_esolver_set_option(opts,e));e->block_size=2;
    OK(lis_esolver_set_linear_options((char*)linear,e));
    for(j=0;j<2;j++){
        OK(lis_vector_set_all(1.0,x));
        rc=lis_gesolve(A,B,x,&lam,e);
        OK(lis_esolver_get_status(e,&status));
        printf("RESULT api=%lld status=%lld repeat=%lld\n",(long long)rc,(long long)status,(long long)j);
        if(ilutp){
            CHECK(rc==LIS_SUCCESS && status==LIS_SUCCESS && e->nev>=2);
            CHECK(fabs((double)lam-1.0)<1e-7 && isfinite((double)e->resid[0]) && e->resid[0]<1e-8);
        }else{
            CHECK(status==LIS_ERR_ILL_ARG || status==LIS_ERR_NOT_IMPLEMENTED);
            CHECK(rc==LIS_SUCCESS || rc==status);
        }
    }
    failure=0;
done:
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    if(B)lis_matrix_destroy(B);
    lis_finalize();
    if(!failure)puts("PASS ownership regression");
    return failure;
}
