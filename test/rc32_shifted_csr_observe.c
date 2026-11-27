/* Selective shifted-CSR regression. Adapted from recovered RC32 observation driver.
 * No near-nullspace prototype, no LU experimental switches.
 * GNU linker wrapping observes the actual factory inputs without changing them.
 */
#include "lis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL line=%d: %s\n",__LINE__,#c); goto done; } } while (0)
#define OK(c) CHECK((c)==LIS_SUCCESS)
static int calls, errors, expected_type, expected_precon, patched;
extern LIS_INT __real_lis_esolver_rebuild_precon(LIS_MATRIX,LIS_SOLVER,LIS_PRECON *);
LIS_INT __wrap_lis_esolver_rebuild_precon(LIS_MATRIX A,LIS_SOLVER s,LIS_PRECON *p)
{
    ++calls;
    if (s->options[LIS_OPTIONS_PRECON]!=expected_precon || A->matrix_type!=expected_type)
        ++errors;
    printf("OBS rebuild precon=%lld matrix=%lld expected=%d\n",
           (long long)s->options[LIS_OPTIONS_PRECON],(long long)A->matrix_type,expected_type);
    return __real_lis_esolver_rebuild_precon(A,s,p);
}
static LIS_INT diag(LIS_MATRIX *A,int n,int metric)
{
    LIS_INT rc,i;
    if ((rc=lis_matrix_create(LIS_COMM_WORLD,A)) ||
        (rc=lis_matrix_set_size(*A,0,n))) return rc;
    for (i=0;i<n;i++) {
        double d=metric ? (i==0?1.0:(i==1?1.0/1.01:1.0/(3.0+i*.15))) : 1.0;
        if ((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,*A))) return rc;
    }
    if ((rc=lis_matrix_set_type(*A,LIS_MATRIX_CSR))) return rc;
    return lis_matrix_assemble(*A);
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,B=NULL;
    LIS_ESOLVER e=NULL;
    LIS_VECTOR x=NULL,ax=NULL,bx=NULL;
    LIS_SCALAR lam;
    LIS_REAL nr,na,nb;
    LIS_INT rc,status,i,j;
    int fail=1;
    char opts[512];
    const char *pre,*transform,*linear;
    double t,elapsed=0,worst=0;
    if (argc!=4) return 2;
    pre=argv[1]; transform=argv[2]; patched=!strcmp(argv[3],"patched");
    if (!strcmp(pre,"none")) {expected_precon=LIS_PRECON_TYPE_NONE; linear="-i bicg -p none";}
    else if (!strcmp(pre,"jacobi")) {expected_precon=LIS_PRECON_TYPE_JACOBI; linear="-i bicg -p jacobi";}
    else if (!strcmp(pre,"mfamg")) {expected_precon=LIS_PRECON_TYPE_MFAMG; linear="-i bicgstab -p mfamg";}
    else if (!strcmp(pre,"mfamg_ilutp")) {expected_precon=LIS_PRECON_TYPE_MFAMG_ILUTP; linear="-i fgmres -p mfamg_ilutp";}
    else if (!strcmp(pre,"ilutp")) {expected_precon=LIS_PRECON_TYPE_ILUT; linear="-i fgmres -p ilut -ilutp 1";}
    else return 2;
    expected_type=(!strcmp(pre,"none") || !strcmp(pre,"jacobi")) ? LIS_MATRIX_OPERATOR : LIS_MATRIX_CSR;
    OK(lis_initialize(&argc,&argv));
    OK(diag(&A,12,0)); OK(diag(&B,12,1));
    OK(lis_vector_duplicate(A,&x)); OK(lis_vector_duplicate(A,&ax)); OK(lis_vector_duplicate(A,&bx));
    OK(lis_esolver_create(&e));
    snprintf(opts,sizeof opts,"-e gli -etransform %s -shift 0.98 -ss 12 -nev 2 -etol 1e-8 -emaxiter 500 -eprint none",transform);
    OK(lis_esolver_set_option(opts,e)); e->block_size=2;
    snprintf(opts,sizeof opts,"%s -tol 1e-13 -maxiter 500",linear);
    OK(lis_esolver_set_linear_options(opts,e));
    for (j=0;j<2;j++) {
        OK(lis_vector_set_all(1.0,x));
        t=lis_wtime(); rc=lis_gesolve(A,B,x,&lam,e); elapsed+=lis_wtime()-t;
        OK(lis_esolver_get_status(e,&status));
        printf("RESULT pre=%s transform=%s repeat=%lld api=%lld status=%lld nev=%lld\n",pre,transform,(long long)j,(long long)rc,(long long)status,(long long)e->nev);
        CHECK(rc==LIS_SUCCESS && status==LIS_SUCCESS && e->nev>=2);
        for (i=0;i<2;i++) {
            OK(lis_esolver_get_specific_evalue(e,i,&lam));
            CHECK(isfinite((double)lam) && fabs((double)lam-(i?1.01:1.0))<1e-7);
            OK(lis_esolver_get_specific_evector(e,i,x));
            OK(lis_matvec(A,x,ax)); OK(lis_matvec(B,x,bx));
            OK(lis_vector_nrm2(ax,&na)); OK(lis_vector_nrm2(bx,&nb));
            OK(lis_vector_axpy(-lam,bx,ax)); OK(lis_vector_nrm2(ax,&nr));
            CHECK(na+fabs((double)lam)*nb>0);
            nr/=na+fabs((double)lam)*nb;
            CHECK(isfinite((double)nr) && nr<1e-8);
            if(nr>worst)worst=nr;
        }
    }
    if (patched) CHECK(calls>=2 && errors==0);
    printf("PASS pre=%s transform=%s calls=%d wall=%.9g max_physical_residual=%.17g\n",pre,transform,calls,elapsed,worst);
    fail=0;
done:
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(ax)lis_vector_destroy(ax);
    if(bx)lis_vector_destroy(bx);
    if(A)lis_matrix_destroy(A);
    if(B)lis_matrix_destroy(B);
    lis_finalize();
    return fail;
}
