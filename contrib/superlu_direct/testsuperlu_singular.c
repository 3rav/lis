#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include "lislib.h"
#include "lis_superlu_direct.h"

/* Verify that both factorization failure paths release partial state. */
static int run_case(int mode)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=1,i;
    const int n=3;
    double a[3][3]={{4.,1.,0.},{2.,2.,0.},{0.,0.,5.}};
    if(mode==0) for(i=0;i<n;i++) a[1][i]=0.; /* zero row */
    if(mode==1) { a[1][0]=8.; a[1][1]=2.; } /* dependent rows */
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto cleanup;
    err=lis_matrix_set_size(A,0,n);if(err)goto cleanup;
    for(i=0;i<n;i++) for(int j=0;j<n;j++) if(a[i][j]) {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,j,a[i][j],A);
        if(err)goto cleanup;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto cleanup;
    err=lis_matrix_assemble(A);if(err)goto cleanup;
    err=lis_solver_create(&solver);if(err)goto cleanup;
    solver->A=A;
    err=lis_solver_set_option("-p superlu -i fgmres -scale none",solver);
    if(err)goto cleanup;
    err=lis_precon_create(solver,&pc);
    /* lis_precon_create destroys pc itself if setup fails. */
    if(err)pc=NULL;
    if(mode<2) {
        if(err!=LIS_BREAKDOWN) {
            fprintf(stderr,"SUPERLU_SINGULAR failure mode=%d status=%d expected=%d\n",
                    mode,(int)err,(int)LIS_BREAKDOWN);
            err=LIS_BREAKDOWN;
        } else err=LIS_SUCCESS;
    } else if(err) {
        fprintf(stderr,"SUPERLU_SINGULAR valid matrix status=%d\n",(int)err);
    }
cleanup:
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(A)lis_matrix_destroy(A);
    return err!=LIS_SUCCESS;
}
int main(int argc,char **argv)
{
    int i,failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_superlu_direct_register())return 1;
    for(i=0;i<100;i++) {
        failed|=run_case(0);
        failed|=run_case(1);
        failed|=run_case(2);
    }
    lis_finalize();
    puts(failed?"SUPERLU_SINGULAR FAIL":
         "SUPERLU_SINGULAR PASS zero_row rank_deficient reuse_100");
    return failed;
}
