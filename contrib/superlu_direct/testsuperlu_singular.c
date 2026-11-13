#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include "lislib.h"
#include "lis_superlu_direct.h"
/* Singular rank-deficient matrix: setup must fail and release partial LU. */
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err;
    err=lis_initialize(&argc,&argv);if(err)return 1;
    err=lis_superlu_direct_register();if(err)return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)return 1;
    err=lis_matrix_set_size(A,0,3);if(err)return 1;
    lis_matrix_set_value(LIS_INS_VALUE,0,0,1,A);
    lis_matrix_set_value(LIS_INS_VALUE,0,1,2,A);
    lis_matrix_set_value(LIS_INS_VALUE,1,0,2,A);
    lis_matrix_set_value(LIS_INS_VALUE,1,1,4,A);
    lis_matrix_set_value(LIS_INS_VALUE,2,2,1,A);
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)return 1;
    err=lis_matrix_assemble(A);if(err)return 1;
    err=lis_solver_create(&solver);if(err)return 1;
    solver->A=A;
    err=lis_solver_set_option("-p superlu -scale none",solver);
    if(err)return 1;
    err=lis_precon_create(solver,&pc);
    printf("SUPERLU_SINGULAR status=%d expected=%d\n",
           (int)err,(int)LIS_BREAKDOWN);
    /* On failure lis_precon_create already calls the user destructor. */
    if(err==LIS_SUCCESS && pc) lis_precon_destroy(pc);
    solver->precon=NULL;
    lis_solver_destroy(solver);
    lis_matrix_destroy(A);
    lis_finalize();
    return err==LIS_BREAKDOWN ? 0:1;
}
