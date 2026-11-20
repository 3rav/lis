#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include "lis_matrix_materialize.h"
#include <math.h>
#include <stdio.h>

static LIS_INT row_cb(void *unused,LIS_INT row,LIS_INT *nnz,
                      LIS_INT *cols,LIS_SCALAR *vals)
{
    (void)unused;
    if(!nnz || row<0 || row>=3) return LIS_ERR_ILL_ARG;
    if(cols==NULL && vals==NULL){*nnz=1;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<1)return LIS_ERR_ILL_ARG;
    cols[0]=row;
    vals[0]=(LIS_SCALAR)(2+row);
    *nnz=1;
    return LIS_SUCCESS;
}
static int run_solve(LIS_MATRIX_AUTO_SESSION session,LIS_SOLVER solver,
                     LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,
                     LIS_INT expected_policy,LIS_INT expected_action)
{
    LIS_INT err,policy=77,action=88,status,i;
    lis_vector_set_all((LIS_SCALAR)0,x);
    err=lis_matrix_auto_session_solve(session,A,b,x,&policy,&action);
    if(err || policy!=expected_policy || action!=expected_action)
    {
        fprintf(stderr,"AUTO_ROW_ONLY err=%d policy=%d action=%d\n",
                (int)err,(int)policy,(int)action);
        return 1;
    }
    err=lis_solver_get_status(solver,&status);
    if(err || status!=LIS_SUCCESS)return 1;
    for(i=0;i<3;i++)
        if(fabs((double)x->value[i]-1.0)>1e-8)return 1;
    return 0;
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_MATRIX_AUTO_SESSION session=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err,i;
    int failed=1;
    if(lis_initialize(&argc,&argv))return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,3);if(err)goto done;
    err=lis_matrix_set_user(A,NULL,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,row_cb);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<3;i++)b->value[i]=(LIS_SCALAR)(2+i);
    err=lis_matrix_auto_session_create(&session);if(err)goto done;
    err=lis_matrix_auto_session_get_solver(session,&solver);if(err)goto done;
    err=lis_solver_set_option("-i gmres -p none -scale none -restart 5 -maxiter 60 -tol 1e-11 -print none",solver);
    if(err)goto done;
    if(run_solve(session,solver,A,b,x,
        LIS_MATRIX_POLICY_DIRECT,LIS_MATRIX_POLICY_ACTION_DIRECT))goto done;
    err=lis_solver_set_option("-p jacobi",solver);if(err)goto done;
    if(run_solve(session,solver,A,b,x,
        LIS_MATRIX_POLICY_CACHED_CSR,LIS_MATRIX_MATERIALIZE_REBUILD))goto done;
    err=lis_solver_set_option("-p none",solver);if(err)goto done;
    if(run_solve(session,solver,A,b,x,
        LIS_MATRIX_POLICY_DIRECT,LIS_MATRIX_POLICY_ACTION_DIRECT))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"AUTO_ROW_ONLY FAILED\n");
    if(session)lis_matrix_auto_session_destroy(session);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    printf("AUTO_ROW_ONLY %s\n",failed?"FAILED":"PASSED");
    return failed;
}
