/*
 * CalculiX 2.19 storage ABI -> LIS_MATRIX_USER integration validation.
 * The build includes the unmodified CalculiX.h from 3rav/CalculiX at
 * commit 84d5dd9ca1bdaf6ad6e92c5c0d50b8f1965975ba.
 *
 * CalculiX ccx/src/spooles.c (inputformat 0/1) defines icol[j] as
 * the number of strict-lower entries in column j, irow[k] as a
 * one-based row number, ad[j] as diagonal, au[0:nzs] as lower
 * and au[nzs:2*nzs] as opposite upper values for nonsymmetric storage.
 *
 * This is a real-header and real-storage-ABI integration test, NOT
 * proof that the unmodified ccx executable calls LIS. The fixture is
 * intentionally constructed, not extracted from an .inp run.
 */
#include "lis.h"
#include "lis_matrix_split_view.h"
#include "lis_matrix_materialize.h"
#include "CalculiX.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CCX_N 4
#define CCX_NNZ 5

typedef struct {
    ITG icol[CCX_N];
    ITG irow[CCX_NNZ];
    double ad[CCX_N];
    double au[2*CCX_NNZ];
} CCX_STORAGE;

static void ccx_reference(const CCX_STORAGE *ccx, int nonsymmetric,
                          int transpose, const double x[CCX_N],
                          double y[CCX_N])
{
    ITG j,k,p=0;
    for(j=0;j<CCX_N;j++) y[j]=ccx->ad[j]*x[j];
    for(j=0;j<CCX_N;j++) {
        for(k=0;k<ccx->icol[j];k++,p++) {
            const ITG i=ccx->irow[p]-1; /* exact ccx FORTRAN -> C origin */
            const double lo=ccx->au[p];
            const double up=nonsymmetric?ccx->au[CCX_NNZ+p]:lo;
            if(!transpose) {
                y[i]+=lo*x[j];
                y[j]+=up*x[i];
            } else {
                y[i]+=up*x[j];
                y[j]+=lo*x[i];
            }
        }
    }
}

static int close_vec(LIS_VECTOR x,const double ref[CCX_N],double tol)
{
    LIS_INT i;
    for(i=0;i<CCX_N;i++) {
        if(!isfinite((double)x->value[i]) ||
           fabs((double)x->value[i]-ref[i])>tol) {
            fprintf(stderr,"vector[%lld]: actual=%.17g expected=%.17g\\n",
                    (long long)i,(double)x->value[i],ref[i]);
            return 0;
        }
    }
    return 1;
}

static int solve_sequence(const char *label,CCX_STORAGE *ccx,int nonsymmetric,
                          LIS_MATRIX A,LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view,
                          const char *options,LIS_INT want_policy)
{
    const double truth[CCX_N]={1.25,-1.5,0.5,2.0};
    const double probe[CCX_N]={-0.75,1.0,2.25,-1.5};
    double reference[CCX_N],ref_transpose[CCX_N];
    LIS_MATRIX_AUTO_SESSION session=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR b=NULL,x=NULL,v=NULL,y=NULL;
    LIS_INT step,i,err,status,policy=-99,action=-99;
    int failed=1;
    err=lis_matrix_auto_session_create(&session);if(err)goto done;
    err=lis_matrix_auto_session_get_solver(session,&solver);if(err)goto done;
    err=lis_solver_set_option((char *)options,solver);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&v);if(err)goto done;
    err=lis_vector_duplicate(A,&y);if(err)goto done;

    for(step=0;step<2;step++) {
        /* Emulate a new Newton tangent: same sparsity, new borrowed values. */
        if(step==1) {
            ccx->ad[0]+=0.75;
            ccx->ad[3]+=0.5;
            ccx->au[0]-=0.125;
            if(nonsymmetric)ccx->au[CCX_NNZ+0]+=0.2;
            err=lis_matrix_split_triangular_shared_view_set_values_real64(
                view,ccx->ad,ccx->au,
                nonsymmetric?ccx->au+CCX_NNZ:NULL);
            if(err)goto done;
        }
        /* Independent CCX column-storage operator, not LIS SpMV, creates b. */
        ccx_reference(ccx,nonsymmetric,0,truth,reference);
        ccx_reference(ccx,nonsymmetric,1,probe,ref_transpose);
        for(i=0;i<CCX_N;i++) {
            b->value[i]=(LIS_SCALAR)reference[i];
            v->value[i]=(LIS_SCALAR)probe[i];
        }
        err=lis_matvech(A,v,y);
        if(err || !close_vec(y,ref_transpose,1e-11))goto done;
        err=lis_matvec(A,v,y);
        ccx_reference(ccx,nonsymmetric,0,probe,ref_transpose);
        if(err || !close_vec(y,ref_transpose,1e-11))goto done;

        err=lis_vector_set_all((LIS_SCALAR)0,x);if(err)goto done;
        err=lis_matrix_auto_session_solve(session,A,b,x,&policy,&action);
        if(err) {
            fprintf(stderr,"%s step=%lld session error=%lld\\n",
                    label,(long long)step,(long long)err);
            goto done;
        }
        err=lis_solver_get_status(solver,&status);
        if(err || status!=LIS_SUCCESS || policy!=want_policy) {
            fprintf(stderr,"%s step=%lld status=%lld policy=%lld want=%lld\\n",
                    label,(long long)step,(long long)status,
                    (long long)policy,(long long)want_policy);
            goto done;
        }
        if(want_policy==LIS_MATRIX_POLICY_DIRECT &&
           action!=LIS_MATRIX_POLICY_ACTION_DIRECT)goto done;
        if(want_policy==LIS_MATRIX_POLICY_CACHED_CSR &&
           action==LIS_MATRIX_POLICY_ACTION_DIRECT)goto done;
        if(!close_vec(x,truth,1e-8)) goto done;

        /* Independent Ax-b residual check; catches false solver SUCCESS. */
        err=lis_matvec(A,x,y);if(err)goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1,b,y);if(err)goto done;
        {
            LIS_REAL norm=0,rhsnorm=0;
            err=lis_vector_nrm2(y,&norm);if(err)goto done;
            err=lis_vector_nrm2(b,&rhsnorm);if(err)goto done;
            if(!isfinite((double)norm) || norm>1e-9*(1+rhsnorm)) {
                fprintf(stderr,"%s step=%lld residual=%g rhs=%g\\n",
                        label,(long long)step,(double)norm,(double)rhsnorm);
                goto done;
            }
            printf("CCX_USER case=%s step=%lld policy=%s action=%lld residual=%.3e PASS\\n",
                   label,(long long)step,
                   policy==LIS_MATRIX_POLICY_DIRECT?"DIRECT":"CACHED_CSR",
                   (long long)action,(double)norm);
        }
    }
    failed=0;
done:
    if(failed)fprintf(stderr,"CCX_USER FAIL case=%s err=%lld\\n",
                      label,(long long)err);
    if(session)lis_matrix_auto_session_destroy(session);
    if(y)lis_vector_destroy(y);
    if(v)lis_vector_destroy(v);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    return failed;
}

int main(int argc,char **argv)
{
    CCX_STORAGE ccx_sym={
        {2,2,1,0}, {2,4,3,4,4},
        {9.,10.,11.,12.},
        {-0.5,0.2,-0.7,0.3,-0.4, 0,0,0,0,0}
    };
    CCX_STORAGE ccx_unsym={
        {2,2,1,0}, {2,4,3,4,4},
        {9.,10.,11.,12.},
        {-0.5,0.2,-0.7,0.3,-0.4,
         -0.2,-0.1,-0.8,0.1,-0.6}
    };
    int64_t ptr[CCX_N+1],row[CCX_NNZ];
    ITG k,j,total=0;
    LIS_MATRIX A=NULL,B=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW vA=NULL,vB=NULL;
    LIS_INT err;
    int failed=1;
    err=lis_initialize(&argc,&argv);if(err)return 1;
    /* The original CCX icol[] is a vector of per-column counts, not ptr[]. */
    ptr[0]=1;
    for(j=0;j<CCX_N;j++) {
        if(ccx_sym.icol[j]<0)goto done;
        total+=ccx_sym.icol[j];
        ptr[j+1]=1+(int64_t)total;
    }
    if(total!=CCX_NNZ)goto done;
    for(k=0;k<CCX_NNZ;k++)
        row[k]=(int64_t)ccx_sym.irow[k];

    err=lis_matrix_create_split_triangular_topology_i64(
        CCX_N,CCX_NNZ,ptr,row,1,&topology);
    if(err)goto done;
    err=lis_matrix_create_split_triangular_shared_view_real64(
        LIS_COMM_WORLD,topology,ccx_sym.ad,ccx_sym.au,NULL,&A,&vA);
    if(err)goto done;
    err=lis_matrix_create_split_triangular_shared_view_real64(
        LIS_COMM_WORLD,topology,ccx_unsym.ad,
        ccx_unsym.au,ccx_unsym.au+CCX_NNZ,&B,&vB);
    if(err)goto done;
    lis_matrix_split_triangular_topology_destroy(topology);
    topology=NULL;

    if(solve_sequence("SPD-CG-JACOBI",&ccx_sym,0,A,vA,
           "-i cg -p jacobi -tol 1e-11 -maxiter 200 -print none -scale none",
           LIS_MATRIX_POLICY_DIRECT))goto done;
    if(solve_sequence("UNSYM-GMRES-ILUT",&ccx_unsym,1,B,vB,
           "-i gmres -p ilut -restart 20 -tol 1e-11 -maxiter 200 -print none -scale none",
           LIS_MATRIX_POLICY_DIRECT))goto done;
    if(solve_sequence("UNSYM-GMRES-SSOR",&ccx_unsym,1,B,vB,
           "-i gmres -p ssor -restart 20 -tol 1e-11 -maxiter 200 -print none -scale none",
           LIS_MATRIX_POLICY_CACHED_CSR))goto done;
    failed=0;
done:
    if(A)lis_matrix_destroy(A);
    if(B)lis_matrix_destroy(B);
    if(vA)lis_matrix_split_triangular_shared_view_destroy(vA);
    if(vB)lis_matrix_split_triangular_shared_view_destroy(vB);
    if(topology)lis_matrix_split_triangular_topology_destroy(topology);
    lis_finalize();
    printf("CCX_USER_VALIDATION=%s\\n",failed?"FAIL":"PASS");
    return failed;
}
