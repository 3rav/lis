/* Standalone exercise of the optional CCX->USER checked bridge.
 * Includes the original CCX CalculiX.h, does not modify ccx Makefiles.
 */
#ifndef HAVE_CONFIG_H
#error "Compile with the same LIS config and integer-width macros"
#endif
#include "ccx_lis_user_adapter.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static void reference(ITG n,const ITG *icol,const ITG *row,
                      const double *diag,const double *au,
                      ITG offset,const double *x,double *b)
{
    ITG i,j,k,p=0;
    for(i=0;i<n;i++)b[i]=diag[i]*x[i];
    for(j=0;j<n;j++)for(k=0;k<icol[j];k++,p++){
        ITG irow=row[p]-1;
        b[irow]+=au[p]*x[j];
        b[j]+=au[offset+p]*x[irow];
    }
}
static int check(const char *label,const double *actual,
                 const double *expected,ITG n)
{
    ITG i;
    for(i=0;i<n;i++)if(!isfinite(actual[i])||
                        fabs(actual[i]-expected[i])>1e-8){
        fprintf(stderr,"%s[%lld]=%.17g want %.17g\n",
            label,(long long)i,actual[i],expected[i]);
        return 1;
    }
    return 0;
}
int main(int argc,char **argv)
{
    const ITG n=4,nnz=5,offset=7;
    ITG icol[4]={2,2,1,0},irow[5]={2,4,3,4,4};
    ITG changed_icol[4]={1,2,1,0},changed_irow[4]={2,3,4,4};
    ITG bad_icol[4]={-1,2,1,0};
    ITG huge_icol[4]={2147483647,2,1,0};
    ITG bad_row[5]={1,4,3,4,4};
    ITG unsorted_row[5]={4,2,3,4,4};
    double diag[4]={9,10,11,12};
    double au[12]={-.5,.2,-.7,.3,-.4,0,0,-.2,-.1,-.8,.1,-.6};
    double diag_new[4]={8,9,10,11};
    double au_new[10]={-.2,.3,-.25,-.4,0,-.1,-.5,.15,-.2,0};
    /* Borrowed arrays must remain alive through every error path and
     * eventual context destruction, not only within a local test block. */
    double sym_diag[4]={9,10,11,12};
    double sym_au[10]={-.5,.2,-.7,.3,-.4, -.5,.2,-.7,.3,-.4};
    double sym_rhs[4];
    const ITG empty_col[4]={0,0,0,0};
    const double diagonal_values[4]={2,3,4,5};
    const double exact_diagonal[4]={1,-2,3,-4};
    const double rhs_diagonal[4]={2,-6,12,-20};
    const double truth[4]={1,-1,2,.5};
    double b[4],solution[4]={0,0,0,0},bcopy[4];
    double rr=99;
    CCX_LIS_USER_CONTEXT *ctx=NULL,*bad=NULL,*diagonal=NULL,*direct=NULL,*symmetric=NULL;
    LIS_INT err=0,policy=0,action=0;
    int failed=1;
    err=lis_initialize(&argc,&argv);
    if(err)return 1;

    /* Bad CCX input must fail before materialization or array access. */
#define EXPECT_REJECT(tag, expr) do { \
    err=(expr); \
    if(err==LIS_SUCCESS || bad!=NULL) { \
        fprintf(stderr,"CCX_USER_ADAPTER bad %s unexpectedly succeeded\n",tag); \
        goto done; \
    } \
} while(0)
    EXPECT_REJECT("negative counts",ccx_lis_user_open(
        n,nnz,bad_icol,irow,diag,au,12,offset,1,NULL,&bad));
    EXPECT_REJECT("overflowed counts",ccx_lis_user_open(
        n,nnz,huge_icol,irow,diag,au,12,offset,1,NULL,&bad));
    EXPECT_REJECT("diagonal row index",ccx_lis_user_open(
        n,nnz,icol,bad_row,diag,au,12,offset,1,NULL,&bad));
    EXPECT_REJECT("unsorted rows",ccx_lis_user_open(
        n,nnz,icol,unsorted_row,diag,au,12,offset,1,NULL,&bad));
    EXPECT_REJECT("short au capacity",ccx_lis_user_open(
        n,nnz,icol,irow,diag,au,11,offset,1,NULL,&bad));
    EXPECT_REJECT("overlapping upper",ccx_lis_user_open(
        n,nnz,icol,irow,diag,au,12,1,1,NULL,&bad));
#undef EXPECT_REJECT

    err=ccx_lis_user_open(
        n,nnz,icol,irow,diag,au,12,offset,1,
        "-i gmres -p ssor -restart 10 -tol 1e-11 -maxiter 100 -scale none -print none",
        &ctx);
    if(err||!ctx)goto done;
    reference(n,icol,irow,diag,au,offset,truth,b);
    for(ITG i=0;i<n;i++)bcopy[i]=b[i];
    /* Same native CCX backing arrays, but an ILUT-eligible solver:
     * the bridge must not request materialization or copy numeric values. */
    err=ccx_lis_user_open(
        n,nnz,icol,irow,diag,au,12,offset,1,
        "-i gmres -p ilut -restart 10 -tol 1e-11 -maxiter 100 -scale none -print none",
        &direct);
    if(err || !direct)goto done;
    err=ccx_lis_user_solve(direct,b,solution,&rr,&policy,&action);
    if(err || policy!=LIS_MATRIX_POLICY_DIRECT ||
       action!=LIS_MATRIX_POLICY_ACTION_DIRECT || rr>1e-10 ||
       check("unsymmetric direct",solution,truth,n))goto done;
    ccx_lis_user_close(direct);direct=NULL;
    printf("CCX_ADAPTER unsymmetric ILUT direct USER PASS\n");

    /* Symmetric CCX inputformat=0: the borrowed lower values must be
     * reflected by the generic USER view; no upper array is passed. */
    {
        reference(n,icol,irow,sym_diag,sym_au,5,truth,sym_rhs);
        err=ccx_lis_user_open(n,nnz,icol,irow,sym_diag,sym_au,10,0,0,
            "-i cg -p jacobi -tol 1e-11 -maxiter 100 -scale none -print none",
            &symmetric);
        if(err || !symmetric)goto done;
        err=ccx_lis_user_solve(symmetric,sym_rhs,solution,&rr,
                               &policy,&action);
        if(err || policy!=LIS_MATRIX_POLICY_DIRECT ||
           action!=LIS_MATRIX_POLICY_ACTION_DIRECT || rr>1e-10 ||
           check("symmetric direct",solution,truth,n))goto done;
        ccx_lis_user_close(symmetric);symmetric=NULL;
        printf("CCX_ADAPTER symmetric CG/Jacobi direct USER PASS\n");
    }

    err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
    if(err||policy!=LIS_MATRIX_POLICY_CACHED_CSR||
       action!=LIS_MATRIX_MATERIALIZE_REBUILD||
       rr>1e-10||check("initial",solution,truth,n)||
       check("RHS unchanged",b,bcopy,n))goto done;
    printf("CCX_ADAPTER first CSR_REBUILD residual=%.3e PASS\n",rr);
    err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
    if(err||action!=LIS_MATRIX_MATERIALIZE_REUSE||rr>1e-10||
       check("reuse",solution,truth,n))goto done;
    printf("CCX_ADAPTER no-change CSR_REUSE PASS\n");

    /* Reject a non-finite RHS without exposing an incomplete solution.
     * A failed solve must also leave policy/action and metrics untouched. */
    {
        double bad_rhs[4]={b[0],b[1],NAN,b[3]};
        double untouched[4]={991,992,993,994};
        double untouched_metric=123;
        LIS_INT untouched_policy=55,untouched_action=66;
        err=ccx_lis_user_solve(ctx,bad_rhs,untouched,&untouched_metric,
                               &untouched_policy,&untouched_action);
        if(err!=LIS_ERR_ILL_ARG || untouched_metric!=123 ||
           untouched_policy!=55 || untouched_action!=66 ||
           untouched[0]!=991 || untouched[1]!=992 ||
           untouched[2]!=993 || untouched[3]!=994)goto done;
        err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
        if(err || action!=LIS_MATRIX_MATERIALIZE_REUSE ||
           check("recovered RHS",solution,truth,n))goto done;
    }
    printf("CCX_ADAPTER invalid RHS atomic rejection PASS\n");

    /* Reject invalid numeric candidates before publishing a new borrowed
     * matrix generation. This protects the previous CSR and preconditioner
     * because the previous borrowed storage is not modified. */
    {
        double invalid_diag[4];
        double invalid_au[12];
        memcpy(invalid_diag,diag,sizeof(diag));
        memcpy(invalid_au,au,sizeof(au));
        invalid_diag[2]=NAN;
        err=ccx_lis_user_update(ctx,invalid_diag,au,12,offset);
        if(err!=LIS_ERR_ILL_ARG)goto done;
        invalid_diag[2]=diag[2];
        invalid_au[offset+3]=INFINITY;
        err=ccx_lis_user_update(ctx,diag,invalid_au,12,offset);
        if(err!=LIS_ERR_ILL_ARG)goto done;
        err=ccx_lis_user_rebind(ctx,nnz,icol,irow,
                                invalid_diag,invalid_au,12,offset);
        if(err!=LIS_ERR_ILL_ARG)goto done;
        err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
        if(err||action!=LIS_MATRIX_MATERIALIZE_REUSE||
           check("rejected nonfinite update",solution,truth,n))goto done;
    }
    printf("CCX_ADAPTER nonfinite tangent rejection and rollback PASS\n");

    /* Invalid numeric update must preserve the old operator/cache. */
    err=ccx_lis_user_update(ctx,diag,au,11,offset);
    if(err==LIS_SUCCESS)goto done;
    err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
    if(err||action!=LIS_MATRIX_MATERIALIZE_REUSE||
       check("failed-update rollback",solution,truth,n))goto done;

    diag[0]+=1;
    au[offset+1]-=.25;
    err=ccx_lis_user_update(ctx,diag,au,12,offset);
    if(err)goto done;
    reference(n,icol,irow,diag,au,offset,truth,b);
    err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
    if(err||action!=LIS_MATRIX_MATERIALIZE_VALUES||rr>1e-10||
       check("numeric refresh",solution,truth,n))goto done;
    printf("CCX_ADAPTER Newton numeric CSR_VALUES PASS\n");

    err=ccx_lis_user_rebind(ctx,nnz,bad_icol,irow,
                            diag,au,12,offset);
    if(err==LIS_SUCCESS)goto done;
    err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
    if(err||action!=LIS_MATRIX_MATERIALIZE_REUSE||
       check("failed-topology rollback",solution,truth,n))goto done;

    /* Replace sparse topology with a genuinely different nnz/count map. */
    err=ccx_lis_user_rebind(ctx,4,changed_icol,changed_irow,
                            diag_new,au_new,10,5);
    if(err)goto done;
    reference(n,changed_icol,changed_irow,diag_new,au_new,5,truth,b);
    err=ccx_lis_user_solve(ctx,b,solution,&rr,&policy,&action);
    if(err||action!=LIS_MATRIX_MATERIALIZE_REBUILD||rr>1e-10||
       check("contact topology",solution,truth,n))goto done;
    printf("CCX_ADAPTER contact-like structure CSR_REBUILD PASS\n");

    /* Diagonal-only stiffness systems are legal, with nnz=0 and au=NULL.
     * Exercise an entirely direct route without a CSR copy. */
    {
        err=ccx_lis_user_open(n,0,empty_col,NULL,diagonal_values,NULL,0,0,0,
            "-i cg -p jacobi -tol 1e-12 -maxiter 100 -scale none -print none",
            &diagonal);
        if(err || !diagonal)goto done;
        err=ccx_lis_user_solve(diagonal,rhs_diagonal,solution,&rr,
                               &policy,&action);
        if(err || policy!=LIS_MATRIX_POLICY_DIRECT ||
           action!=LIS_MATRIX_POLICY_ACTION_DIRECT || rr>1e-10 ||
           check("diagonal only",solution,exact_diagonal,n))goto done;
        printf("CCX_ADAPTER diagonal-only direct USER PASS\n");
    }
    failed=0;
done:
    if(failed)fprintf(stderr,
        "CCX_USER_ADAPTER FAIL err=%lld policy=%lld action=%lld\n",
        (long long)err,(long long)policy,(long long)action);
    ccx_lis_user_close(ctx);
    ccx_lis_user_close(direct);
    ccx_lis_user_close(symmetric);
    ccx_lis_user_close(diagonal);
    ccx_lis_user_close(bad);
    lis_finalize();
    printf("CCX_USER_ADAPTER=%s\n",failed?"FAIL":"PASS");
    return failed;
}
