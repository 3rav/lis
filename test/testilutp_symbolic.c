/* RC32 audit: public sequence lifecycle, changed pattern, fresh-factor equivalence.
 * Test-only linker wrappers observe calls and preserve original behavior.
 */
#include "lis.h"
#include "lis_precon.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 6
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL line=%d: %s\n",__LINE__,#c); goto done; } } while (0)
#define OK(c) CHECK((c)==LIS_SUCCESS)
static int creates,refreshes,refresh_errors;
extern LIS_INT __real_lis_precon_create(LIS_SOLVER,LIS_PRECON *);
extern LIS_INT __real_lis_precon_update_ilut_frozen(LIS_SOLVER,LIS_PRECON);
LIS_INT __wrap_lis_precon_create(LIS_SOLVER s,LIS_PRECON *p)
{
    creates++;
    return __real_lis_precon_create(s,p);
}
LIS_INT __wrap_lis_precon_update_ilut_frozen(LIS_SOLVER s,LIS_PRECON p)
{
    LIS_INT err;
    refreshes++;
    err=__real_lis_precon_update_ilut_frozen(s,p);
    if(err)refresh_errors++;
    return err;
}
static LIS_INT make_matrix(LIS_MATRIX *A,int pattern)
{
    LIS_INT err,i,j;
    if((err=lis_matrix_create(LIS_COMM_WORLD,A)) ||
       (err=lis_matrix_set_size(*A,0,N)))return err;
    for(i=0;i<N;i++)for(j=0;j<N;j++) {
        /* Symmetric structural pattern, nonsymmetric numerical values. */
        if(i==j || (i+j==(pattern?N-1:1)) ||
           (!pattern && (i+j==5 || i+j==9))) {
            /* pattern 0: pairs 0-1, 0-5, 1-4, 2-3, 4-5;
             * pattern 1: pairs 0-5, 1-4, 2-3. */
            if((err=lis_matrix_set_value(LIS_INS_VALUE,i,j,
                 i==j?8.0+i:0.1*(i+1)+0.2*(j+1),*A)))return err;
        }
    }
    if((err=lis_matrix_set_type(*A,LIS_MATRIX_CSR)))return err;
    return lis_matrix_assemble(*A);
}
static LIS_INT setup(LIS_SOLVER *s)
{
    LIS_INT err;
    if((err=lis_solver_create(s)))return err;
    if((err=lis_solver_set_option("-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 10000 -ilutp_symbolic 1 -ilutp_md 1 -ilutp_reuse 1 -tol 1e-12 -maxiter 40 -restart 20 -print none",*s)))return err;
    return lis_solver_set_sequence_precon_refresh_policy(*s,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN);
}
static int solve_compare(LIS_MATRIX A,LIS_SOLVER s,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR truth,LIS_VECTOR scratch)
{
    LIS_SOLVER fresh=NULL;
    LIS_VECTOR y=NULL;
    LIS_REAL r,bnorm,errnorm,tnorm,delta;
    LIS_INT status,i,j;
    LIS_PRECON p,q;
    void *user=NULL;
    int failed=1;
    if(lis_matvec(A,truth,b)||lis_vector_set_all(0,x))goto done;
    if(lis_solve_sequence(A,b,x,s)||lis_solver_get_status(s,&status)||status!=LIS_SUCCESS)goto done;
    p=s->sequence_precon;
    if(!p||!p->ilutp_symbolic_factor||!p->ilutp_symbolic_valid||!p->ilutp_reuse)goto done;
    if(lis_precon_set_user_data(p,&failed)||lis_precon_get_user_data(p,&user)||user!=&failed||lis_precon_set_user_data(p,NULL))goto done;
    if(lis_matvec(A,x,scratch)||lis_vector_axpy(-1,b,scratch)||lis_vector_nrm2(scratch,&r)||lis_vector_nrm2(b,&bnorm))goto done;
    if(!(bnorm>0)||!isfinite(r)||r/bnorm>1e-11)goto done;
    if(lis_vector_copy(x,scratch)||lis_vector_axpy(-1,truth,scratch)||lis_vector_nrm2(scratch,&errnorm)||lis_vector_nrm2(truth,&tnorm))goto done;
    if(!(tnorm>0)||!isfinite(errnorm)||errnorm/tnorm>1e-10)goto done;
    if(setup(&fresh)||lis_vector_duplicate(A,&y)||lis_vector_set_all(0,y))goto done;
    if(lis_solve_sequence(A,b,y,fresh)||lis_solver_get_status(fresh,&status)||status!=LIS_SUCCESS)goto done;
    q=fresh->sequence_precon;
    if(!q||!q->ilutp_symbolic_factor||!q->ilutp_symbolic_valid)goto done;
    for(i=0;i<A->n;i++) {
        if(p->ilut_perm[i]!=q->ilut_perm[i] || !isfinite((double)p->D->value[i]) || !isfinite((double)q->D->value[i]) || fabs(p->D->value[i]-q->D->value[i])>1e-12 || p->L->nnz[i]!=q->L->nnz[i] || p->U->nnz[i]!=q->U->nnz[i])goto done;
        for(j=0;j<p->L->nnz[i];j++)if(p->L->index[i][j]!=q->L->index[i][j] || !isfinite((double)p->L->value[i][j]) || !isfinite((double)q->L->value[i][j]) || fabs(p->L->value[i][j]-q->L->value[i][j])>1e-12)goto done;
        for(j=0;j<p->U->nnz[i];j++)if(p->U->index[i][j]!=q->U->index[i][j] || !isfinite((double)p->U->value[i][j]) || !isfinite((double)q->U->value[i][j]) || fabs(p->U->value[i][j]-q->U->value[i][j])>1e-12)goto done;
    }
    if(lis_vector_axpy(-1,x,y)||lis_vector_nrm2(y,&delta))goto done;
    if(!isfinite(delta)||delta/tnorm>1e-10)goto done;
    printf("EQUIVALENCE residual=%.17g known_error=%.17g fresh_delta=%.17g\n",(double)(r/bnorm),(double)(errnorm/tnorm),(double)(delta/tnorm));
    failed=0;
done:
    if(fresh)lis_solver_destroy(fresh);
    if(y)lis_vector_destroy(y);
    return failed;
}
/* A 6 x 8 grid has cycles and nontrivial elimination fill even after
 * minimum-degree ordering. Opposite diagonal signs give an indefinite
 * symmetric matrix when symmetric_values is true; the other case uses
 * nonsymmetric values on the same symmetric sparsity graph. Strict
 * diagonal dominance keeps the no-pivot reference well conditioned. */
static int test_grid(int symmetric_values)
{
    const LIS_INT n=48, width=8;
    LIS_MATRIX A=NULL;
    LIS_SOLVER s=NULL;
    LIS_VECTOR b=NULL,x=NULL,t=NULL,r=NULL;
    LIS_INT i,j,k,fill=0;
    int before_create,before_refresh,failed=1;
    OK(lis_matrix_create(LIS_COMM_WORLD,&A));
    OK(lis_matrix_set_size(A,0,n));
    for(i=0;i<n;i++)for(j=0;j<n;j++) {
        if(i==j) {
            OK(lis_matrix_set_value(LIS_INS_VALUE,i,j,
                (symmetric_values && i%2 ? -1.0 : 1.0)*(6.0+0.03*i),A));
        } else if((i/width==j/width && (i-j==1 || j-i==1)) ||
                  i-j==width || j-i==width) {
            double v=symmetric_values ? -0.25-0.01*((i+j)%7) :
                                      -0.15-0.01*((2*i+j)%11);
            OK(lis_matrix_set_value(LIS_INS_VALUE,i,j,v,A));
        }
    }
    OK(lis_matrix_set_type(A,LIS_MATRIX_CSR));
    OK(lis_matrix_assemble(A));
    OK(setup(&s));
    OK(lis_vector_duplicate(A,&b)); OK(lis_vector_duplicate(A,&x));
    OK(lis_vector_duplicate(A,&t)); OK(lis_vector_duplicate(A,&r));
    for(i=0;i<n;i++)t->value[i]=0.5+0.1*(i%13);
    CHECK(!solve_compare(A,s,b,x,t,r));
    for(i=0;i<n;i++)
        fill+=s->sequence_precon->L->nnz[i]+s->sequence_precon->U->nnz[i];
    CHECK(fill>A->nnz-n);
    before_create=creates; before_refresh=refreshes;
    for(i=0;i<n;i++)for(k=A->ptr[i];k<A->ptr[i+1];k++) {
        j=A->index[k];
        A->value[k]*=1.0+0.02*((i+j)%9);
    }
    OK(lis_matrix_mark_values_changed(A));
    CHECK(!solve_compare(A,s,b,x,t,r));
    CHECK(creates-before_create==1 && refreshes-before_refresh==1);
    printf("PASS grid %s n=%d factor_offdiag=%d input_offdiag=%d\\n",
           symmetric_values?"indefinite":"nonsymmetric",
           (int)n,(int)fill,(int)(A->nnz-n));
    failed=0;
done:
    if(s)lis_solver_destroy(s);
    if(b)lis_vector_destroy(b);
    if(x)lis_vector_destroy(x);
    if(t)lis_vector_destroy(t);
    if(r)lis_vector_destroy(r);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,B=NULL;
    LIS_SOLVER s=NULL;
    LIS_VECTOR b=NULL,x=NULL,t=NULL,r=NULL;
    LIS_INT i,status,rc;
    int before_create,before_refresh,before_errors,failed=1;
    OK(lis_initialize(&argc,&argv));
    OK(make_matrix(&A,0)); OK(make_matrix(&B,1));
    OK(setup(&s));
    CHECK(lis_solver_set_option("-ilutp_symbolic 2",s)==LIS_ERR_ILL_ARG);
    CHECK(s->ilutp_symbolic==1);
    OK(lis_vector_duplicate(A,&b)); OK(lis_vector_duplicate(A,&x));
    OK(lis_vector_duplicate(A,&t)); OK(lis_vector_duplicate(A,&r));
    for(i=0;i<N;i++)t->value[i]=0.3+i;
    CHECK(!solve_compare(A,s,b,x,t,r));
    before_create=creates; before_refresh=refreshes;
    for(i=0;i<A->nnz;i++)A->value[i]*=1.2;
    OK(lis_matrix_mark_values_changed(A));
    CHECK(!solve_compare(A,s,b,x,t,r));
    /* Only the fresh reference solver creates factors. */
    CHECK(creates-before_create==1 && refreshes-before_refresh==1);
    puts("PASS sequence numeric refresh");
    before_create=creates; before_refresh=refreshes;
    CHECK(!solve_compare(B,s,b,x,t,r));
    CHECK(creates-before_create==2 && refreshes==before_refresh);
    puts("PASS sequence changed matrix pattern rebuild");
    /* Force a numeric failure on the same CSR pattern, then restore
     * a valid matrix and verify no partially refreshed factors survive. */
    before_errors=refresh_errors; before_create=creates;
    for(i=0;i<B->nnz;i++)B->value[i]=0;
    OK(lis_matrix_mark_values_changed(B));
    OK(lis_vector_set_all(1,b)); OK(lis_vector_set_all(0,x));
    rc=lis_solve_sequence(B,b,x,s);
    OK(lis_solver_get_status(s,&status));
    CHECK(rc!=LIS_SUCCESS || status!=LIS_SUCCESS);
    CHECK(refresh_errors>before_errors && creates>before_create);
    CHECK(s->sequence_precon==NULL && s->precon==NULL);
    for(i=0;i<N;i++) {
        LIS_INT k;
        for(k=B->ptr[i];k<B->ptr[i+1];k++)
            B->value[k]=i==B->index[k]?8.0+i:0.1*(i+1)+0.2*(B->index[k]+1);
    }
    OK(lis_matrix_mark_values_changed(B));
    CHECK(!solve_compare(B,s,b,x,t,r));
    puts("PASS sequence failed refresh cleanup and recovery");
    /* An option change invalidates the sequence snapshot. */
    OK(lis_solver_set_option("-ilutp_reuse 0",s));
    before_create=creates; before_refresh=refreshes;
    OK(lis_vector_set_all(0,x));
    OK(lis_solve_sequence(B,b,x,s));
    CHECK(creates>before_create && refreshes==before_refresh && s->sequence_precon->ilutp_reuse==NULL);
    CHECK(s->sequence_precon->ilutp_symbolic_factor && s->sequence_precon->ilutp_symbolic_valid);
    puts("PASS sequence option snapshot rebuild");
    CHECK(!test_grid(0));
    CHECK(!test_grid(1));
    failed=0;
done:
    if(s)lis_solver_destroy(s);
    if(b)lis_vector_destroy(b);
    if(x)lis_vector_destroy(x);
    if(t)lis_vector_destroy(t);
    if(r)lis_vector_destroy(r);
    if(A)lis_matrix_destroy(A);
    if(B)lis_matrix_destroy(B);
    lis_finalize();
    return failed;
}
