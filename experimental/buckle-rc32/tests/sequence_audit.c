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
    if((err=lis_solver_set_option("-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 10000 -tol 1e-12 -maxiter 40 -restart 20 -print none",*s)))return err;
    return lis_solver_set_sequence_precon_refresh_policy(*s,LIS_SEQUENCE_PRECON_REFRESH_FROZEN_PATTERN);
}
static int solve_compare(LIS_MATRIX A,LIS_SOLVER s,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR truth,LIS_VECTOR scratch)
{
    LIS_SOLVER fresh=NULL;
    LIS_VECTOR y=NULL;
    LIS_REAL r,bnorm,errnorm,tnorm,delta;
    LIS_INT status;
    int failed=1;
    if(lis_matvec(A,truth,b)||lis_vector_set_all(0,x))goto done;
    if(lis_solve_sequence(A,b,x,s)||lis_solver_get_status(s,&status)||status!=LIS_SUCCESS)goto done;
    if(lis_matvec(A,x,scratch)||lis_vector_axpy(-1,b,scratch)||lis_vector_nrm2(scratch,&r)||lis_vector_nrm2(b,&bnorm))goto done;
    if(!(bnorm>0)||!isfinite(r)||r/bnorm>1e-11)goto done;
    if(lis_vector_copy(x,scratch)||lis_vector_axpy(-1,truth,scratch)||lis_vector_nrm2(scratch,&errnorm)||lis_vector_nrm2(truth,&tnorm))goto done;
    if(!(tnorm>0)||!isfinite(errnorm)||errnorm/tnorm>1e-10)goto done;
    if(setup(&fresh)||lis_vector_duplicate(A,&y)||lis_vector_set_all(0,y))goto done;
    if(lis_solve_sequence(A,b,y,fresh)||lis_solver_get_status(fresh,&status)||status!=LIS_SUCCESS)goto done;
    if(lis_vector_axpy(-1,x,y)||lis_vector_nrm2(y,&delta))goto done;
    if(!isfinite(delta)||delta/tnorm>1e-10)goto done;
    printf("EQUIVALENCE residual=%.17g known_error=%.17g fresh_delta=%.17g\n",(double)(r/bnorm),(double)(errnorm/tnorm),(double)(delta/tnorm));
    failed=0;
done:
    if(fresh)lis_solver_destroy(fresh);
    if(y)lis_vector_destroy(y);
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
