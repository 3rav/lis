/* Standalone, opt-in adapter for CalculiX's original ITG/icol/irow/ad/au
 * format. Keep CCX's normal SPOOLES build and selection unchanged.
 * Critical: lis.h defines the layout of public LIS structs according to
 * configured CFLAGS; silently compiling without its 64-bit macro can lead
 * to misaligned pointers/segfaults even if liblis itself is correct.
 */
#ifndef HAVE_CONFIG_H
#error "Compile the CCX USER adapter against liblis's generated lis_config.h with matching CFLAGS"
#endif
#include "ccx_lis_user_adapter.h"
#include "lis_matrix_split_view.h"
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct CCX_LIS_USER_CONTEXT {
    ITG n, nnz, nzs3;
    int asymmetric;
    LIS_MATRIX A;
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view;
    LIS_MATRIX_AUTO_SESSION session;
    LIS_SOLVER solver;  /* borrowed from session */
    LIS_VECTOR b,x,r;
};

static LIS_INT checked_arrays(ITG nnz, const double *ad,
    const double *au, size_t au_capacity, ITG nzs3, int asymmetric)
{
    uint64_t nn=(uint64_t)nnz, offset;
    if(nnz<0 || !ad || (nnz>0 && !au) ||
       (asymmetric!=0 && asymmetric!=1))return LIS_ERR_ILL_ARG;
    if(nn>au_capacity)return LIS_ERR_ILL_ARG;
    if(asymmetric) {
        if(nzs3<0 || (uint64_t)nzs3<nn)return LIS_ERR_ILL_ARG;
        offset=(uint64_t)nzs3;
        if(offset>au_capacity || nn>au_capacity-offset)
            return LIS_ERR_ILL_ARG;
    }
    return LIS_SUCCESS;
}

/* Deliberately O(n+nnz): this is a checked validation adapter, not a
 * zero-overhead production interface. It rejects corrupt CCX tangents
 * before exposing borrowed arrays to a Krylov method. */
static LIS_INT checked_finite_values(ITG n,ITG nnz,
    const double *ad,const double *au,ITG nzs3,int asymmetric)
{
    ITG i;
    for(i=0;i<n;i++)if(!isfinite(ad[i]))return LIS_ERR_ILL_ARG;
    for(i=0;i<nnz;i++) {
        if(!isfinite(au[i]))return LIS_ERR_ILL_ARG;
        if(asymmetric && !isfinite(au[(size_t)nzs3+(size_t)i]))
            return LIS_ERR_ILL_ARG;
    }
    return LIS_SUCCESS;
}

static LIS_INT build_matrix(ITG n, ITG nnz,
    const ITG *icol, const ITG *irow,
    const double *ad, const double *au, size_t au_capacity,
    ITG nzs3, int asymmetric,
    LIS_MATRIX *matrix,LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW *view)
{
    int64_t *ptr=NULL,*idx=NULL;
    uint64_t count=0;
    ITG col,k,previous,used;
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology=NULL;
    LIS_INT err=LIS_ERR_ILL_ARG;
    if(!matrix||!view)return err;
    *matrix=NULL;*view=NULL;
    /* Detect mismatched CalculiX/LIS integer ABIs before constructing
     * objects; this makes the caller's LONGLONG/_LONG__LONG pairing explicit.
     * In particular, the 64-bit CCX header must not be used with 32-bit
     * LIS_INT (or vice versa). */
    if(sizeof(ITG)!=sizeof(LIS_INT))return LIS_ERR_ILL_ARG;
    if(n<=0 || nnz<0 || nnz>=INT64_MAX || !icol ||
       (nnz>0 && !irow) ||
       checked_arrays(nnz,ad,au,au_capacity,nzs3,asymmetric))
        return err;
    /* The i64 topology entrypoint itself rejects out-of-range LIS_INT,
     * but check before allocating/copying a potentially huge CCX pattern. */
    if(sizeof(LIS_INT)==sizeof(int32_t)) {
        if((uint64_t)n>INT32_MAX || (uint64_t)nnz>INT32_MAX-1)
            return LIS_ERR_ILL_ARG;
    } else if(sizeof(LIS_INT)!=sizeof(int64_t))
        return LIS_ERR_NOT_IMPLEMENTED;
    err=checked_finite_values(n,nnz,ad,au,nzs3,asymmetric);
    if(err)return err;
    err=LIS_ERR_ILL_ARG;
    if((uint64_t)n>=SIZE_MAX/sizeof(int64_t) ||
       (uint64_t)nnz>SIZE_MAX/sizeof(int64_t))
        return LIS_ERR_ILL_ARG;
    ptr=(int64_t *)malloc(((size_t)n+1)*sizeof(*ptr));
    idx=(int64_t *)malloc(((size_t)nnz? (size_t)nnz:1)*sizeof(*idx));
    if(!ptr||!idx){err=LIS_ERR_OUT_OF_MEMORY;goto done;}
    ptr[0]=1;
    for(col=0;col<n;col++) {
        if(icol[col]<0 || (uint64_t)icol[col]>(uint64_t)nnz-count)
            goto done;
        previous=col+1;
        for(used=0;used<icol[col];used++) {
            ITG row=irow[count+(uint64_t)used];
            if(row<=previous || row>n)goto done;
            idx[count+(uint64_t)used]=(int64_t)row;
            previous=row;
        }
        count+=(uint64_t)icol[col];
        ptr[col+1]=(int64_t)(count+1);
    }
    if(count!=(uint64_t)nnz)goto done;
    err=lis_matrix_create_split_triangular_topology_i64(
        (int64_t)n,(int64_t)nnz,ptr,idx,1,&topology);
    if(err)goto done;
    err=lis_matrix_create_split_triangular_shared_view_real64(
        LIS_COMM_WORLD,topology,ad,au,
        (asymmetric && nnz>0)?au+(size_t)nzs3:NULL,matrix,view);
done:
    if(topology)lis_matrix_split_triangular_topology_destroy(topology);
    if(err!=LIS_SUCCESS) {
        /* Do not publish partially created LIS objects to the caller.
         * In particular, rebind() returns immediately on build failure. */
        if(*matrix)lis_matrix_destroy(*matrix);
        if(*view)lis_matrix_split_triangular_shared_view_destroy(*view);
        *matrix=NULL;
        *view=NULL;
    }
    free(idx);free(ptr);
    return err;
}

static LIS_INT alloc_vectors(CCX_LIS_USER_CONTEXT *ctx,
                             LIS_VECTOR *b,LIS_VECTOR *x,LIS_VECTOR *r)
{
    LIS_INT err;
    *b=*x=*r=NULL;
    err=lis_vector_duplicate(ctx->A,b);if(err)return err;
    err=lis_vector_duplicate(ctx->A,x);if(err)return err;
    err=lis_vector_duplicate(ctx->A,r);
    return err;
}

void ccx_lis_user_close(CCX_LIS_USER_CONTEXT *ctx)
{
    if(!ctx)return;
    if(ctx->session)lis_matrix_auto_session_destroy(ctx->session);
    if(ctx->r)lis_vector_destroy(ctx->r);
    if(ctx->x)lis_vector_destroy(ctx->x);
    if(ctx->b)lis_vector_destroy(ctx->b);
    if(ctx->A)lis_matrix_destroy(ctx->A);
    if(ctx->view)lis_matrix_split_triangular_shared_view_destroy(ctx->view);
    free(ctx);
}

LIS_INT ccx_lis_user_open(
    ITG n, ITG nnz, const ITG icol[], const ITG irow[],
    const double ad[], const double au[], size_t au_capacity,
    ITG nzs3, int asymmetric, const char *options,
    CCX_LIS_USER_CONTEXT **out)
{
    CCX_LIS_USER_CONTEXT *ctx=NULL;
    LIS_INT err;
    const char *defaults="-i gmres -p none -tol 1e-11 -maxiter 200 -scale none -print none";
    char *mutable_options=NULL;
    if(!out)return LIS_ERR_ILL_ARG;
    *out=NULL;
    ctx=(CCX_LIS_USER_CONTEXT *)calloc(1,sizeof(*ctx));
    if(!ctx)return LIS_ERR_OUT_OF_MEMORY;
    ctx->n=n;ctx->nnz=nnz;ctx->nzs3=nzs3;ctx->asymmetric=asymmetric;
    err=build_matrix(n,nnz,icol,irow,ad,au,au_capacity,nzs3,
                     asymmetric,&ctx->A,&ctx->view);
    if(err)goto done;
    err=lis_matrix_auto_session_create(&ctx->session);if(err)goto done;
    err=lis_matrix_auto_session_get_solver(ctx->session,&ctx->solver);
    if(err)goto done;
    if(!options)options=defaults;
    mutable_options=(char *)malloc(strlen(options)+1);
    if(!mutable_options){err=LIS_ERR_OUT_OF_MEMORY;goto done;}
    strcpy(mutable_options,options);
    err=lis_solver_set_option(mutable_options,ctx->solver);
    if(err)goto done;
    err=alloc_vectors(ctx,&ctx->b,&ctx->x,&ctx->r);if(err)goto done;
    *out=ctx;
    ctx=NULL;
done:
    free(mutable_options);
    if(ctx)ccx_lis_user_close(ctx);
    return err;
}

LIS_INT ccx_lis_user_update(CCX_LIS_USER_CONTEXT *ctx,
    const double ad[],const double au[],size_t au_capacity,ITG nzs3)
{
    LIS_INT err;
    if(!ctx)return LIS_ERR_ILL_ARG;
    err=checked_arrays(ctx->nnz,ad,au,au_capacity,nzs3,ctx->asymmetric);
    if(err)return err;
    err=checked_finite_values(ctx->n,ctx->nnz,ad,au,nzs3,ctx->asymmetric);
    if(err)return err;
    err=lis_matrix_split_triangular_shared_view_set_values_real64(
        ctx->view,ad,au,(ctx->asymmetric && ctx->nnz>0)?au+(size_t)nzs3:NULL);
    if(!err)ctx->nzs3=nzs3;
    return err;
}

LIS_INT ccx_lis_user_rebind(CCX_LIS_USER_CONTEXT *ctx,
    ITG nnz,const ITG icol[],const ITG irow[],
    const double ad[],const double au[],size_t au_capacity,ITG nzs3)
{
    LIS_MATRIX next=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW next_view=NULL;
    LIS_VECTOR b=NULL,x=NULL,r=NULL;
    LIS_INT err;
    if(!ctx)return LIS_ERR_ILL_ARG;
    err=build_matrix(ctx->n,nnz,icol,irow,ad,au,au_capacity,nzs3,
                     ctx->asymmetric,&next,&next_view);
    if(err)return err;
    err=lis_vector_duplicate(next,&b);if(err)goto failed;
    err=lis_vector_duplicate(next,&x);if(err)goto failed;
    err=lis_vector_duplicate(next,&r);if(err)goto failed;

    /* Crucial: clear persistent solver/preconditioner state before the
       original source matrix and its externally-owned data are released. */
    err=lis_matrix_auto_session_reset(ctx->session);if(err)goto failed;
    if(ctx->r)lis_vector_destroy(ctx->r);
    if(ctx->x)lis_vector_destroy(ctx->x);
    if(ctx->b)lis_vector_destroy(ctx->b);
    if(ctx->A)lis_matrix_destroy(ctx->A);
    if(ctx->view)lis_matrix_split_triangular_shared_view_destroy(ctx->view);
    ctx->A=next;ctx->view=next_view;
    ctx->b=b;ctx->x=x;ctx->r=r;
    ctx->nnz=nnz;ctx->nzs3=nzs3;
    return LIS_SUCCESS;
failed:
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(next)lis_matrix_destroy(next);
    if(next_view)lis_matrix_split_triangular_shared_view_destroy(next_view);
    return err;
}

LIS_INT ccx_lis_user_solve(CCX_LIS_USER_CONTEXT *ctx,
    const double rhs[],double solution[],
    double *relative_residual,LIS_INT *policy,LIS_INT *action)
{
    LIS_INT err,status,i,actual_policy,actual_action;
    LIS_REAL residual,rhsnorm;
    if(!ctx || !rhs || !solution)return LIS_ERR_ILL_ARG;
    /* Validate all input first: a rejected RHS cannot partially replace
       the session's previous right-hand side. */
    for(i=0;i<ctx->n;i++) {
        if(!isfinite(rhs[i]))return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ctx->n;i++) {
        ctx->b->value[i]=(LIS_SCALAR)rhs[i];
    }
    err=lis_vector_set_all((LIS_SCALAR)0,ctx->x);if(err)return err;
    err=lis_matrix_auto_session_solve(ctx->session,ctx->A,ctx->b,
                         ctx->x,&actual_policy,&actual_action);
    if(err)return err;
    err=lis_solver_get_status(ctx->solver,&status);if(err)return err;
    if(status!=LIS_SUCCESS)return status;
    err=lis_matvec(ctx->A,ctx->x,ctx->r);if(err)return err;
    err=lis_vector_axpy((LIS_SCALAR)-1,ctx->b,ctx->r);if(err)return err;
    err=lis_vector_nrm2(ctx->r,&residual);if(err)return err;
    err=lis_vector_nrm2(ctx->b,&rhsnorm);if(err)return err;
    if(!isfinite((double)residual)||!isfinite((double)rhsnorm))
        return LIS_FAILS;
    /* Publish outputs only after all solver and residual checks pass. */
    for(i=0;i<ctx->n;i++)solution[i]=(double)ctx->x->value[i];
    if(relative_residual)
        *relative_residual=(double)residual/(1.0+(double)rhsnorm);
    if(policy)*policy=actual_policy;
    if(action)*action=actual_action;
    return LIS_SUCCESS;
}
