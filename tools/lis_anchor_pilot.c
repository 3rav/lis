/* Experimental retained-anchor policy. BSD license as for LIS. */
#include "lislib.h"
#include "lis_solver_options.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

struct lis_anchor_pilot {
    LIS_MATRIX P;
    LIS_SOLVER solver;
    LIS_PRECON factor;
    LIS_UNSIGNED_INT structure, values;
    int valid, low_memory;
};

static void init_report(lis_anchor_pilot_report *r)
{
    memset(r,0,sizeof(*r));
    r->raw_status=LIS_ERR_ILL_ARG;
    r->relative_residual=r->ratio=NAN;
}

#if defined(USE_MPI) || defined(_COMPLEX) || defined(USE_COMPLEX) || defined(_LONG__DOUBLE) || defined(USE_QUAD_PRECISION)
LIS_INT lis_anchor_pilot_create(LIS_MATRIX P,lis_anchor_pilot **out)
{ (void)P; if(out)*out=NULL; return LIS_ERR_NOT_IMPLEMENTED; }
LIS_INT lis_anchor_pilot_create_low_memory(LIS_MATRIX P,lis_anchor_pilot **out)
{ return lis_anchor_pilot_create(P,out); }
LIS_INT lis_anchor_pilot_solve(lis_anchor_pilot *ctx,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,lis_anchor_pilot_report *r)
{ (void)ctx;(void)A;(void)b;(void)x;if(r)init_report(r);return LIS_ERR_NOT_IMPLEMENTED; }
LIS_INT lis_anchor_fresh_solve(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,lis_anchor_pilot_report *r)
{ (void)A;(void)b;(void)x;if(r)init_report(r);return LIS_ERR_NOT_IMPLEMENTED; }
void lis_anchor_pilot_destroy(lis_anchor_pilot *ctx) { (void)ctx; }
#else
static int valid_matrix(LIS_MATRIX A)
{
#ifdef _OPENMP
    if(omp_get_max_threads()!=1)return 0;
#endif
    return A && A->label==LIS_LABEL_MATRIX && A->matrix_type==LIS_MATRIX_CSR &&
        A->nprocs==1 && A->n>0 && A->n==A->gn && !A->is_scaled && !A->is_splited &&
        lis_matrix_is_assembled(A) && A->ptr && A->index && A->value;
}
static int valid_vector(LIS_MATRIX A,LIS_VECTOR v)
{
    return v && v->label==LIS_LABEL_VECTOR && v->precision==LIS_PRECISION_DEFAULT &&
        v->nprocs==1 && v->value && v->n==A->n && v->gn==A->gn && v->is==A->is && v->ie==A->ie;
}
static int valid_args(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x)
{
    return valid_matrix(A) && valid_vector(A,b) && valid_vector(A,x) && b!=x && b->value!=x->value;
}
static LIS_INT make_solver(LIS_SOLVER *s)
{
    char options[]="-i fgmres -p ilut -scale none -print none -restart 60 -maxiter 100 -tol 1e-8 -ilut_drop 0.005 -ilut_rate 8 -ilu_pivot_tol 1e-9 -ilut_rcm 1 -ilut_row_scale 1 -ilutp 1 -ilutp_threshold 0.5 -ilutp_full 1 -ilutp_full_limit 60000000 -initx_zeros true";
    LIS_INT rc=lis_solver_create(s);
    if(!rc)rc=lis_solver_set_option(options,*s);
    if(!rc)rc=lis_solver_set_num_threads(*s,1);
    return rc;
}
static void unbind(LIS_SOLVER s,LIS_MATRIX P)
{
    /* The caller can release A, b and x after the synchronous operation. */
    if(!s)return;
    s->A=P;s->b=NULL;s->x=NULL;s->xx=NULL;
    (void)lis_solver_set_monitor(s,NULL,NULL);
}
static LIS_INT factor(LIS_MATRIX P,LIS_SOLVER *s,LIS_PRECON *pc,lis_anchor_pilot_report *r)
{
    LIS_INT rc=make_solver(s);
    double start;
    if(rc)return rc;
    start=lis_wtime();
    rc=lis_solver_sequence_prepare_precon(P,*s,pc);
    r->setup_seconds+=lis_wtime()-start;
    if(!rc)r->factorizations++;
    return rc;
}
static int numerical(LIS_INT status)
{
    return status==LIS_SUCCESS || status==LIS_MAXITER || status==LIS_BREAKDOWN;
}
static LIS_INT attempt(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR work,LIS_VECTOR residual,
    LIS_SOLVER s,LIS_PRECON pc,LIS_INT budget,int zero,LIS_INT *iterations,
    lis_anchor_pilot_report *r)
{
    LIS_INT rc;
    LIS_REAL nr,nb,reported;
    if((rc=lis_solver_set_maxiter(s,budget)))return rc;
    s->options[LIS_OPTIONS_INITGUESS_ZEROS]=zero?LIS_TRUE:LIS_FALSE;
    if(zero && (rc=lis_vector_set_all(0,work)))return rc;
    if((rc=lis_solve_kernel(A,b,work,s,pc)))return rc;
    if((rc=lis_solver_get_status(s,&r->raw_status)) ||
       (rc=lis_solver_get_iter(s,iterations)) ||
       (rc=lis_solver_get_residualnorm(s,&reported)))return rc;
    if(!numerical(r->raw_status))return r->raw_status;
    if((rc=lis_matvec(A,work,residual)) || (rc=lis_vector_axpy(-1,b,residual)) ||
       (rc=lis_vector_nrm2(residual,&nr)) || (rc=lis_vector_nrm2(b,&nb)))return rc;
    r->relative_residual=nb>0?nr/nb:nr;
    r->accepted=isfinite((double)reported) && isfinite((double)nr) && isfinite((double)nb) &&
        isfinite(r->relative_residual) && r->relative_residual<=1e-8;
    return LIS_SUCCESS;
}

typedef struct { double r10,r20; } pilot_history;
static void monitor(LIS_SOLVER s,LIS_INT iter,LIS_REAL res,LIS_REAL tol,LIS_INT kind,void *data)
{
    pilot_history *h=data;(void)s;(void)tol;
    if(kind!=LIS_SOLVER_RESIDUAL_PROJECTED)return;
    if(iter==10)h->r10=res;
    if(iter==20)h->r20=res;
}
static LIS_INT fresh_attempt(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR work,LIS_VECTOR residual,lis_anchor_pilot_report *r)
{
    LIS_SOLVER s=NULL;LIS_PRECON pc=NULL;
    LIS_INT rc=factor(A,&s,&pc,r);
    if(!rc)rc=attempt(A,b,work,residual,s,pc,100,1,&r->fresh_iterations,r);
    if(s)lis_solver_destroy(s);
    return rc;
}
static LIS_INT create(LIS_MATRIX P,lis_anchor_pilot **out,int low_memory)
{
    lis_anchor_pilot *ctx;lis_anchor_pilot_report r;LIS_INT rc;
    if(!out)return LIS_ERR_ILL_ARG;
    *out=NULL;
    if(!valid_matrix(P))return LIS_ERR_ILL_ARG;
    ctx=calloc(1,sizeof(*ctx));if(!ctx)return LIS_OUT_OF_MEMORY;
    ctx->P=P;init_report(&r);
    rc=factor(P,&ctx->solver,&ctx->factor,&r);
    if(!rc)rc=lis_matrix_get_generation(P,&ctx->structure,&ctx->values);
    if(rc){lis_anchor_pilot_destroy(ctx);return rc;}
    ctx->valid=1;ctx->low_memory=low_memory;*out=ctx;return LIS_SUCCESS;
}
LIS_INT lis_anchor_pilot_create(LIS_MATRIX P,lis_anchor_pilot **out)
{ return create(P,out,0); }
LIS_INT lis_anchor_pilot_create_low_memory(LIS_MATRIX P,lis_anchor_pilot **out)
{ return create(P,out,1); }
static void retire_anchor(lis_anchor_pilot *ctx)
{
    if(ctx->solver)lis_solver_destroy(ctx->solver);
    ctx->solver=NULL;ctx->factor=NULL;
}
void lis_anchor_pilot_destroy(lis_anchor_pilot *ctx)
{
    if(ctx){retire_anchor(ctx);free(ctx);}
}
static LIS_INT solve(lis_anchor_pilot *ctx,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,lis_anchor_pilot_report *r)
{
    LIS_VECTOR work=NULL,residual=NULL;
    LIS_INT rc=LIS_SUCCESS;
    LIS_UNSIGNED_INT structure,values;
    pilot_history h={NAN,NAN};
    if(!r)return LIS_ERR_ILL_ARG;
    init_report(r);
    if(!valid_args(A,b,x))return LIS_ERR_ILL_ARG;
    if(ctx){
        if(!ctx->valid || !valid_matrix(ctx->P) || A->n!=ctx->P->n)return LIS_ERR_ILL_ARG;
        rc=lis_matrix_get_generation(ctx->P,&structure,&values);
        if(rc)return rc;
        if(structure!=ctx->structure || values!=ctx->values)return LIS_ERR_ILL_ARG;
    }
    if((rc=lis_vector_duplicate(A,&work)) || (rc=lis_vector_duplicate(A,&residual)))goto done;
    if(ctx){
        /* A memory-constrained fallback retired the original factor.
         * Generations were checked above, before any lazy reconstruction. */
        if(!ctx->solver){
            rc=factor(ctx->P,&ctx->solver,&ctx->factor,r);
            if(rc){retire_anchor(ctx);ctx->valid=0;goto done;}
        }
        if((rc=lis_solver_set_monitor(ctx->solver,monitor,&h)))goto done;
        rc=attempt(A,b,work,residual,ctx->solver,ctx->factor,20,1,&r->pilot_iterations,r);
        unbind(ctx->solver,ctx->P);
        if(rc){ctx->valid=0;goto done;}
        if(!r->accepted && r->raw_status==LIS_MAXITER && isfinite(h.r10) && h.r10>0 && isfinite(h.r20) && h.r20>=0){
            r->ratio=h.r20/h.r10;
            if(r->ratio<=0.25){
                r->continued=1;
                rc=attempt(A,b,work,residual,ctx->solver,ctx->factor,80,0,&r->continuation_iterations,r);
                unbind(ctx->solver,ctx->P);
                if(rc){ctx->valid=0;goto done;}
            }
        }
        if(!r->accepted){
            r->fallback=1;
            if(ctx->low_memory)retire_anchor(ctx);
            rc=fresh_attempt(A,b,work,residual,r);
        }
    }else rc=fresh_attempt(A,b,work,residual,r);
    if(!rc && r->accepted)rc=lis_vector_copy(work,x);
    else if(!rc)rc=r->raw_status==LIS_SUCCESS?LIS_BREAKDOWN:r->raw_status;
done:
    if(ctx){unbind(ctx->solver,ctx->P);if(rc && !numerical(rc))ctx->valid=0;}
    if(work)lis_vector_destroy(work);
    if(residual)lis_vector_destroy(residual);
    if(rc)r->accepted=0;
    return rc;
}
LIS_INT lis_anchor_pilot_solve(lis_anchor_pilot *ctx,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,lis_anchor_pilot_report *r)
{
    if(!ctx){if(r)init_report(r);return LIS_ERR_ILL_ARG;}
    return solve(ctx,A,b,x,r);
}
LIS_INT lis_anchor_fresh_solve(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,lis_anchor_pilot_report *r)
{
    return solve(NULL,A,b,x,r);
}
#endif
