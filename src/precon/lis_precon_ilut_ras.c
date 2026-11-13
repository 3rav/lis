/*
 * Opt-in serial restricted additive Schwarz ILUT with diagonal-approximate
 * interface Schur correction. No MPI exchange: distributed matrices are
 * deliberately rejected, not silently treated as serial subdomains.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lislib.h"

#define __FUNC__ "lis_precon_ilut_ras"
#if defined(_COMPLEX) || defined(USE_COMPLEX)
/* Explicitly unsupported: never drop imaginary parts in restricted Schwarz math. */
LIS_INT lis_precon_create_ilut_ras(LIS_SOLVER solver,LIS_PRECON precon)
{
    (void)solver; (void)precon;
    return LIS_ERR_NOT_IMPLEMENTED;
}
LIS_INT lis_psolve_ilut_ras(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    (void)solver; (void)b; (void)x;
    return LIS_ERR_NOT_IMPLEMENTED;
}
void lis_precon_ilut_ras_destroy(struct LIS_ILUT_RAS_STRUCT *s)
{
    (void)s;
}
LIS_INT lis_precon_ilut_ras_interface_size(LIS_PRECON pc,LIS_INT *ngamma)
{
    (void)pc; (void)ngamma;
    return LIS_ERR_NOT_IMPLEMENTED;
}
LIS_INT lis_precon_ilut_ras_factor_entries(LIS_PRECON pc,LIS_INT *entries)
{
    (void)pc; (void)entries;
    return LIS_ERR_NOT_IMPLEMENTED;
}
LIS_INT lis_precon_ilut_ras_probe(LIS_PRECON pc,LIS_VECTOR B,
                         LIS_INT *worst_domain,LIS_REAL *max_ratio,
                         LIS_REAL *min_pivot)
{
    (void)pc; (void)B; (void)worst_domain;
    (void)max_ratio; (void)min_pivot;
    return LIS_ERR_NOT_IMPLEMENTED;
}
#else
struct LIS_ILUT_RAS_PART {
    LIS_INT start,end,nlocal;
    LIS_INT *global;         /* local index -> original row */
    LIS_MATRIX A;
    LIS_SOLVER solver;
    LIS_PRECON pc;
    LIS_VECTOR b,x;
};

struct LIS_ILUT_RAS_STRUCT {
    LIS_MATRIX A;            /* borrowed outer serial CSR */
    LIS_INT n,ndomains,overlap,parallel;
    struct LIS_ILUT_RAS_PART *parts;
    LIS_INT *owner,*order;
};

void lis_precon_ilut_ras_destroy(struct LIS_ILUT_RAS_STRUCT *s)
{
    LIS_INT d;
    if(!s) return;
    if(s->parts)
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        if(p->pc)
        {
            if(p->solver) p->solver->precon=NULL;
            lis_precon_destroy(p->pc);
        }
        if(p->solver) lis_solver_destroy(p->solver);
        if(p->b) lis_vector_destroy(p->b);
        if(p->x) lis_vector_destroy(p->x);
        if(p->A) lis_matrix_destroy(p->A);
        lis_free(p->global);
    }
    lis_free(s->parts);
    lis_free(s->owner);
    lis_free(s->order);
    lis_free(s);
}

/* Read-only internal test/benchmark diagnostic for the separator size. */
/* RAS-only: no separate interface correction. */
LIS_INT
lis_precon_ilut_ras_interface_size(LIS_PRECON pc,LIS_INT *ngamma)
{
    if(!pc || !pc->ilut_ras || !ngamma) return LIS_ERR_ILL_ARG;
    *ngamma=0;
    return LIS_SUCCESS;
}
LIS_INT
lis_precon_ilut_ras_factor_entries(LIS_PRECON pc,LIS_INT *entries)
{
    struct LIS_ILUT_RAS_STRUCT *s;
    LIS_INT d,nl,nu,err;
    LIS_REAL minpivot=0.0,maxpivot=0.0;
    size_t total=0,limit=sizeof(LIS_INT)==4 ?
        (size_t)INT32_MAX : (size_t)INT64_MAX;
    if(!pc || !pc->ilut_ras || !entries) return LIS_ERR_ILL_ARG;
    s=pc->ilut_ras;
    for(d=0;d<s->ndomains;d++)
    {
        if(!s->parts || !s->parts[d].pc) return LIS_ERR_ILL_ARG;
        err=lis_precon_ilut_factor_stats(s->parts[d].pc,
                           &nl,&nu,&minpivot,&maxpivot);
        if(err) return err;
        if(nl<0 || nu<0 || (size_t)nl>limit-total)
            return LIS_ERR_ILL_ARG;
        total+=(size_t)nl;
        if((size_t)nu>limit-total) return LIS_ERR_ILL_ARG;
        total+=(size_t)nu;
    }
    *entries=(LIS_INT)total;
    return LIS_SUCCESS;
}

static LIS_INT
lis_ilut_ras_create_serial_factor(LIS_SOLVER solver,LIS_PRECON *pc)
{
    LIS_INT err;
#ifdef _OPENMP
    int previous_threads=omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    err=lis_precon_create(solver,pc);
#ifdef _OPENMP
    omp_set_num_threads(previous_threads);
#endif
    return err;
}

static LIS_INT
lis_ilut_ras_make_subdomain(LIS_SOLVER outer, struct LIS_ILUT_RAS_STRUCT *s,
                            struct LIS_ILUT_RAS_PART *p)
{
    LIS_MATRIX A=s->A;
    LIS_INT n=s->n, depth,i,k,j,t,count=0,err;
    LIS_INT *selected=NULL,*local=NULL;
    LIS_INT new_start=p->start,new_end=p->end;

    selected=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS selection");
    local=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS global lookup");
    if(!selected || !local) { err=LIS_OUT_OF_MEMORY; goto done; }
    memset(selected,0,(size_t)n*sizeof(LIS_INT));
    for(i=new_start;i<new_end;i++) selected[s->order[i]]=1;
    for(depth=1;depth<=s->overlap;depth++)
    for(i=0;i<n;i++)
    {
        if(selected[i]!=depth) continue;
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            j=A->index[k];
            if(j<0 || j>=n) { err=LIS_ERR_ILL_ARG; goto done; }
            if(!selected[j]) selected[j]=depth+1;
        }
    }
    for(i=0;i<n;i++) if(selected[i]) count++;
    if(count<=0) { err=LIS_ERR_ILL_ARG; goto done; }
    p->global=(LIS_INT *)lis_malloc((size_t)count*sizeof(LIS_INT),"RAS local rows");
    if(!p->global) { err=LIS_OUT_OF_MEMORY; goto done; }
    p->nlocal=count;
    count=0;
    for(i=0;i<n;i++)
    {
        local[i]=-1;
        if(selected[i]) { local[i]=count; p->global[count++]=i; }
    }

    err=lis_matrix_create(LIS_COMM_WORLD,&p->A);
    if(err) goto done;
    err=lis_matrix_set_size(p->A,0,p->nlocal);
    if(err) goto done;
    err=lis_matrix_set_type(p->A,LIS_MATRIX_CSR);
    if(err) goto done;
    for(i=0;i<p->nlocal;i++)
    {
        t=p->global[i];
        for(k=A->ptr[t];k<A->ptr[t+1];k++)
        {
            j=local[A->index[k]];
            if(j<0) continue;
            err=lis_matrix_set_value(LIS_INS_VALUE,i,j,A->value[k],p->A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(p->A);
    if(err) goto done;

    err=lis_solver_create(&p->solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i gmres -p ilut -scale none -print none",p->solver);
    if(err) goto done;
    p->solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
    p->solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
    p->solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
    p->solver->params[LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];
    /* Reuse opt-in ILUT row equilibration from upstream #224 for each
     * independently factored RAS subdomain. */
    p->solver->ilut_row_scale=outer->ilut_row_scale;
    /* The separate opt-in -ilut_rcm policy orders each complete local
     * owner+halo matrix before ILUT. Domain partition RCM is independent
     * (-ilut_ras_rcm); both default to disabled. */
    p->solver->ilut_rcm=outer->ilut_rcm;
    err=lis_solver_set_precon_half(p->solver,outer->precon_half_mode);
    if(err) goto done;
    err=lis_solver_set_matrix(p->A,p->solver);
    if(err) goto done;
    err=lis_ilut_ras_create_serial_factor(p->solver,&p->pc);
    if(err) { p->pc=NULL; goto done; }
    p->solver->precon=p->pc;
    err=lis_vector_duplicate(p->A,&p->b);
    if(err) goto done;
    err=lis_vector_duplicate(p->A,&p->x);

done:
    lis_free(selected);
    lis_free(local);
    return err;
}


LIS_INT
lis_precon_create_ilut_ras(LIS_SOLVER solver, LIS_PRECON precon)
{
    struct LIS_ILUT_RAS_STRUCT *s=NULL;
    LIS_MATRIX A=solver->A;
    LIS_INT n=A->n,d,i,k,err;
    if(A->nprocs!=1 || A->n!=A->gn || A->is!=0 ||
       A->matrix_type!=LIS_MATRIX_CSR || A->ptr==NULL ||
       A->index==NULL || A->value==NULL ||
       solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE)
    {
        LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
            "ILUT RAS currently supports only serial double CSR");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    if(solver->ilut_ras_overlap<0 || solver->ilut_ras_overlap>2 ||
       solver->ilut_ras_domains<2 || solver->ilut_ras_domains>64 ||
       solver->ilut_ras_domains>n)
        return LIS_ERR_ILL_ARG;
    s=(struct LIS_ILUT_RAS_STRUCT *)lis_malloc(sizeof(*s),"ILUT RAS state");
    if(!s) return LIS_OUT_OF_MEMORY;
    memset(s,0,sizeof(*s));
    precon->ilut_ras=s; /* generic preconditioner destruction owns cleanup */
    s->A=A; s->n=n; s->ndomains=solver->ilut_ras_domains;
    s->overlap=solver->ilut_ras_overlap;
    s->parallel=solver->ilut_ras_parallel;
    s->owner=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS owners");
    s->order=(LIS_INT *)lis_malloc((size_t)n*sizeof(LIS_INT),"RAS rank ordering");
    s->parts=(struct LIS_ILUT_RAS_PART *)lis_malloc(
        (size_t)s->ndomains*sizeof(*s->parts),"RAS subdomains");
    /*
     * If one of these allocations fails, the generic preconditioner
     * destructor still visits any allocated parts. Initialize them first.
     */
    if(s->parts)
        memset(s->parts,0,(size_t)s->ndomains*sizeof(*s->parts));
    if(!s->owner || !s->parts || !s->order) return LIS_OUT_OF_MEMORY;
    /*
     * Partition contiguous *RCM ranks*, not raw global indices. Keep the
     * original matrix and its numbering unchanged: the ordering affects
     * only subdomain ownership/extraction and hence Gamma.
     *
     * The quotient/remainder split avoids integer multiplication overflow
     * for large 64-bit LIS_INT matrices.
     */
    if(solver->ilut_ras_rcm)
    {
        err=lis_ilut_rcm_order(A,s->order);
        if(err) return err;
    }
    else for(i=0;i<n;i++) s->order[i]=i;
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        p->start=(n/s->ndomains)*d +
                 ((n%s->ndomains)*d)/s->ndomains;
        p->end=(n/s->ndomains)*(d+1) +
               ((n%s->ndomains)*(d+1))/s->ndomains;
        for(i=p->start;i<p->end;i++)
            s->owner[s->order[i]]=d;
    }
    for(d=0;d<s->ndomains;d++)
    {
        err=lis_ilut_ras_make_subdomain(solver,s,&s->parts[d]);
        if(err) return err;
    }
    return LIS_SUCCESS;
}

/*
 * Each RAS domain owns its factor/workspace. Perform its triangular solve
 * serially INSIDE the outer OpenMP domain loop: nested LIS ILUT OpenMP
 * regions otherwise repartition a domain by the wrong thread count when
 * nested parallelism is disabled, in addition to unnecessary oversubscription.
 *
 * Packed half factors retain their existing dedicated solve path. Local
 * RAS construction does not enable ILUT permutation/row equilibration;
 * reject an unexpected local factor layout instead of silently misusing it.
 */
static LIS_INT
lis_ilut_ras_apply_local(struct LIS_ILUT_RAS_PART *part)
{
    LIS_MATRIX_ILU L,U;
    LIS_VECTOR D,from,to;
    LIS_SCALAR *x;
    LIS_INT n,i,j,col,err;
    LIS_PRECON pc;
    LIS_INT *perm;
    LIS_SCALAR *scales;
    if(!part || !part->solver || !part->pc || !part->b || !part->x)
        return LIS_ERR_ILL_ARG;
    pc=part->pc;
    perm=pc->ilut_perm;
    scales=pc->ilut_row_scale;
    n=part->nlocal;
    if((perm!=NULL)!=(scales!=NULL))
        return LIS_ERR_ILL_ARG;
    from=part->b;
    to=part->x;
    if(perm)
    {
        from=pc->ilut_work_b;
        to=pc->ilut_work_x;
        if(!from || !to || !from->value || !to->value)
            return LIS_ERR_ILL_ARG;
        /* The same row scaling / permutation as ordered ILUT's apply,
         * except the triangular solve is serial for every RAS domain. */
        for(i=0;i<n;i++)
        {
            col=perm[i];
            if(col<0 || col>=n ||
               !isfinite((double)scales[i]) || scales[i]<=0.0)
                return LIS_ERR_ILL_ARG;
            from->value[i]=part->b->value[col]*scales[i];
            if(!isfinite((double)from->value[i])) return LIS_ERR_ILL_ARG;
        }
    }
    if(pc->half)
    {
#ifdef _OPENMP
        int saved_threads=omp_get_max_threads();
        omp_set_num_threads(1);
#endif
        err=lis_half_ilu_psolve(part->solver,from,to);
#ifdef _OPENMP
        omp_set_num_threads(saved_threads);
#endif
        if(err) return err;
    }
    else
    {
        L=pc->L; U=pc->U; D=pc->D;
        if(!L || !U || !D || !D->value || !from->value ||
           !to->value)
            return LIS_ERR_ILL_ARG;
        x=to->value;
        for(i=0;i<n;i++) x[i]=from->value[i];
        for(i=0;i<n;i++)
        {
            for(j=0;j<L->nnz[i];j++)
            {
                col=L->index[i][j];
                if(col<0 || col>=n) return LIS_ERR_ILL_ARG;
                x[i]-=L->value[i][j]*x[col];
            }
            if(!isfinite((double)x[i])) return LIS_ERR_ILL_ARG;
        }
        for(i=n-1;i>=0;i--)
        {
            for(j=0;j<U->nnz[i];j++)
            {
                col=U->index[i][j];
                if(col<0 || col>=n) return LIS_ERR_ILL_ARG;
                x[i]-=U->value[i][j]*x[col];
            }
            x[i]*=D->value[i];
            if(!isfinite((double)x[i])) return LIS_ERR_ILL_ARG;
        }
    }
    if(perm)
    {
        for(i=0;i<n;i++)
        {
            if(!isfinite((double)to->value[i])) return LIS_ERR_ILL_ARG;
            part->x->value[perm[i]]=to->value[i];
        }
    }
    return LIS_SUCCESS;
}

/* Diagnostic only: probe local inverse responses on the actual RHS.
 * Reports the worst local ||M_d^-1 b_d||_2 / ||b_d||_2 and pivot.
 * Never changes the preconditioner or silently regularizes factors.
 */
LIS_INT
lis_precon_ilut_ras_probe(LIS_PRECON pc,LIS_VECTOR B,
                         LIS_INT *worst_domain,LIS_REAL *max_ratio,
                         LIS_REAL *min_pivot)
{
    struct LIS_ILUT_RAS_STRUCT *s;
    LIS_INT d,i,err,nnzl=0,nnzu=0;
    LIS_REAL bn=0.0,xn=0.0,pmin=0.0,pmax=0.0,ratio;
    if(!pc || !pc->ilut_ras || !B ||
       !worst_domain || !max_ratio || !min_pivot)
        return LIS_ERR_ILL_ARG;
    s=pc->ilut_ras;
    if(B->n!=s->n || !B->value) return LIS_ERR_ILL_ARG;
    *max_ratio=0.0;
    *min_pivot=LIS_SCALAR_MAX;
    *worst_domain=-1;
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        for(i=0;i<p->nlocal;i++)
            p->b->value[i]=B->value[p->global[i]];
        err=lis_vector_nrm2(p->b,&bn);
        if(err) return err;
        err=lis_precon_ilut_factor_stats(p->pc,
                    &nnzl,&nnzu,&pmin,&pmax);
        if(err) return err;
        if(pmin<*min_pivot) *min_pivot=pmin;
        if(!(bn>0.0)) continue;
        err=lis_ilut_ras_apply_local(p);
        if(err) return err;
        err=lis_vector_nrm2(p->x,&xn);
        if(err || !isfinite((double)xn))
            return LIS_ERR_ILL_ARG;
        ratio=xn/bn;
        if(!isfinite((double)ratio))
            return LIS_ERR_ILL_ARG;
        if(ratio>*max_ratio)
        {
            *max_ratio=ratio;
            *worst_domain=d;
        }
    }
    return LIS_SUCCESS;
}

LIS_INT
lis_psolve_ilut_ras(LIS_SOLVER solver,LIS_VECTOR B,LIS_VECTOR X)
{
    struct LIS_ILUT_RAS_STRUCT *s;
    LIS_INT d,err,domain_err[64];
    if(!solver || !solver->precon ||
       !(s=solver->precon->ilut_ras) || !B || !X ||
       B->n!=s->n || X->n!=s->n || !B->value || !X->value)
        return LIS_ERR_ILL_ARG;
    err=lis_vector_set_all((LIS_SCALAR)0.0,X);
    if(err) return err;
    /* Independent local factors and disjoint owner-core injections.
     * Each local factor is created/applied serially; only domains run
     * concurrently. The same preconditioner object is not reentrant.
     */
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) if(s->parallel)
#endif
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        LIS_INT i,g,local_err;
        domain_err[d]=LIS_SUCCESS;
        for(i=0;i<p->nlocal;i++)
            p->b->value[i]=B->value[p->global[i]];
        local_err=lis_ilut_ras_apply_local(p);
        if(local_err) { domain_err[d]=local_err; continue; }
        for(i=0;i<p->nlocal;i++)
        {
            g=p->global[i];
            if(s->owner[g]==d) X->value[g]=p->x->value[i];
        }
    }
    for(d=0;d<s->ndomains;d++)
        if(domain_err[d]!=LIS_SUCCESS) return domain_err[d];
    return LIS_SUCCESS;
}

#endif /* real-valued ILUT RAS implementation */
