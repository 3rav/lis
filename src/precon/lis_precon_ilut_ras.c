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
/* Explicitly unsupported: never drop imaginary parts in contact-Schur math. */
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
#else
#define LIS_ILUT_RAS_SCHUR_LIMIT 256

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
    LIS_INT *gamma,*gamma_of,*pivots;
    LIS_INT ngamma;
    LIS_SCALAR *diag_inv,*schur_lu,*schur_rhs,*schur_delta;
    LIS_INT schur_sparse;
    LIS_MATRIX schur_A;
    LIS_SOLVER schur_solver;
    LIS_PRECON schur_pc;
    LIS_VECTOR schur_b,schur_x;
    LIS_VECTOR residual;
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
    lis_free(s->gamma);
    lis_free(s->gamma_of);
    lis_free(s->pivots);
    lis_free(s->diag_inv);
    lis_free(s->schur_lu);
    lis_free(s->schur_rhs);
    lis_free(s->schur_delta);
    if(s->schur_pc)
    {
        if(s->schur_solver) s->schur_solver->precon=NULL;
        lis_precon_destroy(s->schur_pc);
    }
    if(s->schur_solver) lis_solver_destroy(s->schur_solver);
    if(s->schur_b) lis_vector_destroy(s->schur_b);
    if(s->schur_x) lis_vector_destroy(s->schur_x);
    if(s->schur_A) lis_matrix_destroy(s->schur_A);
    if(s->residual) lis_vector_destroy(s->residual);
    lis_free(s);
}

/* Read-only internal test/benchmark diagnostic for the separator size. */
LIS_INT
lis_precon_ilut_ras_interface_size(LIS_PRECON pc,LIS_INT *ngamma)
{
    if(!pc || !pc->ilut_ras || !ngamma) return LIS_ERR_ILL_ARG;
    *ngamma=pc->ilut_ras->ngamma;
    return LIS_SUCCESS;
}


/* Count stored factor entries, including overlapping local ILUT copies,
 * the dense Schur LU or sparse Schur ILUT, and sparse Schur CSR.
 *
 * This is a transparent *entry-count proxy*, NOT peak RSS: it excludes
 * indices, pivots, work vectors, CSR overhead and temporary allocations.
 */
LIS_INT
lis_precon_ilut_ras_factor_entries(LIS_PRECON pc,LIS_INT *entries)
{
    struct LIS_ILUT_RAS_STRUCT *s;
    LIS_INT d,nl=0,nu=0,err;
    LIS_REAL smallest=0.0,largest=0.0;
    size_t total=0,limit=(sizeof(LIS_INT)==4 ?
                          (size_t)INT32_MAX : (size_t)INT64_MAX);
    if(!pc || !pc->ilut_ras || !entries) return LIS_ERR_ILL_ARG;
    s=pc->ilut_ras;
    for(d=0;d<s->ndomains;d++)
    {
        if(!s->parts || !s->parts[d].pc) return LIS_ERR_ILL_ARG;
        err=lis_precon_ilut_factor_stats(s->parts[d].pc,
                            &nl,&nu,&smallest,&largest);
        if(err || nl<0 || nu<0) return err ? err : LIS_ERR_ILL_ARG;
        if((size_t)nl>limit-total) return LIS_ERR_ILL_ARG;
        total+=(size_t)nl;
        if((size_t)nu>limit-total) return LIS_ERR_ILL_ARG;
        total+=(size_t)nu;
    }
    if(s->ngamma>0 && s->schur_sparse)
    {
        if(!s->schur_A || !s->schur_pc ||
           s->schur_A->nnz<0) return LIS_ERR_ILL_ARG;
        if((size_t)s->schur_A->nnz>limit-total)
            return LIS_ERR_ILL_ARG;
        total+=(size_t)s->schur_A->nnz;
        err=lis_precon_ilut_factor_stats(s->schur_pc,
                            &nl,&nu,&smallest,&largest);
        if(err || nl<0 || nu<0) return err ? err : LIS_ERR_ILL_ARG;
        if((size_t)nl>limit-total) return LIS_ERR_ILL_ARG;
        total+=(size_t)nl;
        if((size_t)nu>limit-total) return LIS_ERR_ILL_ARG;
        total+=(size_t)nu;
    }
    else if(s->ngamma>0)
    {
        size_t ng=(size_t)s->ngamma;
        if(ng>limit/ng || ng*ng>limit-total) return LIS_ERR_ILL_ARG;
        total+=ng*ng;
    }
    *entries=(LIS_INT)total;
    return LIS_SUCCESS;
}

/*
 * Existing OpenMP ILUT factorization splits rows according to the active
 * thread count. An ILUT factor built with T threads MUST NOT subsequently
 * be applied with a different number of threads: its block dependencies
 * differ. RAS local factors and the sparse Schur ILUT deliberately use
 * a single thread both at setup and at application.
 *
 * This changes the OpenMP ICV of the calling task only, not a global
 * environment variable. All RAS setup is currently serial.
 */
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
lis_ilut_ras_apply_serial_schur(struct LIS_ILUT_RAS_STRUCT *s)
{
    LIS_INT err;
#ifdef _OPENMP
    int previous_threads=omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    err=lis_psolve(s->schur_solver,s->schur_b,s->schur_x);
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

/*
 * Sparse CSR path for interfaces larger than the historical dense limit.
 * The row accumulator is O(ngamma) workspace, and the assembled operator
 * has a hard fill budget of 4 * nnz(A). This is still a diagonal-approximate
 * Schur operator; sparse storage alone does not fix indefinite KKT pivots.
 */
static int
lis_ilut_ras_int_cmp(const void *a,const void *b)
{
    LIS_INT x=*(const LIS_INT *)a,y=*(const LIS_INT *)b;
    return (x>y)-(x<y);
}

static LIS_INT
lis_ilut_ras_make_sparse_schur(LIS_SOLVER outer,
                              struct LIS_ILUT_RAS_STRUCT *s)
{
    LIS_MATRIX A=s->A;
    LIS_INT ng=s->ngamma,g,i,k,j,q,t,col,count,h,err=LIS_SUCCESS;
    LIS_INT *marks=NULL,*touched=NULL;
    LIS_SCALAR *work=NULL;
    LIS_REAL row_max,mag,drop;
    size_t max_entries,nnz=0;

    if(A->nnz<0 || (size_t)A->nnz>SIZE_MAX/4 ||
       (size_t)ng>SIZE_MAX/sizeof(LIS_SCALAR) ||
       (size_t)ng>SIZE_MAX/sizeof(LIS_INT))
        return LIS_ERR_ILL_ARG;
    max_entries=(size_t)A->nnz*4;
    drop=fabs(outer->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN]);
    if(!isfinite((double)drop) || drop>=1.0) return LIS_ERR_ILL_ARG;
    if(drop<1.0e-4) drop=1.0e-4;

    marks=(LIS_INT *)lis_malloc((size_t)ng*sizeof(LIS_INT),
                               "RAS sparse Schur markers");
    touched=(LIS_INT *)lis_malloc((size_t)ng*sizeof(LIS_INT),
                                 "RAS sparse Schur touched");
    work=(LIS_SCALAR *)lis_malloc((size_t)ng*sizeof(LIS_SCALAR),
                                 "RAS sparse Schur accumulator");
    if(!marks || !touched || !work)
    {
        err=LIS_OUT_OF_MEMORY;
        goto done;
    }
    for(g=0;g<ng;g++) marks[g]=-1;

    err=lis_matrix_create(LIS_COMM_WORLD,&s->schur_A);
    if(err) goto done;
    err=lis_matrix_set_size(s->schur_A,0,ng);
    if(err) goto done;
    err=lis_matrix_set_type(s->schur_A,LIS_MATRIX_CSR);
    if(err) goto done;

    for(g=0;g<ng;g++)
    {
        i=s->gamma[g];
        count=0;
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            j=A->index[k];
            if(j<0 || j>=s->n) { err=LIS_ERR_ILL_ARG; goto done; }
            col=s->gamma_of[j];
            if(col>=0)
            {
                if(marks[col]!=g)
                {
                    marks[col]=g;
                    touched[count++]=col;
                    work[col]=0.0;
                }
                work[col]+=A->value[k];
                if(!isfinite((double)work[col]))
                { err=LIS_ERR_ILL_ARG; goto done; }
            }
            else
            {
                LIS_SCALAR coeff=-A->value[k]*s->diag_inv[j];
                if(!isfinite((double)coeff))
                { err=LIS_ERR_ILL_ARG; goto done; }
                for(q=A->ptr[j];q<A->ptr[j+1];q++)
                {
                    t=A->index[q];
                    if(t<0 || t>=s->n) { err=LIS_ERR_ILL_ARG; goto done; }
                    col=s->gamma_of[t];
                    if(col<0) continue;
                    if(marks[col]!=g)
                    {
                        marks[col]=g;
                        touched[count++]=col;
                        work[col]=0.0;
                    }
                    work[col]+=coeff*A->value[q];
                    if(!isfinite((double)work[col]))
                    { err=LIS_ERR_ILL_ARG; goto done; }
                }
            }
        }
        row_max=0.0;
        for(h=0;h<count;h++)
        {
            mag=fabs(work[touched[h]]);
            if(mag>row_max) row_max=mag;
        }
        qsort(touched,(size_t)count,sizeof(LIS_INT),lis_ilut_ras_int_cmp);
        for(h=0;h<count;h++)
        {
            col=touched[h];
            mag=fabs(work[col]);
            /* Preserve the diagonal even if below the rowwise drop bound. */
            if(col!=g && (mag==0.0 || mag<drop*row_max)) continue;
            if(nnz>=max_entries)
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "RAS sparse Schur fill exceeds 4*nnz(A) budget");
                err=LIS_ERR_ILL_ARG;
                goto done;
            }
            err=lis_matrix_set_value(LIS_INS_VALUE,g,col,
                                     work[col],s->schur_A);
            if(err) goto done;
            nnz++;
        }
    }
    err=lis_matrix_assemble(s->schur_A);
    if(err) goto done;

    err=lis_solver_create(&s->schur_solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i gmres -p ilut -scale none -print none",s->schur_solver);
    if(err) goto done;
    s->schur_solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
    s->schur_solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
    s->schur_solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
    s->schur_solver->params[LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN] =
        outer->params[LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];
    /* Reuse serial ordering/scaling path already validated in PR #224. */
    s->schur_solver->ilut_row_scale=1;
    s->schur_solver->ilut_rcm=outer->ilut_rcm;
    err=lis_solver_set_precon_half(s->schur_solver,
                                   outer->precon_half_mode);
    if(err) goto done;
    err=lis_solver_set_matrix(s->schur_A,s->schur_solver);
    if(err) goto done;
    err=lis_ilut_ras_create_serial_factor(s->schur_solver,&s->schur_pc);
    if(err) { s->schur_pc=NULL; goto done; }
    s->schur_solver->precon=s->schur_pc;
    err=lis_vector_duplicate(s->schur_A,&s->schur_b);
    if(err) goto done;
    err=lis_vector_duplicate(s->schur_A,&s->schur_x);
    if(err) goto done;
    s->schur_sparse=1;
done:
    lis_free(marks);
    lis_free(touched);
    lis_free(work);
    return err;
}

/* Build S~=A_GG - A_GI diag(A_II)^-1 A_IG, a documented
 * approximate block Schur complement, never the exact Schur matrix.
 */
static LIS_INT
lis_ilut_ras_make_schur(LIS_SOLVER outer,struct LIS_ILUT_RAS_STRUCT *s)
{
    LIS_MATRIX A=s->A;
    LIS_INT i,j,k,g,q,t,p,best,ng=0,promoted=0,err;
    LIS_REAL val,peak,schur_norm=0.0;
    LIS_SCALAR tmp;
    size_t entries;

    s->gamma_of=(LIS_INT *)lis_malloc((size_t)s->n*sizeof(LIS_INT),
                                      "RAS gamma map");
    s->diag_inv=(LIS_SCALAR *)lis_malloc((size_t)s->n*sizeof(LIS_SCALAR),
                                       "RAS interior diag");
    if(!s->gamma_of || !s->diag_inv) return LIS_OUT_OF_MEMORY;
    for(i=0;i<s->n;i++) { s->gamma_of[i]=-1; s->diag_inv[i]=0.0; }
    for(i=0;i<s->n;i++)
    for(k=A->ptr[i];k<A->ptr[i+1];k++)
    {
        j=A->index[k];
        if(s->owner[i]!=s->owner[j])
        {
            s->gamma_of[i]=0;
            s->gamma_of[j]=0;
        }
    }
    /*
     * Weak KKT interior diagonals cannot be inverted safely. Promote the
     * affected rows into Gamma before allocating or assembling the Schur
     * complement. Never silently replace their inverse by zero.
     *
     * This is deliberately bounded: promotion can enlarge Gamma and
     * worsen separator quality, and it does not repair a singular Schur.
     */
    for(i=0;i<s->n;i++) if(s->gamma_of[i]<0)
    {
        LIS_REAL row_max=0.0,diag_abs=0.0;
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            LIS_REAL mag=fabs(A->value[k]);
            if(!isfinite((double)mag))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUT RAS Schur input contains nonfinite entry");
                return LIS_ERR_ILL_ARG;
            }
            if(mag>row_max) row_max=mag;
            if(A->index[k]==i) diag_abs=mag;
        }
        if(!(row_max>0.0) || !(diag_abs>1.0e-10*row_max))
        {
            if(promoted>=s->n/5)
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUT RAS Schur weak-diagonal promotions exceed 20 percent");
                return LIS_ERR_ILL_ARG;
            }
            s->gamma_of[i]=0;
            promoted++;
        }
    }
    for(i=0;i<s->n;i++) if(s->gamma_of[i]>=0) ++ng;
    if(ng==0) return LIS_SUCCESS;
    s->ngamma=ng;
    s->gamma=(LIS_INT *)lis_malloc((size_t)ng*sizeof(LIS_INT),"RAS gamma rows");
    s->pivots=(LIS_INT *)lis_malloc((size_t)ng*sizeof(LIS_INT),"RAS pivot perm");
    s->schur_rhs=(LIS_SCALAR *)lis_malloc((size_t)ng*sizeof(LIS_SCALAR),
                                          "RAS schur rhs");
    s->schur_delta=(LIS_SCALAR *)lis_malloc((size_t)ng*sizeof(LIS_SCALAR),
                                            "RAS schur solution");
    if(ng<=LIS_ILUT_RAS_SCHUR_LIMIT)
    {
        entries=(size_t)ng*(size_t)ng;
        s->schur_lu=(LIS_SCALAR *)lis_malloc(entries*sizeof(LIS_SCALAR),
                                            "RAS Schur LU");
    }
    if(!s->gamma || !s->pivots || !s->schur_rhs ||
       !s->schur_delta ||
       (ng<=LIS_ILUT_RAS_SCHUR_LIMIT && !s->schur_lu))
        return LIS_OUT_OF_MEMORY;
    if(s->schur_lu)
        memset(s->schur_lu,0,entries*sizeof(LIS_SCALAR));

    ng=0;
    for(i=0;i<s->n;i++)
    {
        if(s->gamma_of[i]>=0)
        {
            s->gamma_of[i]=ng;
            s->gamma[ng]=i;
            s->pivots[ng]=ng;
            ng++;
        }
        else
        {
            /*
             * diag(A_II)^-1 is meaningful only when every interior row
             * has a usable diagonal relative to its own coefficient
             * scale. Quietly dropping a zero/tiny KKT multiplier pivot
             * would create a misleading approximate Schur operator.
             */
            LIS_REAL row_max=0.0,diag_abs=0.0;
            LIS_SCALAR diag=0.0;
            for(k=A->ptr[i];k<A->ptr[i+1];k++)
            {
                LIS_REAL mag=fabs(A->value[k]);
                if(!isfinite((double)mag))
                {
                    LIS_SETERR(LIS_ERR_ILL_ARG,
                        "ILUT RAS Schur interior contains nonfinite entry");
                    return LIS_ERR_ILL_ARG;
                }
                if(mag>row_max) row_max=mag;
                if(A->index[k]==i) { diag=A->value[k]; diag_abs=mag; }
            }
            if(!(row_max>0.0) || !(diag_abs>1.0e-10*row_max))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUT RAS Schur interior has missing/weak diagonal");
                return LIS_ERR_ILL_ARG;
            }
            s->diag_inv[i]=1.0/diag;
            if(!isfinite((double)s->diag_inv[i]))
            {
                LIS_SETERR(LIS_ERR_ILL_ARG,
                    "ILUT RAS Schur interior reciprocal overflows");
                return LIS_ERR_ILL_ARG;
            }
        }
    }

    if(s->ngamma>LIS_ILUT_RAS_SCHUR_LIMIT)
    {
        err=lis_ilut_ras_make_sparse_schur(outer,s);
        if(err) return err;
        return lis_vector_duplicate(A,&s->residual);
    }

    for(g=0;g<s->ngamma;g++)
    {
        i=s->gamma[g];
        for(k=A->ptr[i];k<A->ptr[i+1];k++)
        {
            j=A->index[k];
            if(s->gamma_of[j]>=0)
                s->schur_lu[(size_t)g*s->ngamma+s->gamma_of[j]]+=A->value[k];
            else if(s->diag_inv[j]!=0.0)
                for(q=A->ptr[j];q<A->ptr[j+1];q++)
                {
                    t=A->index[q];
                    if(s->gamma_of[t]>=0)
                        s->schur_lu[(size_t)g*s->ngamma+s->gamma_of[t]]
                          -=A->value[k]*s->diag_inv[j]*A->value[q];
                }
        }
    }
    /*
     * Reject nonfinite or dangerously small pivots relative to the
     * overall approximate Schur scale, not a hard-coded absolute unit.
     * This is a safety gate, NOT a convergence/conditioning guarantee.
     */
    for(i=0;i<s->ngamma*s->ngamma;i++)
    {
        val=fabs(s->schur_lu[i]);
        if(!isfinite((double)val))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "ILUT RAS approximate Schur has nonfinite coefficients");
            return LIS_ERR_ILL_ARG;
        }
        if(val>schur_norm) schur_norm=val;
    }
    if(!(schur_norm>0.0))
    {
        LIS_SETERR(LIS_ERR_ILL_ARG,
            "ILUT RAS approximate Schur is zero");
        return LIS_ERR_ILL_ARG;
    }

    /* Partial-pivot LU: preserve row-permutation for each RHS. */
    for(k=0;k<s->ngamma;k++)
    {
        best=k;peak=fabs(s->schur_lu[(size_t)k*s->ngamma+k]);
        for(i=k+1;i<s->ngamma;i++)
        {
            val=fabs(s->schur_lu[(size_t)i*s->ngamma+k]);
            if(val>peak) { best=i;peak=val; }
        }
        if(!isfinite((double)peak) || peak<=1.0e-12*schur_norm)
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "ILUT RAS approximate Schur complement has singular pivot");
            return LIS_ERR_ILL_ARG;
        }
        if(best!=k)
        {
            for(j=0;j<s->ngamma;j++)
            {
                tmp=s->schur_lu[(size_t)k*s->ngamma+j];
                s->schur_lu[(size_t)k*s->ngamma+j]=
                    s->schur_lu[(size_t)best*s->ngamma+j];
                s->schur_lu[(size_t)best*s->ngamma+j]=tmp;
            }
            p=s->pivots[k]; s->pivots[k]=s->pivots[best]; s->pivots[best]=p;
        }
        for(i=k+1;i<s->ngamma;i++)
        {
            s->schur_lu[(size_t)i*s->ngamma+k] /=
                s->schur_lu[(size_t)k*s->ngamma+k];
            for(j=k+1;j<s->ngamma;j++)
                s->schur_lu[(size_t)i*s->ngamma+j] -=
                    s->schur_lu[(size_t)i*s->ngamma+k]*
                    s->schur_lu[(size_t)k*s->ngamma+j];
        }
    }

    return lis_vector_duplicate(A,&s->residual);
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
    if(solver->ilut_schur)
    {
        err=lis_ilut_ras_make_schur(solver,s);
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
lis_psolve_ilut_ras(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
    struct LIS_ILUT_RAS_STRUCT *s=solver->precon->ilut_ras;
    LIS_INT d,i,k,gi,err,ng;
    LIS_INT domain_err[64];
    LIS_SCALAR val;
    if(!s || !B || !X || B->n!=s->n || X->n!=s->n)
        return LIS_ERR_ILL_ARG;
    err=lis_vector_set_all((LIS_SCALAR)0.0,X);
    if(err) return err;
    /*
     * Local solves have disjoint solver/work vectors, share only read-only
     * input, and inject exclusively into their nonoverlapping owner cores.
     * Report errors after the join, never return from the OpenMP region.
     * Setup remains serial because LIS solver/preconditioner construction
     * mutates shared dispatch and error-reporting state.
     */
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1) if(s->parallel)
#endif
    for(d=0;d<s->ndomains;d++)
    {
        struct LIS_ILUT_RAS_PART *p=&s->parts[d];
        LIS_INT local_err,local_i,local_g;
        domain_err[d]=LIS_SUCCESS;
        for(local_i=0;local_i<p->nlocal;local_i++)
            p->b->value[local_i]=B->value[p->global[local_i]];
        local_err=lis_ilut_ras_apply_local(p);
        if(local_err)
        {
            domain_err[d]=local_err;
            continue;
        }
        for(local_i=0;local_i<p->nlocal;local_i++)
        {
            local_g=p->global[local_i];
            if(s->owner[local_g]==d)
                X->value[local_g]=p->x->value[local_i];
        }
    }
    for(d=0;d<s->ndomains;d++)
        if(domain_err[d]!=LIS_SUCCESS) return domain_err[d];
    if(s->ngamma==0) return LIS_SUCCESS;
    ng=s->ngamma;
    err=lis_matvec(s->A,X,s->residual);
    if(err) return err;
    for(i=0;i<ng;i++)
        s->schur_rhs[i] =
            B->value[s->gamma[s->pivots[i]]] -
            s->residual->value[s->gamma[s->pivots[i]]];
    if(s->schur_sparse)
    {
        /* One ILUT application to the sparse approximate interface RHS,
         * not an iterative or exact Schur solve. */
        for(i=0;i<ng;i++) s->schur_b->value[i]=s->schur_rhs[i];
        err=lis_ilut_ras_apply_serial_schur(s);
        if(err) return err;
        for(i=0;i<ng;i++)
        {
            s->schur_delta[i]=s->schur_x->value[i];
            if(!isfinite((double)s->schur_delta[i]))
                return LIS_ERR_ILL_ARG;
        }
    }
    else
    {
        /* Small interfaces retain the existing pivoted dense LU path. */
        for(i=0;i<ng;i++)
        {
            val=s->schur_rhs[i];
            for(k=0;k<i;k++) val-=s->schur_lu[(size_t)i*ng+k]*
                                  s->schur_rhs[k];
            s->schur_rhs[i]=val;
        }
        for(i=ng-1;i>=0;i--)
        {
            val=s->schur_rhs[i];
            for(k=i+1;k<ng;k++) val-=s->schur_lu[(size_t)i*ng+k]*
                                  s->schur_delta[k];
            s->schur_delta[i]=val/s->schur_lu[(size_t)i*ng+i];
            if(!isfinite((double)s->schur_delta[i]))
                return LIS_ERR_ILL_ARG;
        }
    }
    for(i=0;i<ng;i++) X->value[s->gamma[i]]+=s->schur_delta[i];
    /* Diagonal-approximate interior back substitution. */
    for(i=0;i<s->n;i++)
    {
        if(s->gamma_of[i]>=0 || s->diag_inv[i]==0.0) continue;
        val=0.0;
        for(k=s->A->ptr[i];k<s->A->ptr[i+1];k++)
        {
            gi=s->gamma_of[s->A->index[k]];
            if(gi>=0) val+=s->A->value[k]*s->schur_delta[gi];
        }
        X->value[i]-=s->diag_inv[i]*val;
    }
    return LIS_SUCCESS;
}

#endif /* real-valued ILUT RAS implementation */
