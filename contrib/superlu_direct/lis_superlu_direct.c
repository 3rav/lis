/*
 * Opt-in SuperLU direct preconditioner for LIS.
 *
 * This translation unit is NOT part of the default LIS build.
 * Register it with lis_superlu_direct_register() after lis_initialize(),
 * then use "-p superlu". Real double, serial CSR only. Factor once, reuse.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lislib.h"
#include "lis_superlu_direct.h"
#include <superlu/slu_ddefs.h>

typedef struct {
    int n, initialized, have_A, have_AC, have_L, have_U;
    int_t nnz;
    double *nzval, *R, *C, *rhs;
    int_t *rowind, *colptr;
    int *perm_c, *perm_r, *etree;
    char equed;
    SuperMatrix A,AC,L,U;
} LIS_SLU_FACTORS;

static void lis_slu_release(LIS_SLU_FACTORS *s)
{
    if(!s) return;
    if(s->have_L) Destroy_SuperNode_Matrix(&s->L);
    if(s->have_U) Destroy_CompCol_Matrix(&s->U);
    if(s->have_AC) Destroy_CompCol_Permuted(&s->AC);
    if(s->have_A) Destroy_SuperMatrix_Store(&s->A);
    free(s->nzval); free(s->rowind); free(s->colptr);
    free(s->R); free(s->C); free(s->rhs);
    free(s->perm_c); free(s->perm_r); free(s->etree);
    free(s);
}

static LIS_INT lis_slu_destroy(LIS_PRECON pc)
{
    if(pc)
    {
        lis_slu_release((LIS_SLU_FACTORS *)pc->user_data);
        pc->user_data=NULL;
    }
    return LIS_SUCCESS;
}

static LIS_INT lis_slu_create(LIS_SOLVER solver,LIS_PRECON pc)
{
#if defined(_COMPLEX) || defined(USE_COMPLEX) || defined(_LONG__DOUBLE)
    (void)solver; (void)pc;
    return LIS_ERR_NOT_IMPLEMENTED;
#else
    LIS_MATRIX m;
    LIS_SLU_FACTORS *s=NULL;
    int i,col,info,relax,panel_size;
    int_t k,nz,info_factor;
    int *next=NULL;
    double rowcnd,colcnd,amax;
    superlu_options_t options;
    SuperLUStat_t stat;
    GlobalLU_t glu;

    if(!solver || !pc || !solver->A) return LIS_ERR_ILL_ARG;
    m=solver->A;
    if(sizeof(LIS_SCALAR)!=sizeof(double) ||
       m->matrix_type!=LIS_MATRIX_CSR || m->nprocs!=1 ||
       m->n<1 || m->n!=m->gn || m->is!=0 || m->n>INT_MAX ||
       (uint64_t)m->nnz>INT_MAX || !m->ptr || !m->index || !m->value)
    {
        fprintf(stderr,"LIS_SUPERLU: only serial real-double CSR with 32-bit SuperLU indices is supported\n");
        return LIS_ERR_NOT_IMPLEMENTED;
    }
    s=(LIS_SLU_FACTORS *)calloc(1,sizeof(*s));
    if(!s) return LIS_OUT_OF_MEMORY;
    pc->user_data=s; /* lis_precon_destroy calls our cleanup on failure */
    s->n=(int)m->n;
    s->nnz=(int_t)m->nnz;
    if((LIS_INT)s->nnz!=m->nnz || m->ptr[0]!=0 ||
       m->ptr[m->n]!=m->nnz) return LIS_ERR_ILL_ARG;

    s->nzval=(double *)malloc((size_t)s->nnz*sizeof(double));
    s->rowind=(int_t *)malloc((size_t)s->nnz*sizeof(int_t));
    s->colptr=(int_t *)calloc((size_t)s->n+1,sizeof(int_t));
    s->perm_c=(int *)malloc((size_t)s->n*sizeof(int));
    s->perm_r=(int *)malloc((size_t)s->n*sizeof(int));
    s->etree=(int *)malloc((size_t)s->n*sizeof(int));
    s->R=(double *)malloc((size_t)s->n*sizeof(double));
    s->C=(double *)malloc((size_t)s->n*sizeof(double));
    s->rhs=(double *)malloc((size_t)s->n*sizeof(double));
    next=(int *)malloc((size_t)s->n*sizeof(int));
    if(!s->nzval || !s->rowind || !s->colptr || !s->perm_c ||
       !s->perm_r || !s->etree || !s->R || !s->C || !s->rhs || !next)
    {
        free(next);
        return LIS_OUT_OF_MEMORY;
    }

    /* CSR -> CSC, with no assumption that A is symmetric. */
    for(i=0;i<s->n;i++)
    {
        if(m->ptr[i]<0 || m->ptr[i+1]<m->ptr[i] ||
           m->ptr[i+1]>m->nnz)
        {
            free(next);
            return LIS_ERR_ILL_ARG;
        }
        for(k=m->ptr[i];k<m->ptr[i+1];k++)
        {
            LIS_INT c=m->index[k];
            if(c<0 || c>=s->n || !isfinite((double)m->value[k]))
            {
                free(next);
                return LIS_ERR_ILL_ARG;
            }
            s->colptr[c+1]++;
        }
    }
    for(i=0;i<s->n;i++)s->colptr[i+1]+=s->colptr[i];
    for(i=0;i<s->n;i++)next[i]=(int)s->colptr[i];
    for(i=0;i<s->n;i++)
        for(k=m->ptr[i];k<m->ptr[i+1];k++)
        {
            col=(int)m->index[k];
            nz=next[col]++;
            s->rowind[nz]=(int_t)i;
            s->nzval[nz]=(double)m->value[k];
        }
    free(next);
    next=NULL;
    dCreate_CompCol_Matrix(&s->A,s->n,s->n,s->nnz,
                           s->nzval,s->rowind,s->colptr,
                           SLU_NC,SLU_D,SLU_GE);
    s->have_A=1;

    /* SuperLU equilibration; keep both scalings for forward/transpose. */
    dgsequ(&s->A,s->R,s->C,&rowcnd,&colcnd,&amax,&info);
    if(info!=0)
    {
        fprintf(stderr,"LIS_SUPERLU: equilibration failed (info=%d)\n",info);
        return LIS_BREAKDOWN;
    }
    dlaqgs(&s->A,s->R,s->C,rowcnd,colcnd,amax,&s->equed);
    set_default_options(&options);
    options.ColPerm=COLAMD;
    options.DiagPivotThresh=1.0;
    options.Equil=NO; /* already equilibrated by dlaqgs */
    options.PrintStat=NO;
    options.IterRefine=NOREFINE;
    get_perm_c(options.ColPerm,&s->A,s->perm_c);
    sp_preorder(&options,&s->A,s->perm_c,s->etree,&s->AC);
    s->have_AC=1;
    relax=sp_ienv(2);
    panel_size=sp_ienv(1);
    memset(&glu,0,sizeof(glu));
    StatInit(&stat);
    info_factor=0;
    dgstrf(&options,&s->AC,relax,panel_size,s->etree,NULL,0,
           s->perm_c,s->perm_r,&s->L,&s->U,&glu,&stat,&info_factor);
    info=(int)info_factor;
    StatFree(&stat);
    if(s->L.Store)s->have_L=1;
    if(s->U.Store)s->have_U=1;
    if(info!=0)
    {
        fprintf(stderr,"LIS_SUPERLU: numerical factorization failed (info=%d)\n",info);
        return LIS_BREAKDOWN;
    }
    s->initialized=1;
    pc->A=solver->A;
    pc->is_copy=LIS_FALSE;
    fprintf(stderr,
        "LIS_SUPERLU: FACTORIZED n=%d nnz=%lld equil=%c (COLAMD, pivot=1.0)\n",
        s->n,(long long)s->nnz,s->equed);
    return LIS_SUCCESS;
#endif
}

static LIS_INT lis_slu_apply_inner(LIS_SOLVER solver,LIS_VECTOR b,
                                   LIS_VECTOR x,trans_t op)
{
#if defined(_COMPLEX) || defined(USE_COMPLEX) || defined(_LONG__DOUBLE)
    (void)solver;(void)b;(void)x;(void)op;
    return LIS_ERR_NOT_IMPLEMENTED;
#else
    LIS_SLU_FACTORS *s;
    SuperMatrix B;
    SuperLUStat_t stat;
    int i,info=0,has_rows,has_cols;
    if(!solver || !solver->precon || !b || !x) return LIS_ERR_ILL_ARG;
    s=(LIS_SLU_FACTORS *)solver->precon->user_data;
    if(!s || !s->initialized || b->n!=s->n || x->n!=s->n)
        return LIS_ERR_ILL_ARG;
    has_rows=(s->equed=='R' || s->equed=='B');
    has_cols=(s->equed=='C' || s->equed=='B');
    for(i=0;i<s->n;i++)
    {
        double v=(double)b->value[i];
        if(op==NOTRANS ? has_rows : has_cols)
            v*=op==NOTRANS ? s->R[i] : s->C[i];
        if(!isfinite(v))return LIS_BREAKDOWN;
        s->rhs[i]=v;
    }
    dCreate_Dense_Matrix(&B,s->n,1,s->rhs,s->n,SLU_DN,SLU_D,SLU_GE);
    StatInit(&stat);
    dgstrs(op,&s->L,&s->U,s->perm_c,s->perm_r,&B,&stat,&info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&B);
    if(info!=0)return LIS_BREAKDOWN;
    for(i=0;i<s->n;i++)
    {
        double v=s->rhs[i];
        if(op==NOTRANS ? has_cols : has_rows)
            v*=op==NOTRANS ? s->C[i] : s->R[i];
        if(!isfinite(v))return LIS_BREAKDOWN;
        x->value[i]=(LIS_SCALAR)v;
    }
    return LIS_SUCCESS;
#endif
}

static LIS_INT lis_slu_apply(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    return lis_slu_apply_inner(solver,b,x,NOTRANS);
}
static LIS_INT lis_slu_applyh(LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR x)
{
    return lis_slu_apply_inner(solver,b,x,TRANS);
}

LIS_INT lis_superlu_direct_register(void)
{
    return lis_precon_register_ex("superlu",lis_slu_create,lis_slu_apply,
                                  lis_slu_applyh,lis_slu_destroy);
}
