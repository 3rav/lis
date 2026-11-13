/*
 * Standalone sparse approximate-Schur proof of concept.
 *
 * This file is NOT wired into LIS, and is NOT a replacement for
 * experimental serial RAS/Schur in PR #225.
 *
 * Builds a sparse CSR approximation
 *     S = A_GG - A_GI diag(A_II)^(-1) A_IG,
 * verifies its action without a dense intermediate, and exercises a
 * diagonal-preconditioned Richardson iteration on a 384-DOF interface.
 *
 * cc -O2 -std=c99 -Wall -Wextra -pedantic -o sparse_schur_poc \
 *     test/experimental_sparse_schur.c -lm && ./sparse_schur_poc
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n;
    size_t *ptr;
    int *col;
    double *val;
    size_t nnz, capacity;
} csr;

static void csr_free(csr *a)
{
    free(a->ptr);
    free(a->col);
    free(a->val);
    memset(a, 0, sizeof(*a));
}

static int csr_init(csr *a, int n)
{
    memset(a, 0, sizeof(*a));
    a->n=n;
    a->ptr=(size_t *)calloc((size_t)n+1,sizeof(size_t));
    return a->ptr!=NULL;
}

static int csr_append(csr *a,int col,double value)
{
    if(!isfinite(value) || col<0 || col>=a->n) return 0;
    if(a->nnz==a->capacity) {
        size_t cap=a->capacity ? 2*a->capacity : 128;
        int *new_col;
        double *new_val;
        if(cap<=a->capacity || cap>SIZE_MAX/sizeof(double) ||
           cap>SIZE_MAX/sizeof(int)) return 0;
        new_col=(int *)realloc(a->col,cap*sizeof(int));
        if(!new_col) return 0;
        a->col=new_col;
        new_val=(double *)realloc(a->val,cap*sizeof(double));
        if(!new_val) return 0;
        a->val=new_val;
        a->capacity=cap;
    }
    a->col[a->nnz]=col;
    a->val[a->nnz]=value;
    ++a->nnz;
    return 1;
}

/* Build an SPD two-block system with an interface of n_gamma unknowns.
 * Gamma rows: tridiagonal block with additional one-to-one interior link.
 * Interior rows: positive diagonal and one link to Gamma.
 * S has at most three nonzeros per row, independent of n_gamma.
 */
static int make_example(csr *a,int n_gamma)
{
    int i;
    if(!csr_init(a,2*n_gamma)) return 0;
    for(i=0;i<n_gamma;i++) {
        if(i>0 && !csr_append(a,i-1,-0.5)) return 0;
        if(!csr_append(a,i,4.0)) return 0;
        if(i+1<n_gamma && !csr_append(a,i+1,-0.5)) return 0;
        if(!csr_append(a,n_gamma+i,-0.6)) return 0;
        a->ptr[i+1]=a->nnz;
    }
    for(i=0;i<n_gamma;i++) {
        if(!csr_append(a,i,-0.6)) return 0;
        if(!csr_append(a,n_gamma+i,2.0)) return 0;
        a->ptr[n_gamma+i+1]=a->nnz;
    }
    return 1;
}

/* Sparse row accumulation: storage O(nnz(S)+n_gamma), not O(n_gamma^2).
 * In a real implementation, integer index overflow, row work and
 * cost of interface extraction also need large-scale investigation.
 */
static int integer_cmp(const void *ap,const void *bp)
{
    int a=*(const int *)ap,b=*(const int *)bp;
    return (a>b)-(a<b);
}

static int make_schur_sparse(const csr *a,int ng,const double *inv,csr *s)
{
    int i;
    int *marks=NULL,*touched=NULL;
    double *acc=NULL;
    int ok=0;
    if(!csr_init(s,ng)) return 0;
    marks=(int *)malloc((size_t)ng*sizeof(int));
    touched=(int *)malloc((size_t)ng*sizeof(int));
    acc=(double *)calloc((size_t)ng,sizeof(double));
    if(!marks || !touched || !acc) goto done;
    for(i=0;i<ng;i++) marks[i]=-1;
    for(i=0;i<ng;i++) {
        size_t k,q;
        int count=0,j;
        for(k=a->ptr[i];k<a->ptr[i+1];k++) {
            int col=a->col[k];
            if(col<ng) {
                if(marks[col]!=i) {
                    marks[col]=i;
                    acc[col]=0.0;
                    touched[count++]=col;
                }
                acc[col]+=a->val[k];
            } else {
                int interior=col-ng;
                double coeff=-a->val[k]*inv[interior];
                if(!isfinite(coeff)) goto done;
                for(q=a->ptr[col];q<a->ptr[col+1];q++) {
                    j=a->col[q];
                    if(j>=ng) continue;
                    if(marks[j]!=i) {
                        marks[j]=i;
                        acc[j]=0.0;
                        touched[count++]=j;
                    }
                    acc[j]+=coeff*a->val[q];
                }
            }
        }
        qsort(touched,(size_t)count,sizeof(int),integer_cmp);
        for(j=0;j<count;j++) {
            int col=touched[j];
            if(!isfinite(acc[col])) goto done;
            /* Exact zeros can be omitted from the sparse pattern. */
            if(acc[col]!=0.0 && !csr_append(s,col,acc[col])) goto done;
        }
        s->ptr[i+1]=s->nnz;
    }
    ok=1;
done:
    free(marks);
    free(touched);
    free(acc);
    return ok;
}

static void csr_matvec(const csr *a,const double *x,double *b)
{
    int i;
    for(i=0;i<a->n;i++) {
        size_t k;
        double v=0.0;
        for(k=a->ptr[i];k<a->ptr[i+1];k++)
            v+=a->val[k]*x[a->col[k]];
        b[i]=v;
    }
}

static int test_poc(void)
{
    const int ng=384;
    csr a={0},s={0};
    double *inv=NULL,*exact=NULL,*b=NULL,*x=NULL,*res=NULL,*diag=NULL;
    double bn=0.0,err=0.0,max_action_error=0.0;
    int i,iter=0,ok=0;
    inv=(double *)malloc((size_t)ng*sizeof(double));
    exact=(double *)malloc((size_t)ng*sizeof(double));
    b=(double *)malloc((size_t)ng*sizeof(double));
    x=(double *)calloc((size_t)ng,sizeof(double));
    res=(double *)malloc((size_t)ng*sizeof(double));
    diag=(double *)malloc((size_t)ng*sizeof(double));
    if(!inv || !exact || !b || !x || !res || !diag ||
       !make_example(&a,ng)) goto done;

    for(i=0;i<ng;i++) {
        size_t k;
        double diag_ii=0.0,row_max=0.0;
        for(k=a.ptr[ng+i];k<a.ptr[ng+i+1];k++) {
            double v=fabs(a.val[k]);
            if(v>row_max) row_max=v;
            if(a.col[k]==ng+i) diag_ii=a.val[k];
        }
        if(!(fabs(diag_ii)>1e-10*row_max) || !isfinite(diag_ii))
            goto done;
        inv[i]=1.0/diag_ii;
        if(!isfinite(inv[i])) goto done;
    }
    if(!make_schur_sparse(&a,ng,inv,&s)) goto done;

    for(i=0;i<ng;i++) exact[i]=1.0+0.1*sin(0.17*(double)i);
    csr_matvec(&s,exact,b);
    for(i=0;i<ng;i++) {
        size_t k,q;
        double direct=0.0;
        for(k=a.ptr[i];k<a.ptr[i+1];k++) {
            int j=a.col[k];
            if(j<ng) direct+=a.val[k]*exact[j];
            else {
                double inner=0.0;
                for(q=a.ptr[j];q<a.ptr[j+1];q++)
                    if(a.col[q]<ng) inner+=a.val[q]*exact[a.col[q]];
                direct-=a.val[k]*inv[j-ng]*inner;
            }
        }
        if(fabs(b[i]-direct)>max_action_error)
            max_action_error=fabs(b[i]-direct);
        for(k=s.ptr[i];k<s.ptr[i+1];k++)
            if(s.col[k]==i) diag[i]=s.val[k];
        if(!(diag[i]>0.0) || !isfinite(diag[i])) goto done;
        bn+=b[i]*b[i];
    }
    bn=sqrt(bn);
    if(!(bn>0.0)) goto done;

    /* Diagonal-preconditioned Richardson: for the intentionally
     * diagonally dominant illustrative system only.
     * Not an indefinite-KKT or contact-Newton solver.
     */
    for(iter=0;iter<150;iter++) {
        double rn=0.0;
        csr_matvec(&s,x,res);
        for(i=0;i<ng;i++) {
            res[i]=b[i]-res[i];
            rn+=res[i]*res[i];
        }
        err=sqrt(rn)/bn;
        if(!isfinite(err)) goto done;
        if(err<1.0e-10) break;
        for(i=0;i<ng;i++)
            x[i]+=res[i]/diag[i];
    }
    printf("SPARSE_SCHUR_POC ngamma=%d nnz=%zu dense_entries=%zu "
           "action_max_abs_error=%.3e iter=%d rel_residual=%.3e\n",
           ng,s.nnz,(size_t)ng*(size_t)ng,max_action_error,iter,err);
    ok=(ng>256 && s.nnz<(size_t)ng*(size_t)ng/8 &&
        max_action_error<1e-12 && err<1e-10);
done:
    if(!ok) fprintf(stderr,"SPARSE_SCHUR_POC FAILED\n");
    csr_free(&a);
    csr_free(&s);
    free(inv); free(exact); free(b); free(x); free(res); free(diag);
    return ok ? 0 : 1;
}

int main(void)
{
    return test_poc();
}
