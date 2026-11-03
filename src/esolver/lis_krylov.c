/* Native thick-restarted symmetric Krylov/Rayleigh-Ritz eigensolver.
 * BSD-3-Clause: Copyright (c) 2026 LIS contributors.
 * Operator images are reused through every restart. No external eigensolver.
 */
#include "lis_krylov.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const lis_krylov_problem *p;
    lis_krylov_stats *stats;
    int n, ncv;
    double *v, *bv, *w, *t, *q, *diag;
    double *x, *bx, *rhs, *scratch;
    int *order;
} workspace;

static double inner(int n, const double *a, const double *b)
{
    double sum=0., compensation=0.;
    int i;
    for (i=0;i<n;i++) {
        double value=a[i]*b[i]-compensation;
        double next=sum+value;
        compensation=(next-sum)-value;
        sum=next;
    }
    return sum;
}
static double stable_norm(int n, const double *x)
{
    double scale=0., sumsq=1.;
    int i;
    for (i=0;i<n;i++) {
        double a=fabs(x[i]);
        if (!isfinite(a)) return INFINITY;
        if (a == 0.) continue;
        if (scale < a) {
            double ratio=scale/a;
            sumsq=1.+sumsq*ratio*ratio;
            scale=a;
        } else {
            double ratio=a/scale;
            sumsq+=ratio*ratio;
        }
    }
    return scale==0. ? 0. : scale*sqrt(sumsq);
}
static int metric(workspace *s, const double *x, double *y)
{
    const lis_krylov_problem *p=s->p;
    int i;
    if (p->transform==LIS_KRYLOV_STANDARD) {
        for(i=0;i<s->n;i++)y[i]=x[i];
    } else {
        s->stats->metric_count++;
        if ((p->transform==LIS_KRYLOV_SHIFT_INVERT_A ?
             p->apply_a : p->apply_b)(p->ctx,s->n,x,y)) return -1;
    }
    for(i=0;i<s->n;i++)if(!isfinite(y[i]))return -2;
    return 0;
}
static int op(workspace *s, const double *x, double *y)
{
    const lis_krylov_problem *p=s->p;
    int i,err;
    s->stats->op_count++;
    if (p->transform==LIS_KRYLOV_STANDARD)
        err=p->apply_a(p->ctx,s->n,x,y);
    else {
        err=metric(s,x,s->rhs);
        if(err)return err;
        s->stats->solve_count++;
        err=p->solve(p->ctx,s->n,s->rhs,y);
    }
    if(err)return -1;
    for(i=0;i<s->n;i++)if(!isfinite(y[i]))return -2;
    return 0;
}
static int orthonormalize(workspace *s, int m, double *x, double *bx)
{
    int pass,j,i;
    double square, scale, initial_square;
    {int err=metric(s,x,bx);if(err)return err;}
    initial_square=inner(s->n,x,bx);
    if(!isfinite(initial_square) || initial_square <= 0.)return -2;
    for (pass=0;pass<2;pass++) {
        for(j=0;j<m;j++) {
            const double *vj=s->v+(size_t)j*s->n;
            const double *bvj=s->bv+(size_t)j*s->n;
            double dot=inner(s->n,bvj,x);
            if (!isfinite(dot)) return -2;
            for(i=0;i<s->n;i++) {
                x[i]-=dot*vj[i];
            }
        }
        /* Reapply the metric to avoid cancellation in B*x when a
         * residual direction is almost in the retained subspace. */
        {int err=metric(s,x,bx);if(err)return err;}
    }
    square=inner(s->n,x,bx);
    if(!isfinite(square))return -2;
    if(square<=1e-27 || square < initial_square*1e-24)return 1;
    scale=1./sqrt(square);
    for(i=0;i<s->n;i++) {
        x[i]*=scale;
        bx[i]*=scale;
    }
    for(j=0;j<m;j++)
        if(fabs(inner(s->n,s->bv+(size_t)j*s->n,x))>1e-8)
            return 1;
    return 0;
}
static int add_vector(workspace *s, int m, double *x)
{
    int status=orthonormalize(s,m,x,s->bx);
    if(status)return status;
    memcpy(s->v+(size_t)m*s->n,x,(size_t)s->n*sizeof(double));
    memcpy(s->bv+(size_t)m*s->n,s->bx,(size_t)s->n*sizeof(double));
    {int err=op(s,x,s->w+(size_t)m*s->n);if(err)return err;}
    return 0;
}
static void seed_vector(int n, double *x, unsigned seed)
{
    unsigned a=seed ^ 0x6d2b79f5u;
    int i;
    for(i=0;i<n;i++) {
        a ^= a<<13; a ^= a>>17; a ^= a<<5;
        x[i]=((double)(a & 0xffffu)/32768.)-1.;
    }
}
/* Symmetric Jacobi diagonalization. Columns of q hold eigenvectors. */
static int diagonalize(int m, double *t, double *q, double *diag, int stride)
{
    int sweep,p,r,i;
    for(i=0;i<m;i++)for(p=0;p<m;p++)q[(size_t)i*stride+p]=(i==p);
    for(sweep=0;sweep<80;sweep++) {
        double max=0., norm=0.;
        for(p=0;p<m;p++)for(r=p+1;r<m;r++) {
            double a=fabs(t[(size_t)p*stride+r]);
            if(a>max)max=a;
        }
        for(p=0;p<m;p++)norm+=fabs(t[(size_t)p*stride+p]);
        if(max<=2e-15*(norm+DBL_MIN))break;
        for(p=0;p<m;p++)for(r=p+1;r<m;r++) {
            double apr=t[(size_t)p*stride+r];
            double app=t[(size_t)p*stride+p], arr=t[(size_t)r*stride+r];
            double tau, angle,c,s;
            if(fabs(apr)<=1e-16*(fabs(app)+fabs(arr)+DBL_MIN))continue;
            tau=(arr-app)/(2.*apr);
            angle=copysign(1.,tau)/(fabs(tau)+hypot(1.,tau));
            c=1./sqrt(1.+angle*angle);
            s=angle*c;
            t[(size_t)p*stride+p]=app-angle*apr;
            t[(size_t)r*stride+r]=arr+angle*apr;
            t[(size_t)p*stride+r]=0.;
            t[(size_t)r*stride+p]=0.;
            for(i=0;i<m;i++) {
                double aip,air;
                if(i!=p&&i!=r) {
                    aip=t[(size_t)i*stride+p];air=t[(size_t)i*stride+r];
                    t[(size_t)i*stride+p]=t[(size_t)p*stride+i]=c*aip-s*air;
                    t[(size_t)i*stride+r]=t[(size_t)r*stride+i]=s*aip+c*air;
                }
                aip=q[(size_t)i*stride+p];air=q[(size_t)i*stride+r];
                q[(size_t)i*stride+p]=c*aip-s*air;
                q[(size_t)i*stride+r]=s*aip+c*air;
            }
        }
    }
    if(sweep==80)return -1;
    for(i=0;i<m;i++)diag[i]=t[(size_t)i*stride+i];
    return 0;
}
static int earlier(double a, double b, const char *which)
{
    if (which[0]=='L' && which[1]=='M') return fabs(a)>fabs(b);
    if (which[0]=='S' && which[1]=='M') return fabs(a)<fabs(b);
    if (which[0]=='L' && which[1]=='A') return a>b;
    return a<b;
}
static void sort_order(int m, const double *diag, int *order,const char *which)
{
    int i,j;
    for(i=0;i<m;i++)order[i]=i;
    for(i=1;i<m;i++) {
        int key=order[i];
        j=i-1;
        while(j>=0 && earlier(diag[key],diag[order[j]],which)) {
            order[j+1]=order[j];j--;
        }
        order[j+1]=key;
    }
}
static void rotated(workspace *s, const double *basis, int m,
                    const double *q, int idx, double *out)
{
    int i,j;
    for(i=0;i<s->n;i++)out[i]=0.;
    for(j=0;j<m;j++) {
        double c=q[(size_t)j*s->ncv+idx];
        const double *col=basis+(size_t)j*s->n;
        for(i=0;i<s->n;i++)out[i]+=c*col[i];
    }
}
static int projected_eigenpair(workspace *s, int m, int idx,
        double mu, double *lambda, double *rel_resid)
{
    const lis_krylov_problem *p=s->p;
    int i;
    double na,nb,nr,den;
    rotated(s,s->v,m,s->q,idx,s->x);
    if (p->transform==LIS_KRYLOV_STANDARD) {
        *lambda=mu;
    } else if (p->transform==LIS_KRYLOV_SHIFT_INVERT_B) {
        if(fabs(mu)<DBL_MIN) return 1;
        *lambda=p->sigma+1./mu;
    } else {
        if(fabs(mu-1.)<DBL_EPSILON) return 1;
        *lambda=p->sigma*mu/(mu-1.);
    }
    if(!isfinite(*lambda))return 1;
    if(p->apply_a(p->ctx,s->n,s->x,s->bx))return -1;
    na=stable_norm(s->n,s->bx);
    if(p->transform==LIS_KRYLOV_STANDARD) {
        memcpy(s->rhs,s->x,(size_t)s->n*sizeof(double));
    } else if(p->apply_b(p->ctx,s->n,s->x,s->rhs))return -1;
    nb=stable_norm(s->n,s->rhs);
    for(i=0;i<s->n;i++)s->bx[i]-=*lambda*s->rhs[i];
    nr=stable_norm(s->n,s->bx);
    den=na+fabs(*lambda)*nb;
    *rel_resid=den>0. ? nr/den : nr;
    if(!isfinite(*rel_resid))return 1;
    return 0;
}
static void restart_basis(workspace *s, int m, int retain)
{
    int i,j,idx;
    double *sets[3]={s->v,s->bv,s->w};
    for(idx=0;idx<3;idx++) {
        for(j=0;j<retain;j++) {
            rotated(s,sets[idx],m,s->q,s->order[j],s->scratch+(size_t)j*s->n);
        }
        for(j=0;j<retain;j++) {
            for(i=0;i<s->n;i++)
                sets[idx][(size_t)j*s->n+i]=s->scratch[(size_t)j*s->n+i];
        }
    }
}
int lis_krylov_symmetric(const lis_krylov_problem *p,
                         const lis_krylov_options *o,
                         double *values,double *vectors,double *residuals,
                         lis_krylov_stats *stats)
{
    workspace s={0};
    size_t count, sub;
    int m=0, iter, ncv, nev, maxiter, rc=LIS_KRYLOV_NO_CONVERGENCE, i,j;
    double tol;
    const char *which;
    if(stats)memset(stats,0,sizeof(*stats));
    if(!p||!o||!values||!vectors||!residuals||!stats || p->n<3 ||
       o->nev<1||o->nev>=p->n || !p->apply_a ||
       p->transform<LIS_KRYLOV_STANDARD||p->transform>LIS_KRYLOV_SHIFT_INVERT_A||
       ((p->transform!=LIS_KRYLOV_STANDARD)&&(!p->apply_b||!p->solve)) ||
       (p->transform==LIS_KRYLOV_SHIFT_INVERT_A && p->sigma==0.) ||
       !isfinite(p->sigma)||!isfinite(o->tol)||o->tol<0.||o->maxiter<0)
        return LIS_KRYLOV_BAD_ARGUMENT;
    which=o->which?o->which:"LM";
    if(strcmp(which,"LM")&&strcmp(which,"SM")&&strcmp(which,"LA")&&strcmp(which,"SA"))
        return LIS_KRYLOV_BAD_ARGUMENT;
    nev=o->nev;maxiter=o->maxiter?o->maxiter:300;
    tol=o->tol?o->tol:1e-10;
    ncv=o->ncv?o->ncv:(nev<((p->n-8)/3)?3*nev+8:p->n);
    if(ncv<=nev||ncv>p->n||p->n>INT_MAX/8 || ncv>INT_MAX/8) return LIS_KRYLOV_BAD_ARGUMENT;
    count=(size_t)p->n*(size_t)ncv;
    sub=(size_t)ncv*(size_t)ncv;
    if (count>SIZE_MAX/sizeof(double)||sub>SIZE_MAX/sizeof(double))return LIS_KRYLOV_BAD_ARGUMENT;
    s.n=p->n;s.ncv=ncv;s.p=p;s.stats=stats;
#define ALLOC(ptr,nn,typ) do { (ptr)=(typ*)calloc((nn),sizeof(typ)); if(!(ptr)){rc=LIS_KRYLOV_NO_MEMORY;goto done;} } while(0)
    ALLOC(s.v,count,double);
    ALLOC(s.bv,count,double);
    ALLOC(s.w,count,double);
    ALLOC(s.scratch,count,double);
    ALLOC(s.t,sub,double);
    ALLOC(s.q,sub,double);
    ALLOC(s.diag,ncv,double);
    ALLOC(s.order,ncv,int);
    ALLOC(s.x,s.n,double);
    ALLOC(s.bx,s.n,double);
    ALLOC(s.rhs,s.n,double);
#undef ALLOC
    seed_vector(s.n,s.x,1u);
    rc=add_vector(&s,m,s.x);
    if(rc){rc=rc==-1?LIS_KRYLOV_CALLBACK_ERROR:LIS_KRYLOV_BREAKDOWN;goto done;}
    m=1;
    rc=LIS_KRYLOV_NO_CONVERGENCE;
    for(iter=0;iter<maxiter;iter++) {
        int oldm=m,expand=0;
        stats->iterations=iter+1;
        for(i=0;i<m;i++)for(j=i;j<m;j++) {
            double a=inner(s.n,s.bv+(size_t)i*s.n,s.w+(size_t)j*s.n);
            double b=inner(s.n,s.bv+(size_t)j*s.n,s.w+(size_t)i*s.n);
            s.t[(size_t)i*ncv+j]=s.t[(size_t)j*ncv+i]=0.5*(a+b);
        }
        if(diagonalize(m,s.t,s.q,s.diag,ncv)) {rc=LIS_KRYLOV_BREAKDOWN;goto done;}
        sort_order(m,s.diag,s.order,which);
        stats->nconv=0;stats->max_residual=0.;
        for(j=0;j<m&&j<nev;j++) {
            double lambda,rel;
            int idx=s.order[j];
            int e=projected_eigenpair(&s,m,idx,s.diag[idx],&lambda,&rel);
            if(e<0){rc=LIS_KRYLOV_CALLBACK_ERROR;goto done;}
            if(e==0 && rel<=tol) {
                int pos=stats->nconv++;
                values[pos]=lambda;
                residuals[pos]=rel;
                if(rel>stats->max_residual)stats->max_residual=rel;
                memcpy(vectors+(size_t)pos*s.n,s.x,(size_t)s.n*sizeof(double));
            }
        }
        if(stats->nconv==nev){rc=LIS_KRYLOV_SUCCESS;goto done;}
        if(m>=ncv) {
            int keep=m/2;
            if(keep<nev)keep=nev;
            if(keep>=m)keep=m-1;
            restart_basis(&s,m,keep);
            stats->restarts++;
            m=keep;
            continue;
        }
        /* Expand with transformed residual directions of targeted Ritz pairs.
         * The new subspace is B-orthogonalized twice, and all operator images
         * are cached for the next projection and for thick restarts.
         */
        for(j=0;j<oldm && j<nev && m<ncv;j++) {
            int idx=s.order[j];
            double mu=s.diag[idx];
            rotated(&s,s.w,oldm,s.q,idx,s.x);
            rotated(&s,s.v,oldm,s.q,idx,s.bx);
            for(i=0;i<s.n;i++)s.x[i]-=mu*s.bx[i];
            /* Tiny transformed residuals are numerically dependent:
             * do not normalize roundoff into a spurious Krylov vector. */
            if(stable_norm(s.n,s.x) < 1e-13 * fmax(1.,fabs(mu)))continue;
            /* add_vector uses bx as metric workspace; x is disjoint. */
            {
                int a=add_vector(&s,m,s.x);
                if(a<0){rc=a==-1?LIS_KRYLOV_CALLBACK_ERROR:LIS_KRYLOV_BREAKDOWN;goto done;}
                if(!a){m++;expand++;}
            }
        }
        if(!expand && m<ncv) {
            int trial;
            for(trial=0;trial<6 && !expand;trial++) {
                seed_vector(s.n,s.x,(unsigned)(iter*7+trial+2));
                {
                    int a=add_vector(&s,m,s.x);
                    if(a<0){rc=a==-1?LIS_KRYLOV_CALLBACK_ERROR:LIS_KRYLOV_BREAKDOWN;goto done;}
                    if(!a){m++;expand++;}
                }
            }
            if(!expand){rc=LIS_KRYLOV_BREAKDOWN;goto done;}
        }
    }
    rc=LIS_KRYLOV_NO_CONVERGENCE;
done:
    free(s.v);free(s.bv);free(s.w);free(s.scratch);
    free(s.t);free(s.q);free(s.diag);free(s.order);
    free(s.x);free(s.bx);free(s.rhs);
    return rc;
}
