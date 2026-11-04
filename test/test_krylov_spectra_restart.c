/* Convergence regression for a clustered generalized spectrum.
 * The legacy half-subspace restart does not reliably resolve this
 * eight-value cluster in 220 cycles at ncv=12; the guarded retention
 * policy must solve the requested first four eigenpairs.  This test
 * does not use the internals of restart_basis or expose a new API.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 96
#define NEV 4
#define NCV 12

typedef struct {
    double a[N], b[N], sigma;
    int operations, solves;
} context;

static int apply_a(void *ptr, int n, const double *x, double *y)
{
    context *c = (context *)ptr;
    int i;
    if (n != N) return 1;
    ++c->operations;
    for (i=0; i<n; ++i) y[i] = c->a[i]*x[i];
    return 0;
}
static int apply_b(void *ptr, int n, const double *x, double *y)
{
    context *c = (context *)ptr;
    int i;
    if (n != N) return 1;
    for (i=0; i<n; ++i) y[i] = c->b[i]*x[i];
    return 0;
}
static int solve(void *ptr, int n, const double *x, double *y)
{
    context *c = (context *)ptr;
    int i;
    if (n != N) return 1;
    ++c->solves;
    for (i=0; i<n; ++i) {
        double d=c->a[i]-c->sigma*c->b[i];
        if (d == 0.0) return 1;
        y[i]=x[i]/d;
    }
    return 0;
}

static int run(int hard_locking, unsigned seed)
{
    context c;
    lis_krylov_problem problem = {0};
    lis_krylov_options opt;
    lis_krylov_stats stats;
    double values[NEV], vectors[N*NEV], residuals[NEV];
    int i, j, k, rc;
    memset(&c,0,sizeof(c));
    c.sigma = 0.90;
    for (i=0;i<N;++i) {
        double eig=i<8 ? 1.0 + 1e-5*i : 2.0 + 0.12*(i-8);
        c.b[i]=1.0+0.02*(i%5);
        c.a[i]=c.b[i]*eig;
    }
    problem.n=N;
    problem.transform=LIS_KRYLOV_SHIFT_INVERT_B;
    problem.sigma=c.sigma;
    problem.ctx=&c;
    problem.apply_a=apply_a;
    problem.apply_b=apply_b;
    problem.solve=solve;
    lis_krylov_options_default(&opt);
    opt.nev=NEV;
    opt.ncv=NCV;
    opt.maxiter=220;
    opt.tol=1e-9;
    opt.which="LM";
    opt.random_seed=seed;
    opt.hard_locking=hard_locking;
    rc=lis_krylov_symmetric(&problem,&opt,values,vectors,residuals,&stats);
    if (rc!=LIS_KRYLOV_SUCCESS || stats.nconv!=NEV ||
        stats.nreturned!=NEV || stats.restarts<1 ||
        stats.solve_count!=c.solves || c.solves>240) {
        fprintf(stderr,"FAIL guarded cluster: hardlock=%d seed=%u rc=%d "
                       "nconv=%d restarts=%d solves=%d\n",
                hard_locking,seed,rc,stats.nconv,stats.restarts,c.solves);
        return 1;
    }
    for(j=0;j<NEV;++j) {
        double nr=0.0,na=0.0,nb=0.0,metric=0.0,lambda=values[j];
        if(fabs(lambda-(1.0+1e-5*j))>1e-6 ||
           !isfinite(residuals[j]) || residuals[j]>opt.tol) {
            fprintf(stderr,"FAIL eigenvalue %d: %.17g residual %.17g\n",
                    j,lambda,residuals[j]);
            return 1;
        }
        for(i=0;i<N;++i) {
            double v=vectors[(size_t)j*N+i];
            double av=c.a[i]*v, bv=c.b[i]*v, r=av-lambda*bv;
            nr+=r*r; na+=av*av; nb+=bv*bv;
            metric+=v*bv;
        }
        if(fabs(metric-1.0)>1e-6 ||
           sqrt(nr)>1e-8*(sqrt(na)+fabs(lambda)*sqrt(nb))) {
            fprintf(stderr,"FAIL physical B metric/residual pair %d\n",j);
            return 1;
        }
        for(k=0;k<j;++k) {
            double dot=0.0;
            for(i=0;i<N;++i)
                dot+=c.b[i]*vectors[(size_t)k*N+i]*vectors[(size_t)j*N+i];
            if(fabs(dot)>1e-5) {
                fprintf(stderr,"FAIL B-orthogonality %d/%d %.17g\n",j,k,dot);
                return 1;
            }
        }
    }
    printf("PASS clustered B-metric restart hardlock=%d seed=%u "
           "cycles=%d restarts=%d solves=%d\n",
           hard_locking,seed,stats.iterations,stats.restarts,c.solves);
    return 0;
}

int main(void)
{
    return run(0,1)||run(0,2)||run(1,1)||run(1,2);
}
