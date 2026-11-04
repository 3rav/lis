/* Synthetic spectral-cluster experiment for restart retention.
 * Compile this exact fixture with baseline and candidate implementations.
 * CSV contains raw operator counts and observable convergence, not claims
 * about absolute performance on production matrices.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define DIM 96
#define NEV 4

typedef struct {
    double a[DIM], b[DIM], sigma;
    unsigned long callbacks_a, callbacks_b, callbacks_solve;
} ctx;

static int op_a(void *v, int n, const double *x, double *y)
{
    ctx *c = (ctx *)v; int i;
    if (n != DIM) return 1;
    ++c->callbacks_a;
    for (i=0; i<n; ++i) y[i] = c->a[i]*x[i];
    return 0;
}
static int op_b(void *v, int n, const double *x, double *y)
{
    ctx *c = (ctx *)v; int i;
    if (n != DIM) return 1;
    ++c->callbacks_b;
    for (i=0; i<n; ++i) y[i] = c->b[i]*x[i];
    return 0;
}
static int op_solve(void *v, int n, const double *x, double *y)
{
    ctx *c = (ctx *)v; int i;
    if (n != DIM) return 1;
    ++c->callbacks_solve;
    for (i=0; i<n; ++i) {
        double d = c->a[i] - c->sigma*c->b[i];
        if (fabs(d)<1e-12) return 1;
        y[i] = x[i]/d;
    }
    return 0;
}

static int run(int family, int transform, int ncv, int locking, int seed)
{
    ctx c;
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_stats s;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    double values[NEV], vecs[DIM*NEV], residuals[NEV];
    clock_t begin, end;
    int i, j, rc;
    memset(&c, 0, sizeof(c));
    for (i=0; i<DIM; ++i) {
        double lambda;
        if (family == 0)
            lambda = 1.0 + 0.1*i; /* well separated */
        else if (family == 1)
            lambda = i < 8 ? 1.0 + 1e-5*i : 2.0 + 0.12*(i-8);
        else if (family == 2)
            lambda = i < 12 ? 1.0 : 2.0 + 0.12*(i-12);
        else if (family == 3)
            lambda = i < 12 ? 1.0 + 1e-4*i : 1.1 + 0.02*(i-12);
        else
            lambda = i < 20 ? 1.0 + 0.002*i : 1.1 + 0.05*(i-20);
        c.b[i] = 1.0 + 0.02*(i%5);
        c.a[i] = c.b[i] * lambda;
    }
    c.sigma = 0.90;
    p.n = DIM;
    p.ctx = &c;
    p.transform = transform;
    p.sigma = transform == LIS_KRYLOV_STANDARD ? 0.0 : c.sigma;
    p.apply_a = op_a;
    p.apply_b = op_b;
    p.solve = op_solve;
    lis_krylov_options_default(&o);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    o.nev = NEV;
    o.ncv = ncv;
    o.maxiter = 220;
    o.tol = 1e-9;
    o.which = transform == LIS_KRYLOV_STANDARD ? "SA" : "LM";
    o.random_seed = (unsigned)seed;
    o.hard_locking = locking;
    begin = clock();
    rc = lis_krylov_symmetric_ex(&p, &o, &ext, NULL,
                                 values, vecs, residuals, &s, &extra);
    end = clock();
    if (rc != LIS_KRYLOV_SUCCESS && rc != LIS_KRYLOV_NO_CONVERGENCE) {
        fprintf(stderr,"FAIL solver family=%d t=%d ncv=%d lock=%d seed=%d rc=%d\n",
                family,transform,ncv,locking,seed,rc);
        return 1;
    }
    if (s.nconv > s.nreturned || s.nreturned > NEV ||
        s.solve_count != (int)c.callbacks_solve) return 1;
    for (j=0; j<s.nreturned; ++j) {
        double nx=0.0, nr=0.0, na=0.0, nb=0.0, denom;
        for(i=0;i<DIM;++i) {
            double v=vecs[(size_t)j*DIM+i];
            double ax=c.a[i]*v;
            double bx=(transform == LIS_KRYLOV_STANDARD ? v : c.b[i]*v);
            double r=ax-values[j]*bx;
            nr+=r*r; na+=ax*ax; nb+=bx*bx; nx+=v*v;
        }
        denom=sqrt(na)+fabs(values[j])*sqrt(nb);
        if (o.residual_scale*sqrt(nx)>denom) denom=o.residual_scale*sqrt(nx);
        if (s.nconv == NEV && sqrt(nr)>2e-7*fmax(denom,1.0)) {
            fprintf(stderr,"FAIL physical residual family %d lambda %g\n", family,values[j]);
            return 1;
        }
        if(!isfinite(values[j]) || !isfinite(residuals[j])) return 1;
    }
    printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lu,%lu,%lu,%.6f\n",
           family,transform,ncv,locking,seed,rc,s.nconv,s.nreturned,
           s.iterations,s.restarts,c.callbacks_a,c.callbacks_b,
           c.callbacks_solve,(double)(end-begin)/CLOCKS_PER_SEC);
    return 0;
}

int main(void)
{
    int f,t,ncv,lock,seed;
    puts("family,transform,ncv,lock,seed,status,nconv,nreturned,iterations,restarts,apply_a,apply_b,solve,cpu_seconds");
    for(f=0; f<5; ++f)
      for(t=0; t<3; ++t)
        for(ncv=12; ncv<=18; ncv+=6)
          for(lock=0; lock<=1; ++lock)
            for(seed=1; seed<=2; ++seed)
              if(run(f,t,ncv,lock,seed)) return 1;
    return 0;
}
