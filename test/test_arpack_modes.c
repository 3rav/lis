/* Numerical regressions against exact generalized diagonal eigenpairs.
 * Build with: make -C test test_arpack_modes
 */
#include "lis_arpack.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 64
typedef struct {
    double k[N], m[N], kg[N], shifted[N];
    int failure, solves;
} diagonal;

static int kapply(void *ctx, int n, const double *x, double *y)
{
    diagonal *d = (diagonal *)ctx; int i;
    if (d->failure) return -1;
    for (i=0;i<n;i++) y[i]=d->k[i]*x[i];
    return 0;
}
static int mapply(void *ctx, int n, const double *x, double *y)
{
    diagonal *d = (diagonal *)ctx; int i;
    for (i=0;i<n;i++) y[i]=d->m[i]*x[i];
    return 0;
}
static int kgapply(void *ctx, int n, const double *x, double *y)
{
    diagonal *d = (diagonal *)ctx; int i;
    for (i=0;i<n;i++) y[i]=d->kg[i]*x[i];
    return 0;
}
static int solve(void *ctx, int n, const double *x, double *y)
{
    diagonal *d = (diagonal *)ctx; int i;
    d->solves++;
    for (i=0;i<n;i++) y[i]=x[i]/d->shifted[i];
    return 0;
}
static int check(int mode, double sigma, const double *expected, int nev)
{
    diagonal d={0};
    lis_arpack_problem p;
    lis_arpack_options o;
    lis_arpack_stats s;
    double ev[8]={0}, vec[N*8]={0}, res[8]={0};
    int i,j, status;
    for(i=0;i<N;i++) {
        double lambda = (i%2 == 0) ? 2.5+0.4*i : -(3.0+0.5*i);
        d.k[i]=10.0 + i;
        d.m[i]=1.0 + 0.025*i;
        d.kg[i]=d.k[i]/lambda;
        if (mode==3) d.k[i]=(4.0+3.0*i)*d.m[i];
        d.shifted[i]=d.k[i]-sigma*(mode==3?d.m[i]:d.kg[i]);
    }
    p.n=N; p.mode=mode; p.sigma=sigma; p.context=&d;
    p.stiffness=kapply; p.second=mode==3?mapply:kgapply; p.solve=solve;
    o.nev=nev; o.ncv=24; o.maxiter=600; o.tol=1e-12; o.which="LM";
    status=lis_arpack_symmetric(&p,&o,ev,vec,res,&s);
    if(status!=LIS_ARPACK_SUCCESS || s.nconv!=nev || s.solve_count!=d.solves) {
        fprintf(stderr,"mode=%d sigma=%g status=%d nconv=%d info=%d solves=%d/%d\n",
                mode,sigma,status,s.nconv,s.arpack_info,d.solves,s.solve_count);
        return 1;
    }
    for(j=0;j<nev;j++) {
        double best=1e10;
        for (i=0;i<nev;i++) {
            double e=fabs(ev[i]-expected[j]);
            if(e<best)best=e;
        }
        if(best>1e-7 || !(res[j]<1e-9)) {
            fprintf(stderr,"mode %d eigenpair %d: lambda=%.15g best=%.3g resid=%.3g\n",
                    mode,j,ev[j],best,res[j]);
            return 1;
        }
    }
    /* ARPACK's B-orthogonality: B=M for modal, B=K for buckling. */
    for(i=0;i<nev;i++) for(j=0;j<nev;j++) {
        int row;
        double dot=0;
        for(row=0;row<N;row++)
            dot+=vec[i*N+row]*vec[j*N+row]*(mode==3?d.m[row]:d.k[row]);
        if (fabs(dot-(i==j?1.0:0.0))>1e-7) {
            fprintf(stderr,"mode %d: B-orthogonality (%d,%d) %.12g\n",mode,i,j,dot);
            return 1;
        }
    }
    printf("PASS mode=%d sigma=%g nev=%d ncv=%d iterations=%d OP=%d B=%d solves=%d maxres=%.3g\n",
        mode,sigma,nev,o.ncv,s.iterations,s.op_count,s.b_count,s.solve_count,res[0]);
    /* Negative tests: fail fast, without returning success or garbage. */
    o.nev=N;
    if(lis_arpack_symmetric(&p,&o,ev,vec,res,&s)!=LIS_ARPACK_BAD_ARGUMENT)return 1;
    o.nev=nev; d.failure=1;
    if(lis_arpack_symmetric(&p,&o,ev,vec,res,&s)!=LIS_ARPACK_CALLBACK_ERROR)return 1;
    return 0;
}
int main(void)
{
    const double freq[]={4,7,10,13};
    const double buckle_pos[]={2.5,3.3,4.1};
    const double buckle_neg[]={-3.5,-4.5,-5.5};
    if(check(3,0.0,freq,4) || check(4,2.0,buckle_pos,3) ||
       check(4,-3.0,buckle_neg,3)) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
