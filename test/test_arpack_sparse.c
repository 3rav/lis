/* Non-diagonal finite-element-like tridiagonal spring chain.
 * K=tridiag(-1,2,-1), M=I; exact lambda_j=2-2cos(j*pi/(n+1)).
 */
#include "lis_arpack.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define N 160
#define NEV 5
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef struct { double piv[N], upper[N], scratch[N]; int factorizations, solves; } chain;
static void factor(chain *c)
{
    int i;
    c->piv[0]=2.0;
    for(i=1;i<N;i++) {
        c->upper[i-1]=-1.0/c->piv[i-1];
        c->piv[i]=2.0 + c->upper[i-1];
    }
    c->factorizations++;
}
static int kapply(void *ctx, int n, const double *x, double *y)
{
    int i; (void)ctx;
    for(i=0;i<n;i++)y[i]=2*x[i]-(i?x[i-1]:0)-(i+1<n?x[i+1]:0);
    return 0;
}
static int mapply(void *ctx, int n, const double *x, double *y)
{
    int i; (void)ctx;
    for(i=0;i<n;i++)y[i]=x[i];
    return 0;
}
static int solve(void *ctx, int n, const double *x, double *y)
{
    chain *c=(chain *)ctx; int i;
    c->scratch[0]=x[0];
    for(i=1;i<n;i++)c->scratch[i]=x[i]+c->scratch[i-1]/c->piv[i-1];
    y[n-1]=c->scratch[n-1]/c->piv[n-1];
    for(i=n-2;i>=0;i--)y[i]=c->scratch[i]/c->piv[i]-c->upper[i]*y[i+1];
    c->solves++;
    return 0;
}
int main(void)
{
    chain c={0};
    lis_arpack_problem p={N,3,0.0,&c,kapply,mapply,solve};
    lis_arpack_options o={NEV,24,500,1e-12,"LM"};
    lis_arpack_stats s={0};
    double d[NEV],v[N*NEV],r[NEV];
    int i,j,rc;
    factor(&c);
    rc=lis_arpack_symmetric(&p,&o,d,v,r,&s);
    if(rc || s.nconv!=NEV || c.factorizations!=1 || c.solves!=s.solve_count)return 1;
    for(i=0;i<NEV;i++) {
        double exact=2.0-2.0*cos((i+1)*M_PI/(N+1));
        double best=1e100;
        for(j=0;j<NEV;j++)if(fabs(d[j]-exact)<best)best=fabs(d[j]-exact);
        if(best>1e-9 || r[i]>1e-8) {
            fprintf(stderr,"spring mode=%d expected=%.12g err=%g residual=%g\n",i,exact,best,r[i]);
            return 1;
        }
    }
    printf("PASS tridiagonal N=%d nev=%d iterations=%d solves=%d factorizations=%d\n",
           N,NEV,s.iterations,s.solve_count,c.factorizations);
    return 0;
}
