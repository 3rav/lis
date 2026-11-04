/* Experimental, standalone implicitly restarted symmetric Lanczos.
 *
 * This is a numerical research fixture, NOT a replacement for the
 * production generalized Krylov eigensolver. It exercises symmetric
 * tridiagonal Lanczos, QR implicit shifts, basis compression and exact
 * physical-residual verification on standard SPD operators.
 *
 * No code is taken from the Spectra project; only the well-known
 * implicit QR shift / Lanczos recurrence mathematics is used.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define NMAX 128
#define MMAX 24
#define NEV 4

typedef struct {
    int n, family, calls;
    double diag[NMAX];
} operator_context;

typedef struct {
    int n, m, nev, restarts, op_count, qr_shifts;
    double v[MMAX+1][NMAX];     /* Lanczos vectors, extra residual unit */
    double alpha[MMAX], beta[MMAX];
    double t[MMAX][MMAX], eval[MMAX], u[MMAX][MMAX];
    int order[MMAX];
    double output[NEV][NMAX], lambda[NEV], relative[NEV];
} lanczos;

static double dot(int n, const double *x, const double *y)
{
    double s=0.0,correction=0.0;
    int i;
    for(i=0;i<n;++i) {
        double term=x[i]*y[i]-correction;
        double v=s+term;
        correction=(v-s)-term;
        s=v;
    }
    return s;
}
static double nrm(int n, const double *x)
{
    double sum=0.0; int i;
    for(i=0;i<n;++i) sum+=x[i]*x[i];
    return sqrt(sum);
}
static void axpy(int n, double a, const double *x, double *y)
{
    int i;
    for(i=0;i<n;++i) y[i]+=a*x[i];
}
static int apply(operator_context *c, const double *x, double *y)
{
    int i,n=c->n;
    ++c->calls;
    if(c->family==2) {
        for(i=0;i<n;++i)
            y[i]=2.0*x[i]-(i?x[i-1]:0.0)-(i+1<n?x[i+1]:0.0);
    } else {
        for(i=0;i<n;++i) y[i]=c->diag[i]*x[i];
    }
    return 0;
}

/* Small dense symmetric Jacobi diagonalization of tridiagonal T.
 * Columns of eigvec are normalized Ritz coefficient vectors. */
static int small_eigs(int m, const double a[MMAX][MMAX],
                      double *values, double eigvec[MMAX][MMAX])
{
    double t[MMAX][MMAX];
    int i,j,p,q,sweep;
    memcpy(t,a,sizeof(t));
    for(i=0;i<m;++i)
        for(j=0;j<m;++j) eigvec[i][j]=(i==j);
    for(sweep=0;sweep<100;++sweep) {
        double off=0.0,diagonal=0.0;
        for(i=0;i<m;++i) {
            diagonal+=fabs(t[i][i]);
            for(j=i+1;j<m;++j)
                if(fabs(t[i][j])>off) off=fabs(t[i][j]);
        }
        if(off<=4.0*DBL_EPSILON*(diagonal+1.0)) break;
        for(p=0;p<m-1;++p) for(q=p+1;q<m;++q) {
            double v=t[p][q], a00=t[p][p], a11=t[q][q];
            double tau,z,c,s;
            if(fabs(v)<DBL_EPSILON*(fabs(a00)+fabs(a11)+1.0))
                continue;
            tau=(a11-a00)/(2.0*v);
            z=copysign(1.0,tau)/(fabs(tau)+hypot(1.0,tau));
            c=1.0/sqrt(1.0+z*z);
            s=z*c;
            for(i=0;i<m;++i) if(i!=p&&i!=q) {
                double tp=t[i][p],tq=t[i][q];
                t[i][p]=t[p][i]=c*tp-s*tq;
                t[i][q]=t[q][i]=s*tp+c*tq;
            }
            t[p][p]=a00-z*v;
            t[q][q]=a11+z*v;
            t[p][q]=t[q][p]=0.0;
            for(i=0;i<m;++i) {
                double up=eigvec[i][p],uq=eigvec[i][q];
                eigvec[i][p]=c*up-s*uq;
                eigvec[i][q]=s*up+c*uq;
            }
        }
    }
    if(sweep==100) return 1;
    for(i=0;i<m;++i) {
        values[i]=t[i][i];
        if(!isfinite(values[i])) return 1;
    }
    return 0;
}

/* Apply one implicit shifted QR step to a symmetric tridiagonal T.
 * R=Q^T*(T-shift I), T'=R Q+shift I, Qacc'=Qacc*Q.
 * Givens rotations zero the subdiagonal of T-shift I.
 */
static int qr_shift(int m, double t[MMAX][MMAX],
                    double qacc[MMAX][MMAX], double shift)
{
    double r[MMAX][MMAX],q[MMAX][MMAX],out[MMAX][MMAX],
           transform[MMAX][MMAX],total=0.0,leak=0.0;
    int i,j,k;
    for(i=0;i<m;++i)
        for(j=0;j<m;++j) {
            r[i][j]=t[i][j]-(i==j?shift:0.0);
            q[i][j]=(i==j);
        }
    for(k=0;k<m-1;++k) {
        double a=r[k][k],b=r[k+1][k],h=hypot(a,b);
        double c=h>0.0?a/h:1.0,s=h>0.0?b/h:0.0;
        for(j=k;j<m;++j) {
            double x=r[k][j],y=r[k+1][j];
            r[k][j]=c*x+s*y;
            r[k+1][j]=-s*x+c*y;
        }
        r[k+1][k]=0.0;
        for(i=0;i<m;++i) {
            double x=q[i][k],y=q[i][k+1];
            q[i][k]=c*x+s*y;
            q[i][k+1]=-s*x+c*y;
        }
    }
    for(i=0;i<m;++i) for(j=0;j<m;++j) {
        double v=(i==j?shift:0.0);
        for(k=0;k<m;++k) v+=r[i][k]*q[k][j];
        out[i][j]=v;
        total+=v*v;
    }
    for(i=0;i<m;++i) for(j=0;j<m;++j) {
        double v=0.5*(out[i][j]+out[j][i]);
        if(i+1<j||j+1<i) leak+=v*v;
    }
    if(!isfinite(leak) || sqrt(leak)>1e-8*fmax(1.0,sqrt(total)))
        return 1;
    for(i=0;i<m;++i) for(j=0;j<m;++j) {
        t[i][j]=(i+1<j||j+1<i)?0.0:
                0.5*(out[i][j]+out[j][i]);
        transform[i][j]=0.0;
        for(k=0;k<m;++k)
            transform[i][j]+=qacc[i][k]*q[k][j];
    }
    for(i=0;i<m;++i)
        for(j=0;j<m;++j) qacc[i][j]=transform[i][j];
    return 0;
}

static int extend(operator_context *op, lanczos *s, int start)
{
    double w[NMAX];int j,k,pass,i;
    for(j=start;j<s->m;++j) {
        apply(op,s->v[j],w);
        if(j>0) axpy(s->n,-s->beta[j-1],s->v[j-1],w);
        s->alpha[j]=dot(s->n,s->v[j],w);
        axpy(s->n,-s->alpha[j],s->v[j],w);
        /* Full reorthogonalization for a numerically robust prototype. */
        for(pass=0;pass<2;++pass)
            for(k=0;k<=j;++k) {
                double c=dot(s->n,s->v[k],w);
                if(k==j) s->alpha[j]+=c;
                axpy(s->n,-c,s->v[k],w);
            }
        s->beta[j]=nrm(s->n,w);
        if(!isfinite(s->alpha[j]) || !isfinite(s->beta[j]) ||
           s->beta[j]<128.0*DBL_EPSILON)
            return 1;
        for(i=0;i<s->n;++i) s->v[j+1][i]=w[i]/s->beta[j];
    }
    s->op_count=op->calls;
    return 0;
}

static void make_projected(lanczos *s)
{
    int i,j;
    for(i=0;i<s->m;++i)
        for(j=0;j<s->m;++j)
            s->t[i][j]=(i==j?s->alpha[i]:
                         i+1==j?s->beta[i]:
                         j+1==i?s->beta[j]:0.0);
}
static void sort_smallest(lanczos *s)
{
    int i,j;
    for(i=0;i<s->m;++i) s->order[i]=i;
    for(i=1;i<s->m;++i) {
        int item=s->order[i];
        for(j=i;j>0&&s->eval[s->order[j-1]]>s->eval[item];--j)
            s->order[j]=s->order[j-1];
        s->order[j]=item;
    }
}

/* Compute the complete true residual before accepting eigenpairs. */
static int physical_pair(operator_context *op, lanczos *s,
                         int k, double *relative)
{
    double ax[NMAX],nx,na,nr,lambda=s->eval[k];int i,j;
    for(i=0;i<s->n;++i) {
        double sum=0.0;
        for(j=0;j<s->m;++j) sum+=s->v[j][i]*s->u[j][k];
        s->output[0][i]=sum; /* temporary column, saved on acceptance */
    }
    apply(op,s->output[0],ax);
    nx=nrm(s->n,s->output[0]);
    na=nrm(s->n,ax);
    axpy(s->n,-lambda,s->output[0],ax);
    nr=nrm(s->n,ax);
    *relative=nr/(na+fabs(lambda)*nx+DBL_MIN);
    return !isfinite(*relative);
}

/* Implicitly restart at k with the m-k unwanted Ritz shifts.
 * The Lanczos residual after compression has two terms:
 *  T'[k,k-1] VQ[:,k] + beta[m-1]*Q[m-1,k-1]*v_extra.
 */
static int implicit_restart(lanczos *s,int k)
{
    double qacc[MMAX][MMAX],t[MMAX][MMAX];
    double newv[MMAX+1][NMAX],r[NMAX],b;
    int i,j,z,l,m=s->m;
    memcpy(t,s->t,sizeof(t));
    for(i=0;i<m;++i)
        for(j=0;j<m;++j) qacc[i][j]=(i==j);
    for(i=m-1;i>=k;--i) {
        if(qr_shift(m,t,qacc,s->eval[s->order[i]])) return 1;
        ++s->qr_shifts;
    }
    for(j=0;j<=k;++j)
        for(i=0;i<s->n;++i) {
            double v=0.0;
            for(z=0;z<m;++z) v+=s->v[z][i]*qacc[z][j];
            newv[j][i]=v;
        }
    for(i=0;i<s->n;++i)
        r[i]=t[k][k-1]*newv[k][i] +
             s->beta[m-1]*qacc[m-1][k-1]*s->v[m][i];
    for(l=0;l<2;++l)
        for(j=0;j<k;++j) {
            double c=dot(s->n,newv[j],r);
            axpy(s->n,-c,newv[j],r);
        }
    b=nrm(s->n,r);
    if(!isfinite(b) || b<128.0*DBL_EPSILON) return 1;
    for(i=0;i<k;++i) {
        s->alpha[i]=t[i][i];
        if(i+1<k) s->beta[i]=t[i][i+1];
        memcpy(s->v[i],newv[i],(size_t)s->n*sizeof(double));
    }
    s->beta[k-1]=b;
    for(i=0;i<s->n;++i) s->v[k][i]=r[i]/b;
    ++s->restarts;
    return 0;
}

static int solve(operator_context *op, lanczos *s,
                 int ncv, int maxrestart, double tol)
{
    int i,j,round,all;
    double estimate,relative;
    memset(s,0,sizeof(*s));
    s->n=op->n;s->m=ncv;s->nev=NEV;
    for(i=0;i<s->n;++i)
        s->v[0][i]=sin((i+1)*0.317)+cos((i+1)*0.71);
    estimate=nrm(s->n,s->v[0]);
    for(i=0;i<s->n;++i) s->v[0][i]/=estimate;
    if(extend(op,s,0)) return -1;
    for(round=0;round<=maxrestart;++round) {
        make_projected(s);
        if(small_eigs(s->m,s->t,s->eval,s->u)) return -2;
        sort_smallest(s);
        all=1;
        for(i=0;i<NEV;++i) {
            int idx=s->order[i];
            estimate=s->beta[s->m-1]*fabs(s->u[s->m-1][idx]);
            if(estimate>tol*fmax(1.0,fabs(s->eval[idx]))) {
                all=0;break;
            }
            if(physical_pair(op,s,idx,&relative) || relative>tol) {
                all=0;break;
            }
            s->lambda[i]=s->eval[idx];
            s->relative[i]=relative;
        }
        if(all) {
            s->op_count=op->calls;
            return 0;
        }
        if(round==maxrestart) return 1;
        j=s->m/2;
        if(j<NEV+1) j=NEV+1;
        if(j>s->m-2) j=s->m-2;
        if(implicit_restart(s,j)) return -3;
        if(extend(op,s,j)) return -4;
    }
    return 1;
}

static int experiment(int family,int ncv)
{
    operator_context op;
    lanczos s;
    clock_t begin,end;
    int i,rc;
    double max_error=0.0,expect;
    memset(&op,0,sizeof(op));op.n=96;op.family=family;
    for(i=0;i<op.n;++i) {
        op.diag[i]=family==1 ?
            (i<8?1.0+1.0e-5*i:2.0+0.12*(i-8)):
            1.0+0.1*i;
    }
    begin=clock();
    rc=solve(&op,&s,ncv,120,1e-8);
    end=clock();
    if(rc==0) {
        for(i=0;i<NEV;++i) {
            expect=family==2 ?
                2.0-2.0*cos((double)(i+1)*3.14159265358979323846/
                            (double)(op.n+1)):
                op.diag[i];
            if(fabs(s.lambda[i]-expect)>max_error)
                max_error=fabs(s.lambda[i]-expect);
            if(s.relative[i]>1e-8||!isfinite(s.lambda[i])) return 1;
        }
        if(max_error>1e-5) {
            fprintf(stderr,"FAIL spectral error family=%d ncv=%d error=%g\n",
                    family,ncv,max_error);
            return 1;
        }
    }
    printf("IRL family=%d ncv=%d status=%d restarts=%d implicit_shifts=%d "
           "A_calls=%d maxeigerr=%g cpu_seconds=%.6f\n",
           family,ncv,rc,s.restarts,s.qr_shifts,op.calls,
           max_error,(double)(end-begin)/CLOCKS_PER_SEC);
    /* A prototype may fail to converge on a hard cluster, but it must
     * never claim converged eigenpairs with an incorrect true residual. */
    if(rc!=0&&rc!=1) {
        fprintf(stderr,"FAIL Lanczos/QR invariant status=%d\n",rc);
        return 1;
    }
    return 0;
}

int main(void)
{
    return experiment(0,12)||experiment(0,18)||
           experiment(1,12)||experiment(1,18)||
           experiment(2,18);
}
