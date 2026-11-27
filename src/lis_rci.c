#include "lis_rci.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

enum phase { IDLE, INITIAL, CG_P, CG_A, BI_P, BI_A, BI_S, BI_T, GM_P, GM_A, VERIFY, END };
struct lis_rci_context {
 lis_rci_options o; lis_rci_info info; enum phase phase; int verify_reason;
 double *mem,*b,*x,*r,*z,*p,*v,*s,*t,*shadow,*ax,*V,*Z,*H,*cs,*sn,*g,*y,*out;
 size_t j; long double tol,rho,alpha,omega; int fresh;
};
static int mul(size_t a,size_t b,size_t *z){if(b && a>SIZE_MAX/b)return 0;*z=a*b;return 1;}
static int add(size_t a,size_t b,size_t *z){if(a>SIZE_MAX-b)return 0;*z=a+b;return 1;}
static int finite_vec(const double *a,size_t n){int ok=1;
#ifdef _OPENMP
#pragma omp parallel for reduction(&:ok) if(n>=4096)
#endif
 for(size_t i=0;i<n;i++)ok &= isfinite(a[i])!=0;return ok;}
static long double dot(const double *a,const double *b,size_t n){long double s=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s) if(n>=4096)
#endif
 for(size_t i=0;i<n;i++)s+=(long double)a[i]*b[i];return s;}
static long double norm(const double *a,size_t n){/* scaled to avoid overflow/underflow even with narrow long double */
 double m=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(max:m) if(n>=4096)
#endif
 for(size_t i=0;i<n;i++)if(fabs(a[i])>m)m=fabs(a[i]);
 if(m==0)return 0;long double s=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s) if(n>=4096)
#endif
 for(size_t i=0;i<n;i++){long double v=(long double)a[i]/m;s+=v*v;}return (long double)m*sqrtl(s);}
static void copy(double *a,const double *b,size_t n){memcpy(a,b,n*sizeof(double));}
static void axpy(double *a,const double *b,long double k,size_t n){
#ifdef _OPENMP
#pragma omp parallel for if(n>=4096)
#endif
 for(size_t i=0;i<n;i++)a[i]=(double)((long double)a[i]+k*b[i]);}
static int end(lis_rci_context*c,int status,lis_rci_request*q){c->phase=END;c->info.status=status;if(q)*q=(lis_rci_request){0};return status;}
static int req(lis_rci_context*c,enum phase p,int status,const double*in,double*out,lis_rci_request*q){
 if(!finite_vec(in,c->o.n))return end(c,LIS_RCI_NONFINITE,q);
 for(size_t i=0;i<c->o.n;i++)out[i]=NAN; /* detect incomplete callbacks */
 c->phase=p;c->out=out;c->info.status=status;*q=(lis_rci_request){in,out,c->o.n};if(status==LIS_RCI_MATVEC)c->info.matvecs++;else c->info.precons++;return status;}
static int verify(lis_rci_context*c,int reason,lis_rci_request*q){c->verify_reason=reason;return req(c,VERIFY,LIS_RCI_MATVEC,c->x,c->ax,q);}
lis_rci_context *lis_rci_create(const lis_rci_options*o){
 if(!o||!o->n||!o->maxiter||o->method<LIS_RCI_CG||o->method>LIS_RCI_FGMRES||!isfinite(o->rtol)||!isfinite(o->atol)||o->rtol<0||o->atol<0||(!o->rtol&&!o->atol))return NULL;
 size_t n=o->n,m=0,total,tmp,mp=0;if(!mul(n,10,&total))return NULL;
 if(o->method>=LIS_RCI_GMRES){if(!o->restart)return NULL;m=o->restart;if(m>o->maxiter)m=o->maxiter;if(m>n)m=n;if(!add(m,1,&mp)||!mul(mp,n,&tmp)||!add(total,tmp,&total)||!mul(m,n,&tmp)||!add(total,tmp,&total)||!mul(mp,m,&tmp)||!add(total,tmp,&total)||!mul(m,4,&tmp)||!add(total,tmp,&total)||!add(total,1,&total))return NULL;}
 if(!mul(total,sizeof(double),&tmp))return NULL;
 lis_rci_context*c=calloc(1,sizeof(*c));if(!c)return NULL;c->mem=calloc(total,sizeof(double));if(!c->mem){free(c);return NULL;}c->o=*o;c->o.restart=m;
 double *d=c->mem;
 c->b=d;d+=n;c->x=d;d+=n;c->r=d;d+=n;c->z=d;d+=n;c->p=d;d+=n;c->v=d;d+=n;c->s=d;d+=n;c->t=d;d+=n;c->shadow=d;d+=n;c->ax=d;d+=n;
 if(m){c->V=d;d+=mp*n;c->Z=d;d+=m*n;c->H=d;d+=mp*m;c->cs=d;d+=m;c->sn=d;d+=m;c->g=d;d+=mp;c->y=d;}
 return c;
}
void lis_rci_destroy(lis_rci_context*c){if(c){free(c->mem);free(c);}}
const double*lis_rci_solution(const lis_rci_context*c){return c?c->x:NULL;}
int lis_rci_get_info(const lis_rci_context*c,lis_rci_info*i){if(!c||!i)return LIS_RCI_BADARG;*i=c->info;return 0;}
static int cycle(lis_rci_context*c,lis_rci_request*q){size_t n=c->o.n,m=c->o.restart;long double rn=norm(c->r,n);c->info.residual=(double)rn;
 if(!isfinite(rn))return end(c,LIS_RCI_NONFINITE,q);if(rn<=c->tol)return end(c,LIS_RCI_DONE,q);if(c->info.iterations>=c->o.maxiter)return end(c,LIS_RCI_MAXITER,q);
 c->fresh=1;c->rho=c->alpha=c->omega=1;
 if(c->o.method==LIS_RCI_CG)return req(c,CG_P,LIS_RCI_PRECON,c->r,c->z,q);
 if(c->o.method==LIS_RCI_BICGSTAB){copy(c->shadow,c->r,n);copy(c->p,c->r,n);return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);}
 if(rn>DBL_MAX)return end(c,LIS_RCI_NONFINITE,q);
 c->j=0;memset(c->H,0,(m+1)*m*sizeof(double));memset(c->g,0,(m+1)*sizeof(double));c->g[0]=(double)rn;
 for(size_t i=0;i<n;i++)c->V[i]=(double)((long double)c->r[i]/rn);
 return req(c,GM_P,LIS_RCI_PRECON,c->V,c->Z,q);
}
int lis_rci_start(lis_rci_context*c,const double*b,const double*x0,lis_rci_request*q){if(!c||!b||!q)return LIS_RCI_BADARG;if(c->phase!=IDLE&&c->phase!=END)return LIS_RCI_STATE;
 size_t n=c->o.n; c->info=(lis_rci_info){.residual=NAN};if(!finite_vec(b,n)||(x0&&!finite_vec(x0,n)))return end(c,LIS_RCI_NONFINITE,q);
 memmove(c->b,b,n*sizeof(double));if(x0)memmove(c->x,x0,n*sizeof(double));else memset(c->x,0,n*sizeof(double));
 c->tol=fmaxl(c->o.atol,(long double)c->o.rtol*norm(c->b,n));if(!isfinite(c->tol))return end(c,LIS_RCI_NONFINITE,q);
 return req(c,INITIAL,LIS_RCI_MATVEC,c->x,c->ax,q);}
static int bi_next(lis_rci_context*c,lis_rci_request*q){size_t n=c->o.n;long double rho=dot(c->shadow,c->r,n);if(!isfinite(rho)||rho==0||c->omega==0)return verify(c,LIS_RCI_BREAKDOWN,q);
 long double beta=(rho/c->rho)*(c->alpha/c->omega);if(!isfinite(beta))return verify(c,LIS_RCI_BREAKDOWN,q);
 for(size_t i=0;i<n;i++)c->p[i]=(double)((long double)c->r[i]+beta*((long double)c->p[i]-c->omega*c->v[i]));c->rho=rho;
 return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);}
static int gm_update(lis_rci_context*c,size_t k,int reason,lis_rci_request*q){size_t n=c->o.n,m=c->o.restart;
 for(size_t ii=k;ii>0;ii--){size_t i=ii-1;long double a=c->g[i];for(size_t j=i+1;j<k;j++)a-=(long double)c->H[j*(m+1)+i]*c->y[j];double d=c->H[i*(m+1)+i];if(d==0||!isfinite(d))return verify(c,LIS_RCI_BREAKDOWN,q);c->y[i]=(double)(a/d);if(!isfinite(c->y[i]))return verify(c,LIS_RCI_BREAKDOWN,q);}
 for(size_t j=0;j<k;j++)axpy(c->x,c->Z+j*n,c->y[j],n);return verify(c,reason,q);}
int lis_rci_advance(lis_rci_context*c,int code,lis_rci_request*q){if(!c||!q)return LIS_RCI_BADARG;if(c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 size_t n=c->o.n,m=c->o.restart;long double a,b,rn; if(code){c->info.callback_code=code;return end(c,LIS_RCI_CALLBACK,q);}if(!finite_vec(c->out,n))return end(c,LIS_RCI_NONFINITE,q);
 switch(c->phase){
 case INITIAL:case VERIFY:
 for(size_t i=0;i<n;i++)c->r[i]=c->b[i]-c->ax[i];if(!finite_vec(c->r,n))return end(c,LIS_RCI_NONFINITE,q);
 rn=norm(c->r,n);c->info.residual=(double)rn;if(rn<=c->tol)return end(c,LIS_RCI_DONE,q);if(c->phase==VERIFY&&c->verify_reason==LIS_RCI_BREAKDOWN)return end(c,LIS_RCI_BREAKDOWN,q);return cycle(c,q);
 case CG_P:
 a=dot(c->r,c->z,n);if(!isfinite(a)||a<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
 if(c->fresh){copy(c->p,c->z,n);c->fresh=0;}else {b=a/c->rho;for(size_t i=0;i<n;i++)c->p[i]=(double)((long double)c->z[i]+b*c->p[i]);}c->rho=a;
 return req(c,CG_A,LIS_RCI_MATVEC,c->p,c->v,q);
 case CG_A:
 a=dot(c->p,c->v,n);if(!isfinite(a)||a<=0)return verify(c,LIS_RCI_BREAKDOWN,q);b=c->rho/a;if(!isfinite(b))return verify(c,LIS_RCI_BREAKDOWN,q);
 axpy(c->x,c->p,b,n);axpy(c->r,c->v,-b,n);c->info.iterations++;
 if(!finite_vec(c->r,n)||!finite_vec(c->x,n))return end(c,LIS_RCI_NONFINITE,q);
 if(norm(c->r,n)<=c->tol||c->info.iterations>=c->o.maxiter)return verify(c,LIS_RCI_MAXITER,q);
 return req(c,CG_P,LIS_RCI_PRECON,c->r,c->z,q);
 case BI_P:
 if(c->fresh){c->rho=dot(c->shadow,c->r,n);c->fresh=0;if(c->rho==0||!isfinite(c->rho))return verify(c,LIS_RCI_BREAKDOWN,q);}return req(c,BI_A,LIS_RCI_MATVEC,c->z,c->v,q);
 case BI_A:
 a=dot(c->shadow,c->v,n);if(a==0||!isfinite(a))return verify(c,LIS_RCI_BREAKDOWN,q);c->alpha=c->rho/a;if(!isfinite(c->alpha))return verify(c,LIS_RCI_BREAKDOWN,q);
 copy(c->s,c->r,n);axpy(c->s,c->v,-c->alpha,n);if(!finite_vec(c->s,n))return end(c,LIS_RCI_NONFINITE,q);
 if(norm(c->s,n)<=c->tol){axpy(c->x,c->z,c->alpha,n);c->info.iterations++;return verify(c,LIS_RCI_MAXITER,q);}return req(c,BI_S,LIS_RCI_PRECON,c->s,c->ax,q);
 case BI_S: return req(c,BI_T,LIS_RCI_MATVEC,c->ax,c->t,q);
 case BI_T:
 a=dot(c->t,c->t,n);b=dot(c->t,c->s,n);if(a==0||!isfinite(a)||!isfinite(b))return verify(c,LIS_RCI_BREAKDOWN,q);c->omega=b/a;if(c->omega==0||!isfinite(c->omega))return verify(c,LIS_RCI_BREAKDOWN,q);
 axpy(c->x,c->z,c->alpha,n);axpy(c->x,c->ax,c->omega,n);copy(c->r,c->s,n);axpy(c->r,c->t,-c->omega,n);c->info.iterations++;
 if(!finite_vec(c->r,n)||!finite_vec(c->x,n))return end(c,LIS_RCI_NONFINITE,q);
 if(norm(c->r,n)<=c->tol||c->info.iterations>=c->o.maxiter)return verify(c,LIS_RCI_MAXITER,q);return bi_next(c,q);
 case GM_P:return req(c,GM_A,LIS_RCI_MATVEC,c->Z+c->j*n,c->v,q);
 case GM_A:{size_t j=c->j;double*h=c->H+j*(m+1);long double before=norm(c->v,n);
 /* Twice-modified Gram-Schmidt. */
 for(int pass=0;pass<2;pass++)for(size_t i=0;i<=j;i++){a=dot(c->V+i*n,c->v,n);h[i]+=(double)a;axpy(c->v,c->V+i*n,-a,n);}
 if(!finite_vec(c->v,n))return end(c,LIS_RCI_NONFINITE,q);rn=norm(c->v,n);if(!isfinite(rn)||rn>DBL_MAX)return end(c,LIS_RCI_NONFINITE,q);h[j+1]=(double)rn;
 int happy=rn<=DBL_EPSILON*before;
 if(!happy)for(size_t i=0;i<n;i++)c->V[(j+1)*n+i]=(double)((long double)c->v[i]/rn);
 for(size_t i=0;i<j;i++){double u=h[i],w=h[i+1];h[i]=c->cs[i]*u+c->sn[i]*w;h[i+1]=-c->sn[i]*u+c->cs[i]*w;}
 double d=hypot(h[j],h[j+1]);if(!isfinite(d)||d==0){if(j)return gm_update(c,j,LIS_RCI_BREAKDOWN,q);return verify(c,LIS_RCI_BREAKDOWN,q);}
 c->cs[j]=h[j]/d;c->sn[j]=h[j+1]/d;h[j]=d;h[j+1]=0;c->g[j+1]=-c->sn[j]*c->g[j];c->g[j]*=c->cs[j];c->info.iterations++;
 if(happy||fabs(c->g[j+1])<=c->tol||j+1==m||c->info.iterations>=c->o.maxiter)return gm_update(c,j+1,happy?LIS_RCI_BREAKDOWN:LIS_RCI_MAXITER,q);
 c->j++;return req(c,GM_P,LIS_RCI_PRECON,c->V+c->j*n,c->Z+c->j*n,q);}
 default:return end(c,LIS_RCI_STATE,q);
 }
}
