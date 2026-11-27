/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

enum phase {
 IDLE, B_NORM, INITIAL, VERIFY, RES_NORM,
 CG_P, CG_RZ, CG_A, CG_PV, CG_NORM,
 BI_P, BI_RHO_FIRST, BI_A, BI_SV, BI_NORM_S, BI_S, BI_T,
 BI_TT, BI_TS, BI_NORM_R, BI_RHO_NEXT,
 GM_P, GM_A, GM_BEFORE, GM_DOT, GM_AFTER, END
};
struct lis_rci_context {
 lis_rci_options o; lis_rci_info info; enum phase phase; int verify_reason;
 double *mem,*b,*x,*r,*z,*p,*v,*s,*t,*shadow,*ax,*V,*Z,*H,*cs,*sn,*g,*y,*out;
 size_t j, orth_i; long double tol,rho,alpha,omega,work,before;
 int fresh, orth_pass, external;
 uint64_t token;
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
/* Request tokens are context-local and never reset on context reuse. */
static int end(lis_rci_context *c,int status,lis_rci_request_ex *q)
{
 c->phase=END; c->info.status=status; c->out=NULL;
 *q=(lis_rci_request_ex){0}; return status;
}
static int request(lis_rci_context *c,enum phase p,int status,
 const double *in,const double *in2,double *out,lis_rci_request_ex *q)
{
 if(!finite_vec(in,c->o.n)||(in2&&!finite_vec(in2,c->o.n)))
  return end(c,LIS_RCI_NONFINITE,q);
 if(c->token==UINT64_MAX)return end(c,LIS_RCI_PROTOCOL,q);
 if(out)for(size_t i=0;i<c->o.n;i++)out[i]=NAN;
 c->phase=p;c->out=out;c->info.status=status;
 *q=(lis_rci_request_ex){++c->token,status,c->o.n,in,in2,out};
 if(status==LIS_RCI_MATVEC)c->info.matvecs++;
 if(status==LIS_RCI_PRECON)c->info.precons++;
 return status;
}
static int req(lis_rci_context *c,enum phase p,int status,
 const double *in,double *out,lis_rci_request_ex *q)
{return request(c,p,status,in,NULL,out,q);}
static int scalar_req(lis_rci_context *c,enum phase p,const double *in,
 const double *in2,lis_rci_request_ex *q)
{return request(c,p,in2?LIS_RCI_DOT:LIS_RCI_NORM,in,in2,NULL,q);}
static int verify(lis_rci_context *c,int reason,lis_rci_request_ex *q)
{c->verify_reason=reason;return req(c,VERIFY,LIS_RCI_MATVEC,c->x,c->ax,q);}
lis_rci_context *lis_rci_create(const lis_rci_options*o){
 if(!o||!o->n||!o->maxiter||o->method<LIS_RCI_CG||o->method>LIS_RCI_FGMRES||!isfinite(o->rtol)||!isfinite(o->atol)||o->rtol<0||o->atol<0||(!o->rtol&&!o->atol))return NULL;
 size_t n=o->n,m=0,total,tmp,mp=0;if(!mul(n,10,&total))return NULL;
 if(o->method>=LIS_RCI_GMRES){if(!o->restart)return NULL;m=o->restart;if(m>o->maxiter)m=o->maxiter;if(m>n)m=n;if(!add(m,1,&mp)||!mul(mp,n,&tmp)||!add(total,tmp,&total)||!mul(m,n,&tmp)||!add(total,tmp,&total)||!mul(mp,m,&tmp)||!add(total,tmp,&total)||!mul(m,4,&tmp)||!add(total,tmp,&total)||!add(total,1,&total))return NULL;}
 if(!mul(total,sizeof(double),&tmp))return NULL;
 lis_rci_context*c=calloc(1,sizeof(*c));if(!c)return NULL;c->mem=calloc(total,sizeof(double));if(!c->mem){free(c);return NULL;}c->o=*o;c->o.restart=m;c->info.status=LIS_RCI_STATE;c->info.residual=NAN;
 double *d=c->mem;
 c->b=d;d+=n;c->x=d;d+=n;c->r=d;d+=n;c->z=d;d+=n;c->p=d;d+=n;c->v=d;d+=n;c->s=d;d+=n;c->t=d;d+=n;c->shadow=d;d+=n;c->ax=d;d+=n;
 if(m){c->V=d;d+=mp*n;c->Z=d;d+=m*n;c->H=d;d+=mp*m;c->cs=d;d+=m;c->sn=d;d+=m;c->g=d;d+=mp;c->y=d;}
 return c;
}
void lis_rci_destroy(lis_rci_context*c){if(c){free(c->mem);free(c);}}
const double*lis_rci_solution(const lis_rci_context*c){return c?c->x:NULL;}
int lis_rci_get_info(const lis_rci_context*c,lis_rci_info*i){if(!c||!i)return LIS_RCI_BADARG;*i=c->info;return 0;}
static int cycle(lis_rci_context *c,long double rn,lis_rci_request_ex *q)
{
 size_t n=c->o.n,m=c->o.restart;
 if(rn<=c->tol)return end(c,LIS_RCI_DONE,q);
 if(c->info.iterations>=c->o.maxiter)return end(c,LIS_RCI_MAXITER,q);
 c->fresh=1;c->rho=c->alpha=c->omega=1;
 if(c->o.method==LIS_RCI_CG)return req(c,CG_P,LIS_RCI_PRECON,c->r,c->z,q);
 if(c->o.method==LIS_RCI_BICGSTAB){
  copy(c->shadow,c->r,n);copy(c->p,c->r,n);
  return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);
 }
 if(rn>DBL_MAX)return end(c,LIS_RCI_NONFINITE,q);
 c->j=0;memset(c->H,0,(m+1)*m*sizeof(double));
 memset(c->g,0,(m+1)*sizeof(double));c->g[0]=(double)rn;
 for(size_t i=0;i<n;i++)c->V[i]=(double)((long double)c->r[i]/rn);
 return req(c,GM_P,LIS_RCI_PRECON,c->V,c->Z,q);
}
static int start_impl(lis_rci_context *c,const double *b,const double *x0,
 int external,lis_rci_request_ex *q)
{
 if(!c||!b||!q)return LIS_RCI_BADARG;
 if(c->phase!=IDLE&&c->phase!=END)return LIS_RCI_STATE;
 size_t n=c->o.n;
 c->info=(lis_rci_info){.residual=NAN};c->external=external;c->verify_reason=0;
 if(!finite_vec(b,n)||(x0&&!finite_vec(x0,n)))return end(c,LIS_RCI_NONFINITE,q);
 memmove(c->b,b,n*sizeof(double));
 if(x0)memmove(c->x,x0,n*sizeof(double));else memset(c->x,0,n*sizeof(double));
 return scalar_req(c,B_NORM,c->b,NULL,q);
}
static int gm_update(lis_rci_context*c,size_t k,int reason,lis_rci_request_ex*q){size_t n=c->o.n,m=c->o.restart;
 for(size_t ii=k;ii>0;ii--){size_t i=ii-1;long double a=c->g[i];for(size_t j=i+1;j<k;j++)a-=(long double)c->H[j*(m+1)+i]*c->y[j];double d=c->H[i*(m+1)+i];if(d==0||!isfinite(d))return verify(c,LIS_RCI_BREAKDOWN,q);c->y[i]=(double)(a/d);if(!isfinite(c->y[i]))return verify(c,LIS_RCI_BREAKDOWN,q);}
 for(size_t j=0;j<k;j++)axpy(c->x,c->Z+j*n,c->y[j],n);return verify(c,reason,q);}
static int advance_impl(lis_rci_context *c,uint64_t token,int code,
 long double scalar,lis_rci_request_ex *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 /* A stale or duplicate completion must leave the outstanding request intact. */
 if(token!=c->token)return LIS_RCI_PROTOCOL;
 if(code){c->info.callback_code=code;return end(c,LIS_RCI_CALLBACK,q);}
 if(c->info.status==LIS_RCI_DOT||c->info.status==LIS_RCI_NORM){
  if(!isfinite(scalar)||(c->info.status==LIS_RCI_NORM&&scalar<0))
   return end(c,LIS_RCI_NONFINITE,q);
 }else if(!finite_vec(c->out,c->o.n))return end(c,LIS_RCI_NONFINITE,q);
 size_t n=c->o.n,m=c->o.restart;long double a,b,rn;
 switch(c->phase){
 case B_NORM:
  c->tol=fmaxl(c->o.atol,(long double)c->o.rtol*scalar);
  if(!isfinite(c->tol))return end(c,LIS_RCI_NONFINITE,q);
  return req(c,INITIAL,LIS_RCI_MATVEC,c->x,c->ax,q);
 case INITIAL:case VERIFY:
  for(size_t i=0;i<n;i++)c->r[i]=c->b[i]-c->ax[i];
  return scalar_req(c,RES_NORM,c->r,NULL,q);
 case RES_NORM:
  c->info.residual=(double)scalar;
  if(scalar<=c->tol)return end(c,LIS_RCI_DONE,q);
  if(c->verify_reason==LIS_RCI_BREAKDOWN)return end(c,LIS_RCI_BREAKDOWN,q);
  return cycle(c,scalar,q);
 case CG_P:return scalar_req(c,CG_RZ,c->r,c->z,q);
 case CG_RZ:
  a=scalar;if(a<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
  if(c->fresh){copy(c->p,c->z,n);c->fresh=0;}
  else{b=a/c->rho;for(size_t i=0;i<n;i++)c->p[i]=(double)((long double)c->z[i]+b*c->p[i]);}
  c->rho=a;return req(c,CG_A,LIS_RCI_MATVEC,c->p,c->v,q);
 case CG_A:return scalar_req(c,CG_PV,c->p,c->v,q);
 case CG_PV:
  if(scalar<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
  b=c->rho/scalar;if(!isfinite(b))return verify(c,LIS_RCI_BREAKDOWN,q);
  axpy(c->x,c->p,b,n);axpy(c->r,c->v,-b,n);c->info.iterations++;
  if(!finite_vec(c->x,n))return end(c,LIS_RCI_NONFINITE,q);
  return scalar_req(c,CG_NORM,c->r,NULL,q);
 case CG_NORM:
  if(scalar<=c->tol||c->info.iterations>=c->o.maxiter)return verify(c,LIS_RCI_MAXITER,q);
  return req(c,CG_P,LIS_RCI_PRECON,c->r,c->z,q);
 case BI_P:
  if(c->fresh)return scalar_req(c,BI_RHO_FIRST,c->shadow,c->r,q);
  return req(c,BI_A,LIS_RCI_MATVEC,c->z,c->v,q);
 case BI_RHO_FIRST:
  c->rho=scalar;c->fresh=0;if(!scalar)return verify(c,LIS_RCI_BREAKDOWN,q);
  return req(c,BI_A,LIS_RCI_MATVEC,c->z,c->v,q);
 case BI_A:return scalar_req(c,BI_SV,c->shadow,c->v,q);
 case BI_SV:
  if(!scalar)return verify(c,LIS_RCI_BREAKDOWN,q);
  c->alpha=c->rho/scalar;if(!isfinite(c->alpha))return verify(c,LIS_RCI_BREAKDOWN,q);
  copy(c->s,c->r,n);axpy(c->s,c->v,-c->alpha,n);
  return scalar_req(c,BI_NORM_S,c->s,NULL,q);
 case BI_NORM_S:
  if(scalar<=c->tol){axpy(c->x,c->z,c->alpha,n);c->info.iterations++;return verify(c,LIS_RCI_MAXITER,q);}
  return req(c,BI_S,LIS_RCI_PRECON,c->s,c->ax,q);
 case BI_S:return req(c,BI_T,LIS_RCI_MATVEC,c->ax,c->t,q);
 case BI_T:return scalar_req(c,BI_TT,c->t,c->t,q);
 case BI_TT:
  if(scalar<=0)return verify(c,LIS_RCI_BREAKDOWN,q);
  c->work=scalar;return scalar_req(c,BI_TS,c->t,c->s,q);
 case BI_TS:
  c->omega=scalar/c->work;if(!c->omega||!isfinite(c->omega))return verify(c,LIS_RCI_BREAKDOWN,q);
  axpy(c->x,c->z,c->alpha,n);axpy(c->x,c->ax,c->omega,n);
  copy(c->r,c->s,n);axpy(c->r,c->t,-c->omega,n);c->info.iterations++;
  if(!finite_vec(c->x,n))return end(c,LIS_RCI_NONFINITE,q);
  return scalar_req(c,BI_NORM_R,c->r,NULL,q);
 case BI_NORM_R:
  if(scalar<=c->tol||c->info.iterations>=c->o.maxiter)return verify(c,LIS_RCI_MAXITER,q);
  return scalar_req(c,BI_RHO_NEXT,c->shadow,c->r,q);
 case BI_RHO_NEXT:
  if(!scalar||!c->omega)return verify(c,LIS_RCI_BREAKDOWN,q);
  b=(scalar/c->rho)*(c->alpha/c->omega);
  if(!isfinite(b))return verify(c,LIS_RCI_BREAKDOWN,q);
  for(size_t i=0;i<n;i++)c->p[i]=(double)((long double)c->r[i]+b*((long double)c->p[i]-c->omega*c->v[i]));
  c->rho=scalar;return req(c,BI_P,LIS_RCI_PRECON,c->p,c->z,q);
 case GM_P:return req(c,GM_A,LIS_RCI_MATVEC,c->Z+c->j*n,c->v,q);
 case GM_A:return scalar_req(c,GM_BEFORE,c->v,NULL,q);
 case GM_BEFORE:
  c->before=scalar;c->orth_i=0;c->orth_pass=0;
  return scalar_req(c,GM_DOT,c->V,c->v,q);
 case GM_DOT:{
  double *h=c->H+c->j*(m+1);size_t i=c->orth_i;
  h[i]+=(double)scalar;axpy(c->v,c->V+i*n,-scalar,n);
  if(++c->orth_i>c->j){c->orth_i=0;c->orth_pass++;}
  if(c->orth_pass<2)return scalar_req(c,GM_DOT,c->V+c->orth_i*n,c->v,q);
  return scalar_req(c,GM_AFTER,c->v,NULL,q);
 }
 case GM_AFTER:{
  size_t j=c->j;double *h=c->H+j*(m+1);rn=scalar;
  if(rn>DBL_MAX)return end(c,LIS_RCI_NONFINITE,q);
  h[j+1]=(double)rn;int happy=rn<=DBL_EPSILON*c->before;
  if(!happy)for(size_t i=0;i<n;i++)c->V[(j+1)*n+i]=(double)((long double)c->v[i]/rn);
  for(size_t i=0;i<j;i++){double u=h[i],w=h[i+1];h[i]=c->cs[i]*u+c->sn[i]*w;h[i+1]=-c->sn[i]*u+c->cs[i]*w;}
  double d=hypot(h[j],h[j+1]);
  if(!isfinite(d)||d==0){if(j)return gm_update(c,j,LIS_RCI_BREAKDOWN,q);return verify(c,LIS_RCI_BREAKDOWN,q);}
  c->cs[j]=h[j]/d;c->sn[j]=h[j+1]/d;h[j]=d;h[j+1]=0;
  c->g[j+1]=-c->sn[j]*c->g[j];c->g[j]*=c->cs[j];c->info.iterations++;
  if(happy||fabs(c->g[j+1])<=c->tol||j+1==m||c->info.iterations>=c->o.maxiter)
   return gm_update(c,j+1,happy?LIS_RCI_BREAKDOWN:LIS_RCI_MAXITER,q);
  c->j++;return req(c,GM_P,LIS_RCI_PRECON,c->V+c->j*n,c->Z+c->j*n,q);
 }
 default:return end(c,LIS_RCI_STATE,q);
 }
}
int lis_rci_start_ex(lis_rci_context *c,const double *b,const double *x0,lis_rci_request_ex *q)
{return start_impl(c,b,x0,1,q);}
int lis_rci_advance_ex(lis_rci_context *c,uint64_t token,int code,long double scalar,lis_rci_request_ex *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(!c->external)return LIS_RCI_STATE;
 return advance_impl(c,token,code,scalar,q);
}
/* Compatibility entry points complete scalar reductions with the local backend. */
static int local_continue(lis_rci_context *c,int status,lis_rci_request_ex *e,lis_rci_request *q)
{
 while(status==LIS_RCI_DOT||status==LIS_RCI_NORM){
  long double s=status==LIS_RCI_DOT?dot(e->in,e->in2,e->n):norm(e->in,e->n);
  status=advance_impl(c,e->token,0,s,e);
 }
 *q=(lis_rci_request){e->in,e->out,e->n};return status;
}
int lis_rci_start(lis_rci_context *c,const double *b,const double *x0,lis_rci_request *q)
{
 if(!q)return LIS_RCI_BADARG;
 lis_rci_request_ex e={0};int s=start_impl(c,b,x0,0,&e);
 if(s==LIS_RCI_BADARG||s==LIS_RCI_STATE)return s;
 return local_continue(c,s,&e,q);
}
int lis_rci_advance(lis_rci_context *c,int code,lis_rci_request *q)
{
 if(!c||!q)return LIS_RCI_BADARG;
 if(c->external||c->phase==IDLE||c->phase==END)return LIS_RCI_STATE;
 lis_rci_request_ex e={0};int s=advance_impl(c,c->token,code,0,&e);
 return local_continue(c,s,&e,q);
}
/* Transport-independent five-word maximum-reduction protocol. */
int lis_rci_consensus_pack(uint64_t token,uint64_t tag,uint32_t rank,int32_t error,uint64_t words[5])
{
 if(!words)return LIS_RCI_BADARG;
 words[0]=token;words[1]=~token;words[2]=tag;words[3]=~tag;
 words[4]=error?((uint64_t)(UINT32_MAX-rank)<<32)|(uint32_t)error:0;
 return LIS_RCI_DONE;
}
int lis_rci_consensus_merge(uint64_t dst[5],const uint64_t src[5])
{
 if(!dst||!src)return LIS_RCI_BADARG;
 for(size_t i=0;i<5;i++)if(src[i]>dst[i])dst[i]=src[i];
 return LIS_RCI_DONE;
}
int lis_rci_consensus_unpack(const uint64_t words[5],lis_rci_consensus_result *result)
{
 if(!words||!result||~words[1]>words[0]||~words[3]>words[2])return LIS_RCI_BADARG;
 uint32_t code=(uint32_t)words[4];
 if(words[4]&&!code)return LIS_RCI_BADARG;
 *result=(lis_rci_consensus_result){
  words[0]!=~words[1],words[2]!=~words[3],
  words[4]?UINT32_MAX-(uint32_t)(words[4]>>32):UINT32_MAX,
  (int32_t)(code<=INT32_MAX?(int64_t)code:(int64_t)code-INT64_C(4294967296))
 };
 return LIS_RCI_DONE;
}
