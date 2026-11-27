/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_krylov_restore.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define N 24
#define NEV 3
#define C(e) do{if(!(e)){fprintf(stderr,"KRYLOV FAIL %d: %s\n",__LINE__,#e);return 1;}}while(0)
typedef struct {uint64_t hash,count;int mode,psd;double sigma;uint64_t fail;} trace;
static void mix(trace *t,const void *p,size_t n)
{const unsigned char *a=p;for(size_t i=0;i<n;i++){t->hash^=a[i];t->hash*=UINT64_C(1099511628211);}}
static int apply(trace *t,uint32_t op,int64_t n,const double *x,double *y)
{
 t->count++;mix(t,&op,sizeof op);mix(t,x,(size_t)n*sizeof(double));
 if(t->count==t->fail)return 71;
 for(int64_t i=0;i<n;i++){
  double B=t->psd&&i>=n-4?0:1+.02*i;
  double A=(1+.11*i)*(1+.02*i);
  y[i]=op==LIS_RCI_OP_A?A*x[i]:op==LIS_RCI_OP_B?B*x[i]:x[i]/(A-t->sigma*B);
 }
 mix(t,y,(size_t)n*sizeof(double));return 0;
}
static int a(void *u,int64_t n,const double *x,double *y){return apply(u,LIS_RCI_OP_A,n,x,y);}
static int b(void *u,int64_t n,const double *x,double *y){return apply(u,LIS_RCI_OP_B,n,x,y);}
static int solve(void *u,int64_t n,const double *x,double *y){return apply(u,LIS_RCI_OP_SHIFTED_SOLVE,n,x,y);}
static int compare(int mode,int psd,int locking,int adaptive,const char *which,int reuse_flag,uint64_t fail)
{
 trace t1={.hash=UINT64_C(1469598103934665603),.mode=mode,.psd=psd,.sigma=.7,.fail=fail},t2=t1;
 lis_krylov_problem64 p={N,mode,t1.sigma,&t1,a,b,solve},p2=p;p2.ctx=&t2;
 lis_krylov_options o;lis_krylov_options_default(&o);
 o.nev=NEV;o.ncv=9;o.maxiter=120;o.tol=1e-9;o.which=which;o.hard_locking=locking;
 o.allow_semidefinite_metric=psd;o.inner_tol=mode?1e-12:0;
 lis_krylov_ext_options ext;lis_krylov_ext_options_default(&ext);
 ext.adaptive_ncv=adaptive;ext.ncv_max=adaptive?16:0;ext.restart_directions=adaptive?2:1;
 lis_krylov_workspace *w1=NULL,*w2=NULL;
 if(reuse_flag){C(!lis_krylov_workspace_create64(&p,&o,&ext,&w1));C(!lis_krylov_workspace_create64(&p2,&o,&ext,&w2));}
 double v1[NEV]={0},v2[NEV]={0},x1[NEV*N]={0},x2[NEV*N]={0},r1[NEV]={0},r2[NEV]={0};
 lis_krylov_stats s1,s2;lis_krylov_ext_stats e1,e2;lis_krylov_ext_stats_default(&e1);lis_krylov_ext_stats_default(&e2);
 int rc1=lis_krylov_symmetric_ex64(&p,&o,&ext,w1,v1,x1,r1,&s1,&e1);
 lis_rci_krylov_restore *K=lis_rci_krylov_restore_create();C(K);
 lis_rci_krylov_restore_request q;
 int st=lis_rci_krylov_restore_start(K,&p2,&o,&ext,w2,NULL,v2,x2,r2,&s2,&e2,&q);
 uint64_t last=0;
 while(st>0){
  C(q.token>last&&q.n==N&&q.operation==(uint32_t)st&&q.sigma==p.sigma);last=q.token;
  C(lis_rci_krylov_restore_advance(K,q.token-1,0,&q)==LIS_RCI_PROTOCOL&&q.token==last);
  int err=apply(&t2,q.operation,q.n,q.in,q.out);
  st=lis_rci_krylov_restore_advance(K,q.token,err,&q);
 }
 int rc2=-1;lis_krylov_budget_stats counts;
 C(lis_rci_krylov_restore_status(K,&rc2,&counts)==st);
 if(rc1!=rc2||t1.hash!=t2.hash||t1.count!=t2.count)
  fprintf(stderr,"DIFF mode=%d psd=%d locking=%d adaptive=%d which=%s rc=%d/%d count=%llu/%llu hash=%llx/%llx\n",mode,psd,locking,adaptive,which,rc1,rc2,(unsigned long long)t1.count,(unsigned long long)t2.count,(unsigned long long)t1.hash,(unsigned long long)t2.hash);
 C(rc1==rc2&&t1.hash==t2.hash&&t1.count==t2.count&&counts.total==t2.count);
 C(s1.nconv==s2.nconv&&s1.nreturned==s2.nreturned&&s1.iterations==s2.iterations&&s1.restarts==s2.restarts&&s1.op_count==s2.op_count&&s1.metric_count==s2.metric_count&&s1.solve_count==s2.solve_count&&s1.max_residual==s2.max_residual);
 C(e1.ncv_initial==e2.ncv_initial&&e1.ncv_final==e2.ncv_final&&e1.ncv_growths==e2.ncv_growths&&e1.basis_vectors==e2.basis_vectors&&e1.projected_solves==e2.projected_solves&&e1.extra_seed_attempts==e2.extra_seed_attempts&&e1.extra_seed_accepted==e2.extra_seed_accepted);
 C(!memcmp(v1,v2,sizeof v1)&&!memcmp(x1,x2,sizeof x1)&&!memcmp(r1,r2,sizeof r1));
 lis_rci_krylov_restore_destroy(K);lis_krylov_workspace_destroy(w1);lis_krylov_workspace_destroy(w2);return 0;
}
static int cancellation(void)
{
 trace t={.hash=1,.sigma=.7};lis_krylov_problem64 p={N,1,t.sigma,&t,a,b,solve};
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=NEV;o.ncv=9;o.inner_tol=1e-12;
 lis_krylov_workspace *w=NULL;C(!lis_krylov_workspace_create64(&p,&o,NULL,&w));
 lis_rci_krylov_restore *K=lis_rci_krylov_restore_create();C(K);
 double v[NEV],x[NEV*N],r[NEV];lis_krylov_stats s;
 lis_rci_krylov_restore_request q;uint64_t last=0;
 for(int stop=0;stop<20;stop++){
  int st=lis_rci_krylov_restore_start(K,&p,&o,NULL,w,NULL,v,x,r,&s,NULL,&q);
  C(st>0&&q.token>last);
  for(int i=0;i<stop;i++){C(st>0);int err=apply(&t,q.operation,q.n,q.in,q.out);st=lis_rci_krylov_restore_advance(K,q.token,err,&q);}
  C(st>0);last=q.token;
  C(lis_rci_krylov_restore_cancel(K,q.token-1)==LIS_RCI_PROTOCOL);
  C(lis_rci_krylov_restore_cancel(K,q.token)==LIS_RCI_CANCELLED&&s.nreturned==0);
 }
 lis_krylov_budget budget;lis_krylov_budget_default(&budget);budget.max_total=7;
 int st=lis_rci_krylov_restore_start(K,&p,&o,NULL,w,&budget,v,x,r,&s,NULL,&q);
 while(st>0){int err=apply(&t,q.operation,q.n,q.in,q.out);st=lis_rci_krylov_restore_advance(K,q.token,err,&q);}
 int native;lis_krylov_budget_stats counts;
 C(st==LIS_RCI_MAXITER&&lis_rci_krylov_restore_status(K,&native,&counts)==st);
 C(native==LIS_KRYLOV_BUDGET_EXHAUSTED&&counts.total==7&&counts.exhausted&&s.nreturned==0);
 st=lis_rci_krylov_restore_start(K,&p,&o,NULL,w,NULL,v,x,r,&s,NULL,&q);
 while(st>0){int err=apply(&t,q.operation,q.n,q.in,q.out);st=lis_rci_krylov_restore_advance(K,q.token,err,&q);}
 C(st==0||st==LIS_RCI_MAXITER);
 lis_rci_krylov_restore_destroy(K);lis_krylov_workspace_destroy(w);return 0;
}
int main(void)
{
 const char *which[]={"LM","SM","LA:ASC","SA:DESC","BE"};
 unsigned cases=0;
 for(int mode=0;mode<3;mode++)for(int lock=0;lock<2;lock++)for(int adapt=0;adapt<2;adapt++)for(int j=0;j<5;j++){
  C(!compare(mode,0,lock,adapt,which[j],(lock+adapt)%2,0));cases++;
 }
 for(int lock=0;lock<2;lock++)for(int adapt=0;adapt<2;adapt++){C(!compare(1,1,lock,adapt,"LM",1,0));cases++;}
 for(uint64_t f=1;f<=12;f++){C(!compare(1,0,0,1,"LM",1,f));cases++;}
 C(!cancellation());
 printf("KRYLOV RESTORE PASS: %u bitwise native differential cases, 20 cancellation boundaries, budget and workspace reuse\n",cases);
 return 0;
}
