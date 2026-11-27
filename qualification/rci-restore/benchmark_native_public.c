/* Baseline-client ABI and public classic dispatch benchmark. */
#define _POSIX_C_SOURCE 200809L
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"BENCH FAIL %d: %s\n",__LINE__,#x);exit(1);}}while(0)
static double now(void){struct timespec t;CHECK(!clock_gettime(CLOCK_MONOTONIC,&t));return t.tv_sec+t.tv_nsec*1e-9;}
static uint64_t hash;
static void mix(const void *p,size_t n){const unsigned char *s=p;for(size_t i=0;i<n;i++){hash^=s[i];hash*=UINT64_C(1099511628211);}}
#define MIX(x) mix(&(x),sizeof(x))
static void init_hash(void){hash=UINT64_C(1469598103934665603);}
static void line(const char *s,LIS_INT n,int repeat,double elapsed)
{printf("BENCH %s %lld %d %016llx %.9f\n",s,(long long)n,repeat,(unsigned long long)hash,elapsed/repeat);}
static LIS_MATRIX matrix(LIS_INT n,int diagonal)
{
 LIS_MATRIX A=NULL;CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,n));
 const double eig[4]={1,2,4,7};
 for(LIS_INT i=0;i<n;i++){
  CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,diagonal?eig[i%4]:2.6+.003*(i%79),A));
  if(!diagonal&&i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.8,A));
  if(!diagonal&&i+1<n)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-.8,A));
 }
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));return A;
}
static void linear(const char *method,LIS_INT n,int repeat)
{
 LIS_MATRIX A=matrix(n,0);LIS_VECTOR b=NULL,x=NULL;LIS_SOLVER s=NULL;
 CHECK(!lis_vector_duplicate(A,&b));CHECK(!lis_vector_duplicate(A,&x));CHECK(!lis_solver_create(&s));
 for(LIS_INT i=0;i<n;i++)b->value[i]=sin(.09*i)+.1;
 char opts[256];const char *extra=!strcmp(method,"gcrodr")?" -gcrodr_keep 3":!strcmp(method,"fgcrodr")?" -fgcrodr_k 3":"";
 snprintf(opts,sizeof opts,"-i %s -p jacobi -restart 12 -maxiter 200 -tol 1e-10 -initx_zeros true%s",method,extra);
 CHECK(!lis_solver_set_option(opts,s));init_hash();
 double elapsed=0;
 for(int k=-1;k<repeat;k++){
  CHECK(!lis_vector_set_all(0,x));double t=now();(void)lis_solve(A,b,x,s);double dt=now()-t;
  LIS_INT status,iter;LIS_REAL resid;
  CHECK(!lis_solver_get_status(s,&status));CHECK(!lis_solver_get_iter(s,&iter));CHECK(!lis_solver_get_residualnorm(s,&resid));
  CHECK(status==LIS_SUCCESS&&isfinite(resid)&&resid<=1e-10);
  if(k>=0){elapsed+=dt;MIX(status);MIX(iter);MIX(resid);mix(x->value,(size_t)n*sizeof(LIS_SCALAR));}
 }
 line(method,n,repeat,elapsed);lis_solver_destroy(s);lis_vector_destroy(x);lis_vector_destroy(b);lis_matrix_destroy(A);
}
static void feast(int repeat)
{
 LIS_MATRIX A=matrix(4,1);LIS_VECTOR x=NULL;LIS_ESOLVER s=NULL;
 CHECK(!lis_vector_duplicate(A,&x));CHECK(!lis_esolver_create(&s));
 CHECK(!lis_esolver_set_option("-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16",s));
 init_hash();double elapsed=0;
 for(int k=-1;k<repeat;k++){
  LIS_SCALAR ev;CHECK(!lis_vector_set_all(1,x));double t=now();CHECK(!lis_esolve(A,x,&ev,s));double dt=now()-t;
  CHECK(s->retcode==LIS_SUCCESS&&s->nev==2);
  if(k>=0){elapsed+=dt;MIX(s->nev);for(LIS_INT j=0;j<s->nev;j++){MIX(s->evalue[j]);MIX(s->resid[j]);MIX(s->iter[j]);mix(s->evector[j]->value,4*sizeof(LIS_SCALAR));}}
 }
 line("feast",4,repeat,elapsed);lis_esolver_destroy(s);lis_vector_destroy(x);lis_matrix_destroy(A);
}
static int op(void *ctx,int64_t n,const double *x,double *y)
{(void)ctx;for(int64_t i=0;i<n;i++)y[i]=(1.+(double)i/n)*x[i];return 0;}
static void krylov(int repeat)
{
 const int64_t n=32768;lis_krylov_problem64 p={n,LIS_KRYLOV_STANDARD,0,NULL,op,NULL,NULL};
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=2;o.ncv=8;o.maxiter=8;o.tol=1e-9;
 double values[2],res[2],*vec=malloc(2*n*sizeof(double));CHECK(vec);init_hash();double elapsed=0;
 for(int k=-1;k<repeat;k++){
  memset(values,0,sizeof values);memset(res,0,sizeof res);memset(vec,0,2*n*sizeof(double));
  lis_krylov_stats st;double t=now();int rc=lis_krylov_symmetric_ex64(&p,&o,NULL,NULL,values,vec,res,&st,NULL);double dt=now()-t;
  CHECK(rc==LIS_KRYLOV_SUCCESS||rc==LIS_KRYLOV_NO_CONVERGENCE);
  if(k>=0){elapsed+=dt;MIX(rc);MIX(st.nreturned);MIX(st.iterations);MIX(st.op_count);mix(values,sizeof values);mix(res,sizeof res);mix(vec,2*n*sizeof(double));}
 }
 line("krylov",(LIS_INT)n,repeat,elapsed);free(vec);
}
int main(int argc,char **argv)
{
 CHECK(!lis_initialize(&argc,&argv));
 LIS_ABI_INFO client,lib;memset(&client,0,sizeof client);memset(&lib,0,sizeof lib);
 lis_abi_get_client_info(&client);CHECK(!lis_abi_get_info(&lib,sizeof lib));
 CHECK(!memcmp(&client,&lib,sizeof client));puts("BASELINE CLIENT ABI INFO MATCH");
 const char *methods[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int m=0;m<6;m++){linear(methods[m],37,100);linear(methods[m],32768,6);}
 feast(12);krylov(3);lis_finalize();return 0;
}
