/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
/* Native LIS bridge exercise beyond the 4096-element RCI OpenMP threshold. */
#include "lis_rci_bridge.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define N 8192
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c);rc=1;goto cleanup;}}while(0)
static void product(const double*x,double*y){for(size_t i=0;i<N;i++)y[i]=(4.0+(i%17)/17.0)*x[i]-(i?x[i-1]:0)-(i+1<N?x[i+1]:0);}
int main(int argc,char**argv){
 LIS_MATRIX A=NULL;lis_rci_bridge*p=NULL;int rc=0;double*mem=NULL;
 if(lis_initialize(&argc,&argv))return 1;
#ifdef _OPENMP
 int threads=0;
#pragma omp parallel reduction(+:threads)
 {threads++;}
 printf("OpenMP actual team=%d\n",threads);
 const char*expect=getenv("OMP_NUM_THREADS");if(expect)CHECK(threads==atoi(expect));
#endif
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,N));
 for(LIS_INT i=0;i<N;i++){CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0+(i%17)/17.0,A));if(i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1,A));if(i+1<N)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1,A));}
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 mem=calloc(4*N,sizeof(double));CHECK(mem!=NULL);double*x=mem,*b=x+N,*exact=b+N,*ax=exact+N;
 for(int m=1;m<=4;m++){
  lis_rci_options o={N,200,30,1e-10,0,(lis_rci_method)m};CHECK(!lis_rci_bridge_create(A,&o,1,&p));
  for(int rep=0;rep<2;rep++){
   for(size_t i=0;i<N;i++){exact[i]=sin(.07*i+rep)+.3*cos(.19*(i+3));x[i]=rep?.125:0;}
   product(exact,b);lis_rci_info info;int s=lis_rci_bridge_solve(p,b,x,&info);product(x,ax);
   long double r=0,bn=0,e=0,en=0;for(size_t i=0;i<N;i++){long double d=(long double)b[i]-ax[i],t=(long double)x[i]-exact[i];r+=d*d;bn+=(long double)b[i]*b[i];e+=t*t;en+=(long double)exact[i]*exact[i];}
   double rr=(double)sqrtl(r/bn),er=(double)sqrtl(e/en);int ok=s==LIS_RCI_DONE&&rr<=1.05e-10&&er<=1e-8&&isfinite(rr)&&isfinite(er);
   printf("LARGE method=%d rhs=%d iter=%zu residual=%.9e error=%.9e %s\n",m,rep,info.iterations,rr,er,ok?"PASS":"FAIL");CHECK(ok);
  }
  lis_rci_bridge_destroy(p);p=NULL;
 }
cleanup:
 free(mem);lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);lis_finalize();return rc;
}
