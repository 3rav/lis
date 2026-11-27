/* Differential native full FEAST qualification. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_esolver.h"
#include "lis_rci_native_feast.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 4
#define C(x) do{if(!(x)){fprintf(stderr,"FEAST FULL FAIL %d: %s\n",__LINE__,#x);rc=1;goto cleanup;}}while(0)
typedef struct {LIS_INT nev,iter[4],iter2[4];LIS_REAL value[4],resid[4],vectors[16];} answer;
static int stop,kind,observed,terminal;
static LIS_INT external(LIS_ESOLVER e)
{
 lis_rci_native_feast *K=lis_rci_native_feast_create();if(!K)return LIS_OUT_OF_MEMORY;
 lis_rci_native_feast_request q;
 int st=lis_rci_native_feast_start(K,e,1,&q);uint64_t last=0;observed=0;
 while(st>0){
  if(q.token<=last||lis_rci_native_feast_advance(K,q.token-1,0,&q)!=NF_PROTOCOL){st=NF_ERROR;break;}
  last=q.token;
  if(++observed==stop){
   if(kind==1)st=lis_rci_native_feast_cancel(K,q.token);
   else if(kind==2)st=lis_rci_native_feast_advance(K,q.token,0,&q);
   else st=lis_rci_native_feast_advance(K,q.token,LIS_ERR_ILL_ARG,&q);
   break;
  }
  LIS_INT err;
  if(q.operation==NF_SHIFT){
   lis_vector_set_all(0,q.out_real);lis_vector_set_all(0,q.out_imag);
   err=lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver);
  }else err=!q.block?lis_matvec(q.A,q.in[0],q.out[0]):lis_matvec_block(q.A,q.count,q.in,q.out);
  st=lis_rci_native_feast_advance(K,q.token,err,&q);
 }
 terminal=st;
 LIS_INT rc=lis_rci_native_feast_native(K);lis_rci_native_feast_destroy(K);return rc;
}
static int solve(int mode,int generalized,int fail,int fault,answer *a)
{
 LIS_MATRIX A=NULL,B=NULL;LIS_VECTOR x=NULL;LIS_ESOLVER e=NULL;
 int rc=0;LIS_SCALAR value=0;memset(a,0,sizeof(*a));
 const double lambda[N]={1,2,4,7};
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 if(generalized){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,N));}
 for(LIS_INT i=0;i<N;i++){
  double b=generalized?1+.2*i:1;
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b*lambda[i],A));
  if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b,B));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 C(!lis_vector_duplicate(A,&x));C(!lis_vector_set_all(1,x));
 C(!lis_esolver_create(&e));
 C(!lis_esolver_set_option("-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16",e));
 stop=fail;kind=fault;terminal=0;observed=0;
 lis_rci_native_feast_test_select(mode==0,mode==2?external:NULL);
 LIS_INT err=B?lis_gesolve(A,B,x,&value,e):lis_esolve(A,x,&value,e);
 if(fail){C(err!=LIS_SUCCESS&&e->nev==0);C(terminal==(fault==1?NF_CANCELLED:fault==2?NF_NONFINITE:NF_ERROR));goto cleanup;}
 C(!err);C(e->nev==2);
 a->nev=e->nev;
 for(LIS_INT i=0;i<e->nev;i++){
  a->value[i]=e->evalue[i];a->resid[i]=e->resid[i];a->iter[i]=e->iter[i];a->iter2[i]=e->iter2[i];
  C(e->evector[i]);memcpy(a->vectors+N*i,e->evector[i]->value,N*sizeof(LIS_SCALAR));
  C(fabs(a->value[i]-(i?4:2))<1e-6&&a->resid[i]<1e-7);
 }
cleanup:
 lis_rci_native_feast_test_select(0,NULL);
 if(e)lis_esolver_destroy(e);if(x)lis_vector_destroy(x);if(B)lis_matrix_destroy(B);if(A)lis_matrix_destroy(A);return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 int rc=0;answer a,b,c,scratch;
 for(int gen=0;gen<=1;gen++){
  C(!solve(0,gen,0,0,&a));C(!solve(1,gen,0,0,&b));C(!solve(2,gen,0,0,&c));
  C(!memcmp(&a,&b,sizeof(a))&&!memcmp(&a,&c,sizeof(a)));
  int n=observed;C(n>0);
  /* First/last plus every boundary within the first 32 requests. */
  for(int f=1;f<=n;f++)if(f<=32||f==n)for(int k=1;k<=3;k++)C(!solve(2,gen,f,k,&scratch));
  printf("NATIVE FULL FEAST PASS generalized=%d requests=%d faults=%d\n",gen,n,3*(n<=32?n:33));
 }
cleanup:
 lis_finalize();return rc;
}
