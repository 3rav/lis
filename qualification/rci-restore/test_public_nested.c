/* Installed public API smoke test and genuinely nested FEAST -> GMRES RCI.
 * AI-assisted contribution, proposed under the LIS terms in COPYING. */
#include "lis_config.h"
#include "lis_rci_native_full.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define C(x) do { if(!(x)){fprintf(stderr,"PUBLIC NESTED FAIL %d: %s\n",__LINE__,#x);return 1;} } while(0)
static const double eigenvalues[4]={1,2,4,7};
static int block_calls,shift_calls,child_calls;
static uint64_t child_last;
static int solve_shift(lis_rci_native_full *child,const lis_rci_native_session_request *q,int fail)
{
 LIS_MATRIX M=NULL; LIS_VECTOR b=NULL,x=NULL; LIS_SOLVER s=NULL;
 LIS_INT err=LIS_ERR_ILL_ARG,api=0,code=0;
 lis_rci_native_session_request r={0};
 double ar=(double)q->shift_real,ai=(double)q->shift_imag;
 int generalized=q->B!=NULL;
 if(q->rows!=4||q->columns!=1||!(q->flags&LIS_RCI_NATIVE_SPLIT))return LIS_ERR_ILL_ARG;
 if(lis_matrix_create(LIS_COMM_WORLD,&M)||lis_matrix_set_size(M,0,8))goto done;
 for(LIS_INT i=0;i<4;i++){
  double mass=generalized?1+.2*i:1, d=mass*(ar-eigenvalues[i]), t=mass*ai;
  if(lis_matrix_set_value(LIS_INS_VALUE,i,i,d,M)||lis_matrix_set_value(LIS_INS_VALUE,i,i+4,-t,M)||
     lis_matrix_set_value(LIS_INS_VALUE,i+4,i,t,M)||lis_matrix_set_value(LIS_INS_VALUE,i+4,i+4,d,M))goto done;
 }
 if(lis_matrix_set_type(M,LIS_MATRIX_CSR)||lis_matrix_assemble(M)||lis_vector_duplicate(M,&b)||lis_vector_duplicate(M,&x)||lis_solver_create(&s))goto done;
 for(LIS_INT i=0;i<4;i++){
  LIS_SCALAR a,z;
  if(lis_vector_get_value((LIS_VECTOR)q->in[0],i,&a)||lis_vector_get_value((LIS_VECTOR)q->in[1],i,&z)||
     lis_vector_set_value(LIS_INS_VALUE,i,a,b)||lis_vector_set_value(LIS_INS_VALUE,i+4,z,b))goto done;
 }
 if(lis_vector_set_all(0,x)||lis_solver_set_option("-i gmres -p none -restart 8 -maxiter 100 -tol 1e-13 -initx_zeros true",s))goto done;
 int st=lis_rci_native_full_start_linear(child,M,b,x,s,&r,sizeof r);
 while(st>0){
  if(r.token<=child_last){lis_rci_native_full_cancel(child,r.token);goto done;}
  child_last=r.token;child_calls++;
  if(fail){lis_rci_native_full_cancel(child,r.token);err=LIS_BREAKDOWN;goto done;}
  LIS_INT cb=lis_rci_native_execute(&r,sizeof r);
  st=lis_rci_native_full_advance(child,r.token,cb,&r,sizeof r);
 }
 if(st!=LIS_RCI_DONE||lis_rci_native_full_status(child,&api,&code)!=LIS_RCI_DONE||api||code)goto done;
 for(LIS_INT i=0;i<4;i++){
  LIS_SCALAR a,z;
  if(lis_vector_get_value(x,i,&a)||lis_vector_get_value(x,i+4,&z)||
     lis_vector_set_value(LIS_INS_VALUE,i,a,(LIS_VECTOR)q->out[0])||lis_vector_set_value(LIS_INS_VALUE,i,z,(LIS_VECTOR)q->out[1]))goto done;
 }
 err=LIS_SUCCESS;
done:
 if(s)lis_solver_destroy(s);if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(M)lis_matrix_destroy(M);
 return err;
}
static int feast(lis_rci_native_full *parent,lis_rci_native_full *child,int generalized,int fail)
{
 LIS_MATRIX A=NULL,B=NULL;LIS_VECTOR x=NULL;LIS_ESOLVER e=NULL;LIS_SCALAR val;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,4));
 if(generalized){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,4));}
 for(LIS_INT i=0;i<4;i++){
  double mass=generalized?1+.2*i:1;
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,mass*eigenvalues[i],A));
  if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,mass,B));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 C(!lis_vector_duplicate(A,&x));C(!lis_vector_set_all(1,x));C(!lis_esolver_create(&e));
 C(!lis_esolver_set_option("-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16",e));
 lis_rci_native_session_request q={0};
 int st=lis_rci_native_full_start_feast(parent,A,B,x,&val,e,&q,sizeof q);
 uint64_t last=0;
 while(st>0){
  C(q.token>last);last=q.token;
  C(lis_rci_native_execute(&q,sizeof q-1)==LIS_ERR_ILL_ARG);
  lis_rci_native_session_request bad=q;bad.version++;
  C(lis_rci_native_execute(&bad,sizeof bad)==LIS_ERR_ILL_ARG);
  LIS_INT cb;
  if(st==LIS_RCI_OP_SHIFTED_SOLVE){shift_calls++;cb=solve_shift(child,&q,fail);}
  else{if(q.flags&LIS_RCI_NATIVE_BLOCK)block_calls++;cb=lis_rci_native_execute(&q,sizeof q);}
  if(!fail)C(!cb);
  st=lis_rci_native_full_advance(parent,q.token,cb,&q,sizeof q);
 }
 LIS_INT api=0,native=0;
 C(lis_rci_native_full_status(parent,&api,&native)==st);
 if(fail){C(st==LIS_RCI_CALLBACK&&native!=LIS_SUCCESS&&e->nev==0);}
 else{
  C(st==LIS_RCI_DONE&&api==LIS_SUCCESS&&native==LIS_SUCCESS&&e->nev==2);
  double a=(double)e->evalue[0],b=(double)e->evalue[1];
  C((fabs(a-2)<1e-6&&fabs(b-4)<1e-6)||(fabs(a-4)<1e-6&&fabs(b-2)<1e-6));
  C(e->resid[0]<1e-7&&e->resid[1]<1e-7);
 }
 lis_esolver_destroy(e);lis_vector_destroy(x);if(B)lis_matrix_destroy(B);lis_matrix_destroy(A);
 return 0;
}
int main(int argc,char **argv)
{
 C(!lis_initialize(&argc,&argv));
 lis_rci_native_full *p=lis_rci_native_full_create(),*c=lis_rci_native_full_create();
 C(p&&c);C(!feast(p,c,0,0));C(!feast(p,c,1,0));C(!feast(p,c,1,1));
 C(block_calls>0&&shift_calls>0&&child_calls>shift_calls);
 lis_rci_native_full_destroy(c);lis_rci_native_full_destroy(p);lis_finalize();
 printf("PUBLIC NESTED RCI PASS blocks=%d shifts=%d child_requests=%d\n",block_calls,shift_calls,child_calls);
 return 0;
}
