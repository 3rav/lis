/* Public setup/cleanup differential qualification, real double. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_esolver.h"
#include "lis_rci_native_full.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define C(x) do{if(!(x)){fprintf(stderr,"FULL LIFECYCLE FAIL %d: %s\n",__LINE__,#x);exit(1);}}while(0)
static lis_rci_native_full *S;static uint64_t last;static int request_count,inside_request,external_mode,hidden_calls;
static LIS_INT execute(lis_rci_native_session_request *q)
{
 LIS_INT rc=0;inside_request=1;
 if(q->operation==LIS_RCI_OP_SHIFTED_SOLVE){
  C(q->flags&LIS_RCI_NATIVE_SPLIT);
  lis_vector_set_all(0,(LIS_VECTOR)q->out[0]);lis_vector_set_all(0,(LIS_VECTOR)q->out[1]);
  rc=lis_efeast_shifted_solve_split(q->A,q->B,q->shift_real,q->shift_imag,(LIS_VECTOR)q->in[0],(LIS_VECTOR)q->in[1],(LIS_VECTOR)q->out[0],(LIS_VECTOR)q->out[1],q->solver);
 }else if(q->flags&LIS_RCI_NATIVE_BLOCK){
  LIS_VECTOR *in=malloc(q->columns*sizeof(*in)),*out=malloc(q->columns*sizeof(*out));C(in&&out);
  for(uint64_t i=0;i<q->columns;i++){in[i]=(LIS_VECTOR)q->in[i];out[i]=(LIS_VECTOR)q->out[i];}
  rc=lis_matvec_block(q->A,(LIS_INT)q->columns,in,out);free(in);free(out);
 }else rc=q->operation==LIS_RCI_OP_PRECON?lis_psolve(q->solver,(LIS_VECTOR)q->in[0],(LIS_VECTOR)q->out[0]):lis_matvec(q->A,(LIS_VECTOR)q->in[0],(LIS_VECTOR)q->out[0]);
 inside_request=0;return rc;
}
static int drive(int st,lis_rci_native_session_request *q,int stop,int kind,LIS_INT *native)
{
 request_count=0;
 while(st>0){
  C(q->token>last);last=q->token;request_count++;
  C(lis_rci_native_full_advance(S,q->token-1,0,q,sizeof *q)==LIS_RCI_PROTOCOL);
  C(lis_rci_native_full_advance(S,q->token,0,q,sizeof *q-1)==LIS_RCI_BADARG);
  if(request_count==stop){
   if(kind==1)st=lis_rci_native_full_cancel(S,q->token);
   else if(kind==2)st=lis_rci_native_full_advance(S,q->token,0,q,sizeof *q);
   else st=lis_rci_native_full_advance(S,q->token,LIS_ERR_ILL_ARG,q,sizeof *q);
   C(st==(kind==1?LIS_RCI_CANCELLED:kind==2?LIS_RCI_NONFINITE:LIS_RCI_CALLBACK));break;
  }
  LIS_INT code=execute(q);C(!code);
  st=lis_rci_native_full_advance(S,q->token,code,q,sizeof *q);
 }
 LIS_INT api=0;C(lis_rci_native_full_status(S,&api,native)==st);if(!stop)C(api==LIS_SUCCESS);
 return st;
}
typedef struct{LIS_INT code,iter;LIS_REAL resid;LIS_SCALAR x[17];} result;
static LIS_INT user_op(void *ctx,const LIS_SCALAR *x,LIS_SCALAR *y)
{
 (void)ctx;if(external_mode&&!inside_request){hidden_calls++;return 91;}
 for(int i=0;i<17;i++)y[i]=(3+.02*i)*x[i]-(i?.5*x[i-1]:0)-(i<16?.5*x[i+1]:0);
 return LIS_SUCCESS;
}
static void linear(const char *method,int mode,int user,int stop,int fault,result *r)
{
 LIS_MATRIX A;LIS_VECTOR b,x;LIS_SOLVER s;memset(r,0,sizeof *r);
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,17));
 if(user){C(!lis_matrix_set_user(A,NULL,user_op,user_op));}
 else{
  for(int i=0;i<17;i++){C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,3+.02*i,A));if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.5,A));if(i<16)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-.5,A));}
  C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 }
 C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_solver_create(&s));
 C(!lis_vector_set_all(1,b));C(!lis_vector_set_all(.125,x));
 char opts[256];snprintf(opts,sizeof opts,"-i %s -p %s -restart 5 -maxiter 200 -tol 1e-10 -initx_zeros false%s",method,user?"none":"jacobi",
 !strcmp(method,"gcrodr")?" -gcrodr_keep 2":!strcmp(method,"fgcrodr")?" -fgcrodr_k 2":"");
 C(!lis_solver_set_option(opts,s));external_mode=mode;hidden_calls=0;
 if(!mode){(void)lis_solve(A,b,x,s);r->code=s->retcode;}
 else{lis_rci_native_session_request q={0};int st=lis_rci_native_full_start_linear(S,A,b,x,s,&q,sizeof q);st=drive(st,&q,stop,fault,&r->code);
  if(!stop)C(st==LIS_RCI_DONE);else C(st<0&&s->retcode!=LIS_SUCCESS);}
 C(!hidden_calls);
 r->iter=s->iter;r->resid=s->resid;memcpy(r->x,x->value,sizeof r->x);
 if(!stop)C(r->code==LIS_SUCCESS&&r->resid<=1e-10);
 external_mode=0;lis_solver_destroy(s);lis_vector_destroy(x);lis_vector_destroy(b);lis_matrix_destroy(A);
}
typedef struct{LIS_INT code,nev,iter[4];LIS_SCALAR value[4],vectors[16];LIS_REAL resid[4];} eigresult;
static void feast(int gen,int mode,int stop,int fault,eigresult *r)
{
 LIS_MATRIX A,B=NULL;LIS_VECTOR x;LIS_ESOLVER e;LIS_SCALAR value;memset(r,0,sizeof *r);
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,4));
 if(gen){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,4));}
 double eig[4]={1,2,4,7};for(int i=0;i<4;i++){double b=gen?1+.2*i:1;C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,eig[i]*b,A));if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b,B));}
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 C(!lis_vector_duplicate(A,&x));C(!lis_vector_set_all(1,x));C(!lis_esolver_create(&e));
 C(!lis_esolver_set_option("-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16",e));
 if(!mode){C(!lis_gesolve(A,B,x,&value,e));r->code=e->retcode;}
 else{lis_rci_native_session_request q;int st=lis_rci_native_full_start_feast(S,A,B,x,&value,e,&q,sizeof q);st=drive(st,&q,stop,fault,&r->code);
  if(!stop)C(st==LIS_RCI_DONE);else C(st<0&&e->nev==0&&e->retcode!=LIS_SUCCESS);}
 r->nev=e->nev;
 for(int i=0;i<e->nev;i++){r->value[i]=e->evalue[i];r->resid[i]=e->resid[i];r->iter[i]=e->iter[i];memcpy(r->vectors+4*i,e->evector[i]->value,4*sizeof(LIS_SCALAR));}
 if(!stop)C(r->code==LIS_SUCCESS&&r->nev==2);
 lis_esolver_destroy(e);lis_vector_destroy(x);if(B)lis_matrix_destroy(B);lis_matrix_destroy(A);
}
int main(int argc,char **argv)
{
 C(!lis_initialize(&argc,&argv));S=lis_rci_native_full_create();C(S);
 const char *methods[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int i=0;i<6;i++)for(int user=0;user<=1;user++){
  /* Native recycling may impose explicit-matrix restrictions. */
  if(user&&i>=4)continue;
  result a,b,scratch;linear(methods[i],0,user,0,0,&a);linear(methods[i],1,user,0,0,&b);
  C(!memcmp(&a,&b,sizeof a));int n=request_count;C(n>1);
  for(int k=1;k<=3;k++){linear(methods[i],1,user,1,k,&scratch);linear(methods[i],1,user,n,k,&scratch);}
  printf("FULL LIFECYCLE PASS %s user=%d requests=%d\n",methods[i],user,n);
 }
 for(int gen=0;gen<=1;gen++){eigresult a,b,scratch;feast(gen,0,0,0,&a);feast(gen,1,0,0,&b);C(!memcmp(&a,&b,sizeof a));int n=request_count;
  for(int k=1;k<=3;k++){feast(gen,1,1,k,&scratch);feast(gen,1,n,k,&scratch);}
  printf("FULL LIFECYCLE PASS feast gen=%d requests=%d\n",gen,n);}
 lis_rci_native_full_destroy(S);lis_finalize();puts("FULL NATIVE LIFECYCLE PASS");return 0;
}
