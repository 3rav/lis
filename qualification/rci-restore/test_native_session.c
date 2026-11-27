/* Exercise unified request ABI against all eight real-double native cores. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_esolver.h"
#include "lis_rci_native_session.h"
#include "lis_rci_native_linear_restore.h"
#include "lis_rci_native_gcro_restore.h"
#include "lis_rci_native_fgcro_restore.h"
#include "lis_rci_native_feast.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define C(x) do{if(!(x)){fprintf(stderr,"SESSION FAIL %d: %s\n",__LINE__,#x);exit(1);}}while(0)
static lis_rci_native_session *S;
static uint64_t last;
static int fault,terminal,requests;
static LIS_INT native_code;
static int consume(lis_rci_native_session_request *q,int st)
{
 int count=0;
 while(st>0){
  C(q->version==LIS_RCI_NATIVE_SESSION_VERSION&&q->size==sizeof(*q)&&q->operation==(uint32_t)st);
  C(q->token>last&&q->rows&&q->columns&&q->in&&q->out);last=q->token;count++;requests++;
  C(lis_rci_native_session_advance(S,q->token-1,0,NULL,q,sizeof(*q))==LIS_RCI_PROTOCOL);
  C(lis_rci_native_session_advance(S,q->token,0,NULL,q,sizeof(*q)-1)==LIS_RCI_BADARG);
  lis_rci_native_session_request copy={0};C(lis_rci_native_session_pending(S,&copy,sizeof copy)==st&&copy.token==q->token);
  if(fault&&count==2){
   if(fault==1)st=lis_rci_native_session_cancel(S,q->token);
   else if(fault==2)st=lis_rci_native_session_advance(S,q->token,0,NULL,q,sizeof(*q));
   else {q->A->value_generation++;st=lis_rci_native_session_advance(S,q->token,0,NULL,q,sizeof(*q));}
   C(st==(fault==1?LIS_RCI_CANCELLED:fault==2?LIS_RCI_NONFINITE:LIS_RCI_NATIVE_STALE));break;
  }
  C(q->flags&LIS_RCI_NATIVE_VECTORS);C(q->index_bytes==sizeof(LIS_INT)&&q->scalar_bytes==sizeof(LIS_SCALAR));
  LIS_INT err=0;
  if(q->operation==LIS_RCI_OP_SHIFTED_SOLVE){
   C(q->flags&LIS_RCI_NATIVE_SPLIT);C(!(q->flags&LIS_RCI_NATIVE_SHIFT_A_MINUS_SIGMA_B));
   lis_vector_set_all(0,(LIS_VECTOR)q->out[0]);lis_vector_set_all(0,(LIS_VECTOR)q->out[1]);
   err=lis_efeast_shifted_solve_split(q->A,q->B,q->shift_real,q->shift_imag,
    (LIS_VECTOR)q->in[0],(LIS_VECTOR)q->in[1],(LIS_VECTOR)q->out[0],(LIS_VECTOR)q->out[1],q->solver);
  }else if(q->flags&LIS_RCI_NATIVE_BLOCK){
   LIS_VECTOR *in=malloc(q->columns*sizeof(*in)),*out=malloc(q->columns*sizeof(*out));C(in&&out);
   for(uint64_t i=0;i<q->columns;i++){in[i]=(LIS_VECTOR)q->in[i];out[i]=(LIS_VECTOR)q->out[i];}
   err=lis_matvec_block(q->A,(LIS_INT)q->columns,in,out);free(in);free(out);
  }else{
   C(q->columns==1);
   err=q->operation==LIS_RCI_OP_PRECON?lis_psolve(q->solver,(LIS_VECTOR)q->in[0],(LIS_VECTOR)q->out[0]):
    lis_matvec(q->A,(LIS_VECTOR)q->in[0],(LIS_VECTOR)q->out[0]);
  }
  C(!err);st=lis_rci_native_session_advance(S,q->token,err,NULL,q,sizeof(*q));
 }
 terminal=st;C(lis_rci_native_session_status(S,&native_code)==st);
 C(lis_rci_native_session_pending(S,q,sizeof(*q))==LIS_RCI_STATE&&!q->token);
 return st;
}
static LIS_INT linear_driver(LIS_SOLVER s,LIS_INT method)
{lis_rci_native_session_request q;int st=lis_rci_native_session_start_linear(S,s,method,&q,sizeof q);consume(&q,st);return native_code;}
static LIS_INT gcro_driver(LIS_SOLVER s){return linear_driver(s,LIS_SOLVER_GCRODR);}
static LIS_INT fgcro_driver(LIS_SOLVER s){return linear_driver(s,LIS_SOLVER_FGCRODR);}
static LIS_INT feast_driver(LIS_ESOLVER e)
{lis_rci_native_session_request q;int st=lis_rci_native_session_start_feast(S,e,&q,sizeof q);consume(&q,st);return native_code;}
static void linear(const char *method)
{
 for(fault=0;fault<=3;fault++){
  LIS_MATRIX A;LIS_VECTOR b,x;LIS_SOLVER s;
  C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,17));
  for(LIS_INT i=0;i<17;i++){
   C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,3+.02*i,A));
   if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.5,A));
   if(i<16)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-.5,A));
  }
  C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
  C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_solver_create(&s));
  C(!lis_vector_set_all(1,b));C(!lis_vector_set_all(0,x));
  char opts[256];snprintf(opts,sizeof opts,"-i %s -p jacobi -restart 5 -maxiter 200 -tol 1e-10 -initx_zeros true%s",
   method,!strcmp(method,"gcrodr")?" -gcrodr_keep 2":!strcmp(method,"fgcrodr")?" -fgcrodr_k 2":"");
  C(!lis_solver_set_option(opts,s));
  terminal=99;int before=requests;
  (void)lis_solve(A,b,x,s);
  C(requests>before);
  if(!fault){C(terminal==LIS_RCI_DONE&&s->retcode==LIS_SUCCESS&&s->resid<=1e-10);}
  else C(terminal<0&&s->retcode!=LIS_SUCCESS);
  lis_solver_destroy(s);lis_vector_destroy(x);lis_vector_destroy(b);lis_matrix_destroy(A);
 }
 printf("SESSION PASS %s\n",method);
}
static void feast(void)
{
 for(fault=0;fault<=3;fault++){
  LIS_MATRIX A;LIS_VECTOR x;LIS_ESOLVER e;LIS_SCALAR value;
  C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,4));
  double values[4]={1,2,4,7};for(int i=0;i<4;i++)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,values[i],A));
  C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
  C(!lis_vector_duplicate(A,&x));C(!lis_vector_set_all(1,x));C(!lis_esolver_create(&e));
  C(!lis_esolver_set_option("-e feast -ss 4 -emaxiter 12 -etol 1e-7 -feast_emin 1.5 -feast_emax 4.5 -feast_ratio 0.30 -feast_points 16",e));
  /* Inner blocking GMRES uses its normal native path, not the same facade. */
  lis_rci_native_linear_test_select(0,NULL);
  terminal=99;int before=requests;(void)lis_esolve(A,x,&value,e);C(requests>before);
  if(!fault)C(terminal==LIS_RCI_DONE&&e->retcode==LIS_SUCCESS&&e->nev==2);
  else C(terminal<0&&e->retcode!=LIS_SUCCESS&&e->nev==0);
  lis_esolver_destroy(e);lis_vector_destroy(x);lis_matrix_destroy(A);
 }
 puts("SESSION PASS feast");
}
static int op(void *ctx,int64_t n,const double *x,double *y){(void)ctx;for(int64_t i=0;i<n;i++)y[i]=(1+.11*i)*x[i];return 0;}
static void krylov(void)
{
 lis_krylov_problem64 p={24,LIS_KRYLOV_STANDARD,0,NULL,op,NULL,NULL};
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=3;o.ncv=9;o.maxiter=120;o.tol=1e-9;
 double val[3],x[72],res[3],vref[3],xref[72],rref[3];lis_krylov_stats stats,ref;
 int expected=lis_krylov_symmetric_ex64(&p,&o,NULL,NULL,vref,xref,rref,&ref,NULL);
 lis_rci_native_generations g={1,2,3,4};
 for(int f=0;f<3;f++){
  lis_rci_native_session_request q={0};int st=lis_rci_native_session_start_krylov(S,&p,&o,NULL,NULL,NULL,val,x,res,&stats,NULL,&g,&q,sizeof q);
  while(st>0){
   C(q.token>last);last=q.token;requests++;
   C(q.flags&LIS_RCI_NATIVE_RAW_DOUBLE);C(q.operation==LIS_RCI_OP_A&&q.rows==24&&q.columns==1);
   C(lis_rci_native_session_advance(S,q.token,0,NULL,&q,sizeof q)==LIS_RCI_BADARG);
   if(f==1){C(lis_rci_native_session_cancel(S,q.token)==LIS_RCI_CANCELLED);break;}
   if(f==2){lis_rci_native_generations bad=g;bad.a_values++;C(lis_rci_native_session_advance(S,q.token,0,&bad,&q,sizeof q)==LIS_RCI_NATIVE_STALE);break;}
   C(!op(NULL,q.rows,q.in[0],q.out[0]));st=lis_rci_native_session_advance(S,q.token,0,&g,&q,sizeof q);
  }
  if(!f){C(lis_rci_native_session_status(S,&native_code)==st&&native_code==expected);
   C(stats.nreturned==ref.nreturned&&stats.iterations==ref.iterations);
   C(!memcmp(val,vref,sizeof val)&&!memcmp(x,xref,sizeof x)&&!memcmp(res,rref,sizeof res));}
 }
 puts("SESSION PASS krylov");
}
int main(int argc,char **argv)
{
 C(!lis_initialize(&argc,&argv));S=lis_rci_native_session_create();C(S);
 lis_rci_native_linear_test_select(0,linear_driver);lis_rci_native_gcro_test_select(0,gcro_driver);lis_rci_native_fgcro_test_select(0,fgcro_driver);
 const char *methods[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int i=0;i<6;i++)linear(methods[i]);
 lis_rci_native_gcro_test_select(0,NULL);lis_rci_native_fgcro_test_select(0,NULL);
 lis_rci_native_feast_test_select(0,feast_driver);feast();lis_rci_native_feast_test_select(0,NULL);krylov();
 lis_rci_native_session_destroy(S);lis_finalize();printf("UNIFIED NATIVE SESSION PASS requests=%d last_token=%llu\n",requests,(unsigned long long)last);return 0;
}
