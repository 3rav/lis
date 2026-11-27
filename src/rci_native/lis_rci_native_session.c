/* Unified qualification-only facade over all eight native continuations. */
#include "lis_rci_native_session.h"
#include "lis_rci_native_linear_restore.h"
#include "lis_rci_native_gcro_restore.h"
#include "lis_rci_native_fgcro_restore.h"
#include "lis_rci_native_feast.h"
#include "lis_rci_krylov_restore.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>
enum { RS_LINEAR=1,RS_GCRO,RS_FGCRO,RS_FEAST,RS_KRYLOV };
struct lis_rci_native_session {
 int kind,status,pending;
 uint64_t token,child_token;
 LIS_INT native;
 LIS_MATRIX A,B;
 LIS_SOLVER solver;
 LIS_ESOLVER esolver;
 lis_rci_native_generations generations;
 size_t capacity;
 const void **in;
 void **out;
 lis_rci_native_session_request request;
 union {
  lis_rci_native_linear_restore *l;
  lis_rci_native_gcro_restore *g;
  lis_rci_native_fgcro_restore *f;
  lis_rci_native_feast *e;
  lis_rci_krylov_restore *k;
 } child;
};
static int rs_buffer(lis_rci_native_session_request *q,size_t n){return q&&n>=sizeof(*q);}
static lis_rci_native_generations rs_generations(LIS_MATRIX A,LIS_MATRIX B)
{
 lis_rci_native_generations g={0,0,0,0};
 if(A){g.a_structure=A->structure_generation;g.a_values=A->value_generation;}
 if(B){g.b_structure=B->structure_generation;g.b_values=B->value_generation;}
 return g;
}
static int rs_equal(lis_rci_native_generations a,lis_rci_native_generations b)
{return a.a_structure==b.a_structure&&a.a_values==b.a_values&&a.b_structure==b.b_structure&&a.b_values==b.b_values;}
static int rs_current(const lis_rci_native_session *S,const lis_rci_native_generations *g)
{
 if(S->kind==RS_KRYLOV)return g&&rs_equal(S->generations,*g);
 if(S->solver&&S->solver->A!=S->A)return 0;
 if(S->esolver&&(S->esolver->A!=S->A||S->esolver->B!=S->B))return 0;
 return rs_equal(S->generations,rs_generations(S->A,S->B));
}
static void rs_release(lis_rci_native_session *S)
{
 switch(S->kind){
 case RS_LINEAR:lis_rci_native_linear_restore_destroy(S->child.l);break;
 case RS_GCRO:lis_rci_native_gcro_restore_destroy(S->child.g);break;
 case RS_FGCRO:lis_rci_native_fgcro_restore_destroy(S->child.f);break;
 case RS_FEAST:lis_rci_native_feast_destroy(S->child.e);break;
 case RS_KRYLOV:lis_rci_krylov_restore_destroy(S->child.k);break;
 }
 S->kind=0;S->pending=0;S->solver=NULL;S->esolver=NULL;S->A=S->B=NULL;
 free(S->in);free(S->out);S->in=NULL;S->out=NULL;S->capacity=0;
 memset(&S->child,0,sizeof S->child);memset(&S->request,0,sizeof S->request);
}
static int rs_init(lis_rci_native_session *S,int kind,size_t capacity)
{
 if(!S||!capacity||capacity>SIZE_MAX/sizeof(void *))return LIS_RCI_BADARG;
 if(S->pending)return LIS_RCI_STATE;
 rs_release(S);S->status=LIS_RCI_STATE;S->native=0;
 S->in=calloc(capacity,sizeof(*S->in));S->out=calloc(capacity,sizeof(*S->out));
 if(!S->in||!S->out){rs_release(S);return LIS_RCI_BADARG;}
 S->capacity=capacity;S->kind=kind;
 int ok=0;
 switch(kind){
 case RS_LINEAR:ok=(S->child.l=lis_rci_native_linear_restore_create())!=NULL;break;
 case RS_GCRO:ok=(S->child.g=lis_rci_native_gcro_restore_create())!=NULL;break;
 case RS_FGCRO:ok=(S->child.f=lis_rci_native_fgcro_restore_create())!=NULL;break;
 case RS_FEAST:ok=(S->child.e=lis_rci_native_feast_create())!=NULL;break;
 case RS_KRYLOV:ok=(S->child.k=lis_rci_krylov_restore_create())!=NULL;break;
 }
 if(!ok){rs_release(S);return LIS_RCI_BADARG;}
 return LIS_RCI_DONE;
}
static int rs_native(lis_rci_native_session *S,int st)
{
 uint64_t a=0,p=0;
 switch(S->kind){
 case RS_LINEAR:lis_rci_native_linear_restore_status(S->child.l,&S->native,&a,&p);break;
 case RS_GCRO:lis_rci_native_gcro_restore_status(S->child.g,&S->native,&a,&p);break;
 case RS_FGCRO:lis_rci_native_fgcro_restore_status(S->child.f,&S->native,&a,&p);break;
 case RS_FEAST:
  S->native=lis_rci_native_feast_native(S->child.e);
  switch(st){
   case NF_DONE:st=LIS_RCI_DONE;break;
   case NF_PROTOCOL:st=LIS_RCI_PROTOCOL;break;
   case NF_CANCELLED:st=LIS_RCI_CANCELLED;break;
   case NF_NONFINITE:st=LIS_RCI_NONFINITE;break;
   case NF_STATE:st=LIS_RCI_STATE;break;
   case NF_BADARG:st=LIS_RCI_BADARG;break;
   case NF_STALE:st=LIS_RCI_NATIVE_STALE;break;
   default:st=S->native==LIS_MAXITER?LIS_RCI_MAXITER:S->native==LIS_BREAKDOWN?LIS_RCI_BREAKDOWN:LIS_RCI_CALLBACK;
  }break;
 case RS_KRYLOV:{int code=0;lis_krylov_budget_stats counts;lis_rci_krylov_restore_status(S->child.k,&code,&counts);S->native=code;break;}
 }
 if(st==LIS_RCI_BADARG&&S->native==0)S->native=S->kind==RS_KRYLOV?LIS_KRYLOV_BAD_ARGUMENT:LIS_ERR_ILL_ARG;
 S->pending=0;memset(&S->request,0,sizeof S->request);
 return S->status=st;
}
static int rs_cancel_child(lis_rci_native_session *S)
{
 int st=LIS_RCI_STATE;
 switch(S->kind){
 case RS_LINEAR:st=lis_rci_native_linear_restore_cancel(S->child.l,S->child_token);break;
 case RS_GCRO:st=lis_rci_native_gcro_restore_cancel(S->child.g,S->child_token);break;
 case RS_FGCRO:st=lis_rci_native_fgcro_restore_cancel(S->child.f,S->child_token);break;
 case RS_FEAST:st=lis_rci_native_feast_cancel(S->child.e,S->child_token);break;
 case RS_KRYLOV:st=lis_rci_krylov_restore_cancel(S->child.k,S->child_token);break;
 }
 return rs_native(S,st);
}
static int rs_abort(lis_rci_native_session *S,int st)
{
 rs_cancel_child(S);
 S->native=S->kind==RS_KRYLOV?LIS_KRYLOV_CALLBACK_ERROR:LIS_BREAKDOWN;
 return S->status=st;
}
static int rs_publish(lis_rci_native_session *S,uint32_t op,uint32_t flags,uint64_t rows,uint64_t columns,
 LIS_REAL zr,LIS_REAL zi,LIS_MATRIX A,LIS_MATRIX B,LIS_SOLVER solver,lis_rci_native_session_request *q)
{
 if(S->token==UINT64_MAX){int st=rs_abort(S,LIS_RCI_PROTOCOL);*q=S->request;return st;}
 S->request=(lis_rci_native_session_request){
  .version=LIS_RCI_NATIVE_SESSION_VERSION,.size=sizeof(*q),.operation=op,.flags=flags,
  .index_bytes=(flags&LIS_RCI_NATIVE_RAW_DOUBLE)?sizeof(int64_t):sizeof(LIS_INT),
  .scalar_bytes=(flags&LIS_RCI_NATIVE_RAW_DOUBLE)?sizeof(double):sizeof(LIS_SCALAR),
  .real_bytes=sizeof(LIS_REAL),.token=++S->token,.rows=rows,.columns=columns,
  .generations=S->generations,.in=S->in,.out=S->out,.shift_real=zr,.shift_imag=zi,.A=A,.B=B,.solver=solver};
 S->pending=1;S->status=(int)op;*q=S->request;return S->status;
}
static int rs_linear_reply(lis_rci_native_session *S,int st,lis_rci_native_linear_request *r,lis_rci_native_session_request *q)
{
 if(st<=0){int code=rs_native(S,st);*q=S->request;return code;}
 S->child_token=r->token;S->in[0]=r->in;S->out[0]=r->out;
 if(!rs_current(S,NULL)){int code=rs_abort(S,LIS_RCI_NATIVE_STALE);*q=S->request;return code;}
 return rs_publish(S,r->operation,LIS_RCI_NATIVE_VECTORS,r->in->n,1,0,0,r->A,NULL,r->solver,q);
}
static int rs_feast_reply(lis_rci_native_session *S,int st,lis_rci_native_feast_request *r,lis_rci_native_session_request *q)
{
 if(st<=0){int code=rs_native(S,st);*q=S->request;return code;}
 S->child_token=r->token;
 if(!rs_current(S,NULL)){int code=rs_abort(S,LIS_RCI_NATIVE_STALE);*q=S->request;return code;}
 if(st==NF_SHIFT){
  S->in[0]=r->in_real;S->in[1]=r->in_imag;S->out[0]=r->out_real;S->out[1]=r->out_imag;
  return rs_publish(S,LIS_RCI_OP_SHIFTED_SOLVE,LIS_RCI_NATIVE_VECTORS|LIS_RCI_NATIVE_SPLIT,
   r->in_real->n,1,r->zr,r->zi,r->A,r->B,r->solver,q);
 }
 if(r->count<=0||(uint64_t)r->count>S->capacity){int code=rs_abort(S,LIS_RCI_BADARG);*q=S->request;return code;}
 for(LIS_INT i=0;i<r->count;i++){S->in[i]=r->in[i];S->out[i]=r->out[i];}
 return rs_publish(S,st==NF_A?LIS_RCI_OP_A:LIS_RCI_OP_B,
  LIS_RCI_NATIVE_VECTORS|(r->block?LIS_RCI_NATIVE_BLOCK:0),r->in[0]->n,r->count,0,0,r->A,NULL,NULL,q);
}
static int rs_krylov_reply(lis_rci_native_session *S,int st,lis_rci_krylov_restore_request *r,lis_rci_native_session_request *q)
{
 if(st<=0){int code=rs_native(S,st);*q=S->request;return code;}
 S->child_token=r->token;S->in[0]=r->in;S->out[0]=r->out;
 return rs_publish(S,r->operation,LIS_RCI_NATIVE_RAW_DOUBLE|
  (r->operation==LIS_RCI_OP_SHIFTED_SOLVE?LIS_RCI_NATIVE_SHIFT_A_MINUS_SIGMA_B:0),r->n,1,r->sigma,0,NULL,NULL,NULL,q);
}
lis_rci_native_session *lis_rci_native_session_create(void)
{lis_rci_native_session *S=calloc(1,sizeof(*S));if(S)S->status=LIS_RCI_STATE;return S;}
void lis_rci_native_session_destroy(lis_rci_native_session *S)
{if(S){if(S->pending)rs_cancel_child(S);rs_release(S);free(S);}}
int lis_rci_native_session_start_linear(lis_rci_native_session *S,LIS_SOLVER solver,LIS_INT method,lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!solver||!rs_buffer(q,bytes))return LIS_RCI_BADARG;
 int kind=RS_LINEAR;
 if(method==LIS_SOLVER_GCRODR)kind=RS_GCRO;
 else if(method==LIS_SOLVER_FGCRODR)kind=RS_FGCRO;
 else if(method!=LIS_SOLVER_CG&&method!=LIS_SOLVER_BICGSTAB&&method!=LIS_SOLVER_GMRES&&method!=LIS_SOLVER_FGMRES)return LIS_RCI_BADARG;
 int st=rs_init(S,kind,1);if(st)return st;
 S->solver=solver;S->A=solver->A;S->generations=rs_generations(S->A,NULL);
 lis_rci_native_linear_request r={0};
 switch(kind){
 case RS_LINEAR:st=lis_rci_native_linear_restore_start_method(S->child.l,solver,method,1,&r);break;
 case RS_GCRO:st=lis_rci_native_gcro_restore_start(S->child.g,solver,1,&r);break;
 default:st=lis_rci_native_fgcro_restore_start(S->child.f,solver,1,&r);break;
 }
 return rs_linear_reply(S,st,&r,q);
}
int lis_rci_native_session_start_feast(lis_rci_native_session *S,LIS_ESOLVER e,lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!e||!rs_buffer(q,bytes))return LIS_RCI_BADARG;
 LIS_INT n=e->options[LIS_EOPTIONS_SUBSPACE];if(n<=0)return LIS_RCI_BADARG;
 int st=rs_init(S,RS_FEAST,(size_t)(n>2?n:2));if(st)return st;
 S->esolver=e;S->A=e->A;S->B=e->B;S->generations=rs_generations(S->A,S->B);
 lis_rci_native_feast_request r={0};st=lis_rci_native_feast_start(S->child.e,e,1,&r);
 return rs_feast_reply(S,st,&r,q);
}
int lis_rci_native_session_start_krylov(lis_rci_native_session *S,const lis_krylov_problem64 *p,const lis_krylov_options *o,
 const lis_krylov_ext_options *ext,lis_krylov_workspace *w,const lis_krylov_budget *budget,double *v,double *x,double *r,
 lis_krylov_stats *stats,lis_krylov_ext_stats *extra,const lis_rci_native_generations *g,lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!g||!rs_buffer(q,bytes))return LIS_RCI_BADARG;
 int st=rs_init(S,RS_KRYLOV,1);if(st)return st;
 S->generations=*g;
 lis_rci_krylov_restore_request req={0};
 st=lis_rci_krylov_restore_start(S->child.k,p,o,ext,w,budget,v,x,r,stats,extra,&req);
 return rs_krylov_reply(S,st,&req,q);
}
int lis_rci_native_session_advance(lis_rci_native_session *S,uint64_t token,LIS_INT code,const lis_rci_native_generations *g,
 lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!rs_buffer(q,bytes))return LIS_RCI_BADARG;
 if(!S->pending)return LIS_RCI_STATE;
 if(token!=S->token)return LIS_RCI_PROTOCOL;
 if(S->kind==RS_KRYLOV&&(!g||code>INT_MAX||code<INT_MIN))return LIS_RCI_BADARG;
 if(!rs_current(S,g)){int st=rs_abort(S,LIS_RCI_NATIVE_STALE);*q=S->request;return st;}
 int st;S->pending=0;
 if(S->kind==RS_FEAST){lis_rci_native_feast_request r={0};st=lis_rci_native_feast_advance(S->child.e,S->child_token,code,&r);return rs_feast_reply(S,st,&r,q);}
 if(S->kind==RS_KRYLOV){lis_rci_krylov_restore_request r={0};st=lis_rci_krylov_restore_advance(S->child.k,S->child_token,(int)code,&r);return rs_krylov_reply(S,st,&r,q);}
 lis_rci_native_linear_request r={0};
 switch(S->kind){
 case RS_LINEAR:st=lis_rci_native_linear_restore_advance(S->child.l,S->child_token,code,&r);break;
 case RS_GCRO:st=lis_rci_native_gcro_restore_advance(S->child.g,S->child_token,code,&r);break;
 default:st=lis_rci_native_fgcro_restore_advance(S->child.f,S->child_token,code,&r);break;
 }
 return rs_linear_reply(S,st,&r,q);
}
int lis_rci_native_session_cancel(lis_rci_native_session *S,uint64_t token)
{
 if(!S)return LIS_RCI_BADARG;
 if(!S->pending)return LIS_RCI_STATE;
 if(token!=S->token)return LIS_RCI_PROTOCOL;
 return rs_cancel_child(S);
}
int lis_rci_native_session_pending(const lis_rci_native_session *S,lis_rci_native_session_request *q,size_t n)
{if(!S||!rs_buffer(q,n))return LIS_RCI_BADARG;*q=S->request;return S->pending?S->status:LIS_RCI_STATE;}
int lis_rci_native_session_status(const lis_rci_native_session *S,LIS_INT *native)
{if(!S||!native)return LIS_RCI_BADARG;*native=S->native;return S->status;}
