/* Full-lifecycle experimental native driver, included in lis_solver.c. */
#include "lis_rci_native_full_internal.h"
#if defined(USE_MPI)||defined(USE_QUAD_PRECISION)||defined(_COMPLEX)||defined(_LONG__DOUBLE)
#error "full lifecycle qualification is real double non-MPI"
#endif
#include "full_linear_continuation.inc"
static int full_buffer(lis_rci_native_session_request *q,size_t n){return q&&n>=sizeof(*q);}
static int full_threads(void){
#ifdef _OPENMP
 return omp_get_max_threads();
#else
 return 1;
#endif
}
static lis_rci_native_generations full_generations(LIS_MATRIX A,LIS_MATRIX B)
{
 lis_rci_native_generations g={A?A->structure_generation:0,A?A->value_generation:0,B?B->structure_generation:0,B?B->value_generation:0};return g;
}
static int full_current(const lis_rci_native_full *S)
{
 lis_rci_native_generations g=full_generations(S->A,S->B),h=S->generations;
 return g.a_structure==h.a_structure&&g.a_values==h.a_values&&g.b_structure==h.b_structure&&g.b_values==h.b_values;
}
static int full_map(LIS_INT rc)
{
 return rc==LIS_SUCCESS?LIS_RCI_DONE:rc==LIS_MAXITER?LIS_RCI_MAXITER:rc==LIS_BREAKDOWN?LIS_RCI_BREAKDOWN:rc==LIS_ERR_ILL_ARG?LIS_RCI_BADARG:LIS_RCI_CALLBACK;
}
int full_native_finish(lis_rci_native_full *S,LIS_INT rc)
{
 S->api_code=rc;S->native_code=rc?rc:(S->mode==FULL_LINEAR?S->solver->retcode:S->esolver->retcode);
 if(S->forced){
  S->native_code=S->forced_code;
  if(S->solver){
   S->solver->retcode=S->forced_code;S->solver->resid=NAN;
   if(S->solver->options[LIS_OPTIONS_SOLVER]==LIS_SOLVER_GCRODR)lis_gcrodr_clear_recycle(S->solver);
   if(S->solver->options[LIS_OPTIONS_SOLVER]==LIS_SOLVER_FGCRODR)lis_fgcrodr_clear_recycle(S->solver);
  }else{S->esolver->retcode=S->forced_code;S->esolver->nev=0;*S->evalue=NAN;}
 }
 S->stage=0;S->pending=0;memset(&S->request,0,sizeof S->request);
 return S->status=S->forced?S->forced:full_map(S->native_code);
}
static void full_force(lis_rci_native_full *S,int status,LIS_INT code){if(!S->forced){S->forced=status;S->forced_code=code;}}
static int full_child_result(lis_rci_native_full *S,int st,lis_rci_native_session_request *q)
{
 if(st<=0){
  lis_rci_native_session_status(S->child,&S->completion);
  if(st==LIS_RCI_NONFINITE||st==LIS_RCI_CANCELLED||st==LIS_RCI_PROTOCOL||st==LIS_RCI_NATIVE_STALE)
   full_force(S,st,LIS_BREAKDOWN);
  if(S->forced)S->completion=S->forced_code;
  return FULL_AGAIN;
 }
 S->child_token=q->token;
 if(S->token==UINT64_MAX||!full_current(S)){
  full_force(S,S->token==UINT64_MAX?LIS_RCI_PROTOCOL:LIS_RCI_NATIVE_STALE,LIS_BREAKDOWN);
  (void)lis_rci_native_session_cancel(S->child,S->child_token);
  S->completion=S->forced_code;return FULL_AGAIN;
 }
 S->request=*q;S->request.token=++S->token;S->request.generations=S->generations;
 S->pending=FULL_CHILD;return S->status=st;
}
int full_native_begin_linear(lis_rci_native_full *S,LIS_SOLVER solver,LIS_INT method)
{
 if(method==LIS_SOLVER_FGCRODR){
  LIS_INT k,update;LIS_REAL keep,reject;
  LIS_INT e=lis_solver_get_fgcrodr_recycle_config(solver,&k,&keep,&reject,&update);
  if(e){S->completion=e;return FULL_AGAIN;}
  if(k==0)method=LIS_SOLVER_FGMRES;
 }
 lis_rci_native_session_request q={0};
 int st=lis_rci_native_session_start_linear(S->child,solver,method,&q,sizeof q);
 return full_child_result(S,st,&q);
}
int full_native_begin_feast(lis_rci_native_full *S,LIS_ESOLVER e)
{
 lis_rci_native_session_request q={0};
 int st=lis_rci_native_session_start_feast(S->child,e,&q,sizeof q);
 return full_child_result(S,st,&q);
}
int full_native_cert(lis_rci_native_full *S,LIS_MATRIX A,LIS_VECTOR in,LIS_VECTOR out)
{
 if(S->token==UINT64_MAX||!full_current(S)){
  full_force(S,S->token==UINT64_MAX?LIS_RCI_PROTOCOL:LIS_RCI_NATIVE_STALE,LIS_BREAKDOWN);S->completion=S->forced_code;return FULL_AGAIN;
 }
 for(LIS_INT i=0;i<in->n;i++)if(!isfinite(in->value[i])){full_force(S,LIS_RCI_NONFINITE,LIS_BREAKDOWN);S->completion=S->forced_code;return FULL_AGAIN;}
 for(LIS_INT i=0;i<out->n;i++)out->value[i]=NAN;
 S->cert_in[0]=in;S->cert_out[0]=out;
 S->request=(lis_rci_native_session_request){.version=LIS_RCI_NATIVE_SESSION_VERSION,.size=sizeof(S->request),
 .operation=LIS_RCI_OP_A,.flags=LIS_RCI_NATIVE_VECTORS,.index_bytes=sizeof(LIS_INT),.scalar_bytes=sizeof(LIS_SCALAR),.real_bytes=sizeof(LIS_REAL),
 .token=++S->token,.rows=in->n,.columns=1,.generations=S->generations,.in=S->cert_in,.out=S->cert_out,.A=A,.solver=S->solver};
 S->pending=FULL_CERT;return S->status=LIS_RCI_OP_A;
}
static int full_drive(lis_rci_native_full *S,lis_rci_native_session_request *q)
{
 int st=FULL_AGAIN;
 while(st==FULL_AGAIN){
  if(S->mode==FULL_FEAST)st=full_native_feast_step(S);
  else st=S->stage==1?full_step_outer(S):full_step_kernel(S);
 }
 if(q)*q=S->request;return st;
}
static int full_reset(lis_rci_native_full *S,int mode,LIS_MATRIX A,LIS_MATRIX B)
{
 if(!S)return LIS_RCI_BADARG;
 if(S->stage||S->pending)return LIS_RCI_STATE;
 free(S->frame);S->frame=NULL;S->mode=mode;S->A=A;S->B=B;S->solver=NULL;S->esolver=NULL;S->evalue=NULL;
 S->forced=0;S->forced_code=0;S->api_code=0;S->native_code=0;S->completion=0;S->status=LIS_RCI_STATE;
 S->threads=full_threads();S->generations=full_generations(A,B);memset(&S->request,0,sizeof S->request);return LIS_RCI_DONE;
}
lis_rci_native_full *lis_rci_native_full_create(void)
{
 lis_rci_native_full *S=calloc(1,sizeof(*S));if(!S)return NULL;
 S->child=lis_rci_native_session_create();if(!S->child){free(S);return NULL;}S->status=LIS_RCI_STATE;return S;
}
void lis_rci_native_full_destroy(lis_rci_native_full *S)
{
 if(S){if(S->pending)(void)lis_rci_native_full_cancel(S,S->token);lis_rci_native_session_destroy(S->child);free(S->frame);free(S);}
}
int lis_rci_native_full_start_linear(lis_rci_native_full *S,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_SOLVER solver,lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!A||!b||!x||!solver||!full_buffer(q,bytes)||b==x)return LIS_RCI_BADARG;
 if(S->stage||S->pending)return LIS_RCI_STATE;
 LIS_INT m=solver->options[LIS_OPTIONS_SOLVER],p=solver->options[LIS_OPTIONS_PRECON];
 if(m!=LIS_SOLVER_CG&&m!=LIS_SOLVER_BICGSTAB&&m!=LIS_SOLVER_GMRES&&m!=LIS_SOLVER_FGMRES&&m!=LIS_SOLVER_GCRODR&&m!=LIS_SOLVER_FGCRODR)return LIS_RCI_BADARG;
 if((A->matrix_type!=LIS_MATRIX_CSR&&A->matrix_type!=LIS_MATRIX_USER)||(p!=LIS_PRECON_TYPE_NONE&&p!=LIS_PRECON_TYPE_JACOBI)||
  (A->matrix_type==LIS_MATRIX_USER&&p!=LIS_PRECON_TYPE_NONE)||solver->options[LIS_OPTIONS_SCALE]!=LIS_SCALE_NONE||
  solver->options[LIS_OPTIONS_STORAGE]||solver->options[LIS_OPTIONS_USE_AT]||solver->options[LIS_OPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE||
  solver->nullspace_dim||solver->left_nullspace_dim||solver->near_nullspace_dim||solver->num_threads||solver->setup)return LIS_RCI_BADARG;
 int st=full_reset(S,FULL_LINEAR,A,NULL);if(st)return st;
 S->solver=solver;struct full_frames_linear *C=calloc(1,sizeof(*C));if(!C)return LIS_RCI_BADARG;
 C->outer.A=A;C->outer.b=b;C->outer.x=x;C->outer.solver=solver;S->frame=C;S->stage=1;
 return full_drive(S,q);
}
int lis_rci_native_full_start_feast(lis_rci_native_full *S,LIS_MATRIX A,LIS_MATRIX B,LIS_VECTOR x,LIS_SCALAR *value,LIS_ESOLVER e,lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!A||!x||!value||!e||!full_buffer(q,bytes))return LIS_RCI_BADARG;
 if(S->stage||S->pending)return LIS_RCI_STATE;
 if(A->matrix_type!=LIS_MATRIX_CSR||(B&&B->matrix_type!=LIS_MATRIX_CSR)||e->options[LIS_EOPTIONS_ESOLVER]!=LIS_ESOLVER_FEAST||
  e->options[LIS_EOPTIONS_STORAGE]||e->options[LIS_EOPTIONS_MODE]||e->options[LIS_EOPTIONS_PRECISION]!=LIS_PRECISION_DOUBLE||
  e->subspace_retry_max||e->transform!=LIS_EIGEN_TRANSFORM_AUTO)return LIS_RCI_BADARG;
 if(e->requested_nev_active&&(e->requested_nev<0||e->requested_nev>e->options[LIS_EOPTIONS_SUBSPACE]))return LIS_RCI_BADARG;
 int st=full_reset(S,FULL_FEAST,A,B);if(st)return st;
 S->esolver=e;S->evalue=value;e->subspace_retry_attempts=0;e->subspace_retry_final_size=e->options[LIS_EOPTIONS_SUBSPACE];
 if(e->requested_nev_active)e->nev=e->requested_nev;
 if(full_native_feast_frame(S,A,B,x,value,e))return LIS_RCI_BADARG;
 return full_drive(S,q);
}
static int full_abort(lis_rci_native_full *S,int status,lis_rci_native_session_request *q)
{
 full_force(S,status,LIS_BREAKDOWN);
 if(S->pending==FULL_CHILD)(void)lis_rci_native_session_cancel(S->child,S->child_token);
 S->pending=0;S->completion=S->forced_code;memset(&S->request,0,sizeof S->request);
 return full_drive(S,q);
}
int lis_rci_native_full_advance(lis_rci_native_full *S,uint64_t token,LIS_INT code,lis_rci_native_session_request *q,size_t bytes)
{
 if(!S||!full_buffer(q,bytes))return LIS_RCI_BADARG;
 if(!S->pending)return LIS_RCI_STATE;
 if(token!=S->token)return LIS_RCI_PROTOCOL;
 if(full_threads()!=S->threads)return LIS_RCI_STATE;
 if(!full_current(S))return full_abort(S,LIS_RCI_NATIVE_STALE,q);
 if(code)full_force(S,LIS_RCI_CALLBACK,code);
 int pending=S->pending;S->pending=0;
 if(pending==FULL_CHILD){
  lis_rci_native_session_request child={0};int st=lis_rci_native_session_advance(S->child,S->child_token,code,NULL,&child,sizeof child);
  st=full_child_result(S,st,&child);
  if(st!=FULL_AGAIN){*q=S->request;return st;}
 }else{
  if(!code){LIS_VECTOR out=(LIS_VECTOR)S->cert_out[0];for(LIS_INT i=0;i<out->n;i++)if(!isfinite(out->value[i])){full_force(S,LIS_RCI_NONFINITE,LIS_BREAKDOWN);break;}}
  S->completion=S->forced?S->forced_code:code;
 }
 memset(&S->request,0,sizeof S->request);return full_drive(S,q);
}
int lis_rci_native_full_cancel(lis_rci_native_full *S,uint64_t token)
{
 if(!S)return LIS_RCI_BADARG;if(!S->pending)return LIS_RCI_STATE;if(token!=S->token)return LIS_RCI_PROTOCOL;
 return full_abort(S,LIS_RCI_CANCELLED,NULL);
}
int lis_rci_native_full_status(const lis_rci_native_full *S,LIS_INT *api,LIS_INT *native)
{if(!S||!api||!native)return LIS_RCI_BADARG;*api=S->api_code;*native=S->native_code;return S->status;}
