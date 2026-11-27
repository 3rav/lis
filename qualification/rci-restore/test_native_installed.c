/* Compile exclusively against installed public headers. */
#define LIS_USE_INSTALLED_CONFIG
#include <lis_rci_native.h>
#include <stdio.h>
#include <math.h>
#define C(x) do{if(!(x)){fprintf(stderr,"INSTALLED FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
static lis_rci_native_context *create(int kind){
 LIS_ABI_INFO a;lis_abi_get_client_info(&a);
 return lis_rci_native_create(kind,LIS_RCI_NATIVE_API_VERSION,sizeof(lis_rci_native_session_request),&a,sizeof a);
}
static int drive(lis_rci_native_context *s,int st,lis_rci_native_session_request q)
{
 uint64_t last=0;
 while(st>0){
  lis_rci_native_session_request p;
  C(q.token>last&&lis_rci_native_pending(s,&p)==st&&p.token==q.token);last=q.token;
  C(lis_rci_native_advance(s,last-1,0,&q)==LIS_RCI_PROTOCOL);
  LIS_INT code=lis_rci_native_execute_default(s,last);
  st=lis_rci_native_advance(s,last,code,&q);
 }
 LIS_INT native;C(lis_rci_native_status(s,&native)==LIS_RCI_DONE&&st==LIS_RCI_DONE&&native==LIS_SUCCESS);
 return 0;
}
static int diag(void *u,int64_t n,const double *x,double *y){(void)u;for(int64_t i=0;i<n;i++)y[i]=(i+1)*x[i];return 0;}
static int ident(void *u,int64_t n,const double *x,double *y){(void)u;for(int64_t i=0;i<n;i++)y[i]=x[i];return 0;}
static int inverse(void *u,int64_t n,const double *x,double *y){(void)u;for(int64_t i=0;i<n;i++)y[i]=x[i]/(i+.5);return 0;}
int main(int argc,char **argv)
{
 C(!lis_initialize(&argc,&argv));LIS_ABI_INFO abi;lis_abi_get_client_info(&abi);
 C(!lis_rci_native_create(LIS_NR_LINEAR,999,sizeof(lis_rci_native_session_request),&abi,sizeof abi));
 C(!lis_rci_native_create(LIS_NR_LINEAR,1,1,&abi,sizeof abi));
 LIS_MATRIX A;LIS_VECTOR b,x;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,12));
 for(int i=0;i<12;i++)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,i+1,A));
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_vector_set_all(1,b));
 const char *m[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int i=0;i<6;i++){
  LIS_SOLVER solver;C(!lis_solver_create(&solver));char opts[200];
  snprintf(opts,sizeof opts,"-i %s -p none -restart 8 -gcrodr_keep 2 -fgcrodr_k 2 -tol 1e-10 -maxiter 300",m[i]);
  C(!lis_solver_set_option(opts,solver));C(!lis_vector_set_all(0,x));
  lis_rci_native_context *s=create(LIS_NR_LINEAR);C(s);
  lis_rci_native_session_request q;int st=lis_rci_native_start_linear(s,A,b,x,solver,&q);C(st>0);C(!drive(s,st,q));
  LIS_SCALAR v[12];C(!lis_vector_get_values(x,0,12,v));for(int j=0;j<12;j++)C(fabs(v[j]-1.0/(j+1))<1e-8);
  lis_rci_native_destroy(s);lis_solver_destroy(solver);printf("INSTALLED PASS %s\n",m[i]);
 }
 LIS_ESOLVER e;C(!lis_esolver_create(&e));C(!lis_esolver_set_option("-e feast -ss 6 -emaxiter 20 -etol 1e-7 -feast_emin 1.5 -feast_emax 3.5 -feast_points 16",e));C(!lis_vector_set_all(1,x));
 lis_rci_native_context *s=create(LIS_NR_FEAST);C(s);LIS_SCALAR ev;lis_rci_native_session_request q;
 int st=lis_rci_native_start_feast(s,A,NULL,x,&ev,e,&q);C(st>0);C(!drive(s,st,q));C(fabs(ev-2)<1e-6||fabs(ev-3)<1e-6);
 lis_rci_native_destroy(s);lis_esolver_destroy(e);puts("INSTALLED PASS FEAST");
 s=create(LIS_NR_KRYLOV);C(s);lis_krylov_problem64 p={12,LIS_KRYLOV_SHIFT_INVERT_B,.5,NULL,diag,ident,inverse};
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=2;o.ncv=6;o.inner_tol=1e-12;o.maxiter=100;
 double v[2],vec[24],res[2];lis_krylov_stats stats;
 st=lis_rci_native_start_krylov(s,&p,&o,NULL,NULL,NULL,v,vec,res,&stats,NULL,&q);C(st>0);
 while(st>0){C(q.representation==LIS_NR_RAW64);int rc=q.operation==LIS_NR_OP_A?diag(NULL,q.n64,q.raw_in,q.raw_out):q.operation==LIS_NR_OP_B?ident(NULL,q.n64,q.raw_in,q.raw_out):inverse(NULL,q.n64,q.raw_in,q.raw_out);st=lis_rci_native_advance(s,q.token,rc,&q);}
 C(st==LIS_RCI_DONE&&stats.nreturned==2);lis_rci_native_destroy(s);puts("INSTALLED PASS KRYLOV");
 lis_vector_destroy(b);lis_vector_destroy(x);lis_matrix_destroy(A);lis_finalize();return 0;
}
