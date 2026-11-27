/* Diagnostic on public symmetric structural matrices, symmetric diagonal equilibration. */
#include <lis_rci_bridge.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef struct { size_t i,j; double v; } entry;
static void product(const entry *a,size_t nz,const double*x,long double*y,size_t n){
 for(size_t i=0;i<n;i++)y[i]=0;
 for(size_t k=0;k<nz;k++){size_t i=a[k].i,j=a[k].j;long double v=a[k].v;
  y[i]+=v*x[j];if(i!=j)y[j]+=v*x[i];}
}
/* Diagnostic use of internal LIS preconditioner API; not a new supported public bridge ABI. */
extern LIS_INT lis_precon_create(LIS_SOLVER,LIS_PRECON*);
extern LIS_INT lis_precon_destroy(LIS_PRECON);
extern LIS_INT lis_psolve_dispatch(LIS_SOLVER,LIS_VECTOR,LIS_VECTOR);
static int strong_solve(LIS_MATRIX A,const entry *a,size_t nz,const double *d,size_t n,const lis_rci_options *o,const double *b,double *y,lis_rci_info *info){
 LIS_SOLVER solver=NULL;LIS_PRECON pre=NULL;LIS_MATRIX P=NULL;LIS_VECTOR u=NULL,v=NULL;lis_rci_context *ctx=NULL;double *buf=NULL;
 int result=LIS_RCI_CALLBACK;
 *info=(lis_rci_info){.residual=NAN,.status=LIS_RCI_CALLBACK};
 /* A separate matrix prevents preconditioner setup from mutating the operator. */
 if(lis_matrix_create(LIS_COMM_WORLD,&P)||lis_matrix_set_size(P,0,(LIS_INT)n))goto done;
 for(size_t k=0;k<nz;k++){size_t i=a[k].i,j=a[k].j;double q=(a[k].v/d[i])/d[j];if(lis_matrix_set_value(LIS_ADD_VALUE,(LIS_INT)i,(LIS_INT)j,q,P)||(i!=j&&lis_matrix_set_value(LIS_ADD_VALUE,(LIS_INT)j,(LIS_INT)i,q,P)))goto done;}
 if(lis_matrix_set_type(P,LIS_MATRIX_CSR)||lis_matrix_assemble(P)||lis_solver_create(&solver))goto done;
 /* SSOR is SPD for the symmetric positive definite operator and omega=1.
    ILUT is used only for the nonsymmetric-capable RCI methods. */
 char cg_opts[]="-i cg -p ssor -ssor_omega 1.0 -scale none";
 char ilut_opts[]="-i bicgstab -p ilut -ilut_drop 1e-5 -ilut_rate 5 -scale none";
 const char*name=o->method==LIS_RCI_CG?"SSOR":"ILUT";
 if(lis_solver_set_option(o->method==LIS_RCI_CG?cg_opts:ilut_opts,solver)||lis_solver_set_matrix(P,solver))goto done;
 double setup=lis_wtime();LIS_PRECON tmp=NULL;LIS_INT e=lis_precon_create(solver,&tmp);
 if(e){printf("PRECON_SETUP method=%d precon=%s error=%lld\n",(int)o->method,name,(long long)e);goto done;}
 pre=tmp;solver->precon=pre;
 printf("PRECON_SETUP method=%d precon=%s seconds=%.6f\n",(int)o->method,name,lis_wtime()-setup);fflush(stdout);
 if(lis_vector_duplicate(A,&u)||lis_vector_duplicate(A,&v))goto done;
 ctx=lis_rci_create(o);buf=malloc(n*sizeof(*buf));if(!ctx||!buf)goto done;
 lis_rci_request q;int s=lis_rci_start(ctx,b,y,&q);
 while(s==LIS_RCI_MATVEC||s==LIS_RCI_PRECON){
  memcpy(buf,q.in,q.n*sizeof(*buf));
  e=lis_vector_set_values2(LIS_INS_VALUE,0,(LIS_INT)n,buf,u);
  if(!e){if(s==LIS_RCI_MATVEC)e=lis_matvec(A,u,v);else e=lis_psolve_dispatch(solver,u,v);}
  if(!e)e=lis_vector_get_values(v,0,(LIS_INT)n,q.out);
  s=lis_rci_advance(ctx,e?1:0,&q);
 }
 lis_rci_get_info(ctx,info);memcpy(y,lis_rci_solution(ctx),n*sizeof(*y));result=s;
done:
 free(buf);lis_rci_destroy(ctx);if(v)lis_vector_destroy(v);if(u)lis_vector_destroy(u);
 if(solver)solver->precon=NULL;if(pre)lis_precon_destroy(pre);if(solver)lis_solver_destroy(solver);if(P)lis_matrix_destroy(P);
 return result;
}

int main(int argc,char **argv){
 if(argc!=4){fprintf(stderr,"usage: %s matrix.mtx n triangle_entries\n",argv[0]);return 2;}
 char *endp=NULL;size_t n=(size_t)strtoull(argv[2],&endp,10);if(!endp||*endp||!n||n>20000)return 2;
 size_t nz=(size_t)strtoull(argv[3],&endp,10);if(!endp||*endp||!nz||nz>1000000)return 2;
 if(LIS_RCI_CHECK_CLIENT_ABI()!=LIS_RCI_DONE)return 2;
 int rc=2,init=0,failed=0;FILE *f=NULL;entry *a=NULL;double *v=NULL;long double *w=NULL;LIS_MATRIX A=NULL;
 #define CHECK(c) do{if(!(c)){fprintf(stderr,"INFRA_FAILURE line=%d expression=%s\n",__LINE__,#c);goto done;}}while(0)
 f=fopen(argv[1],"r");CHECK(f!=NULL);
 char line[1024];CHECK(fgets(line,sizeof(line),f)!=NULL);CHECK(strncmp(line,"%%MatrixMarket matrix coordinate real symmetric",46)==0);
 do{CHECK(fgets(line,sizeof(line),f)!=NULL);}while(line[0]=='%');
 size_t nr,nc,nn;CHECK(sscanf(line,"%zu%zu%zu",&nr,&nc,&nn)==3&&nr==n&&nc==n&&nn==nz);
 a=calloc(nz,sizeof(*a));v=calloc(5*n,sizeof(*v));w=calloc(2*n,sizeof(*w));CHECK(a&&v&&w);
 double *d=v,*exact=d+n,*b=exact+n,*y=b+n,*x=y+n;long double *orig=w,*ax=w+n;
 for(size_t k=0;k<nz;k++){size_t i,j;double q;CHECK(fscanf(f,"%zu%zu%lf",&i,&j,&q)==3&&i>=1&&i<=n&&j>=1&&j<=n&&isfinite(q));
  a[k]=(entry){i-1,j-1,q};if(i==j)d[i-1]+=q;}
 CHECK(fscanf(f,"%1023s",line)==EOF);fclose(f);f=NULL;
 for(size_t i=0;i<n;i++){CHECK(isfinite(d[i])&&d[i]>0);d[i]=sqrt(d[i]);exact[i]=sin(.017*(i+1))+.3*cos(.071*(i+3))+0.2;}
 product(a,nz,exact,orig,n);long double bn=0,xn=0;for(size_t i=0;i<n;i++){b[i]=(double)(orig[i]/d[i]);bn+=orig[i]*orig[i];xn+=(long double)exact[i]*exact[i];}
 CHECK(bn>0&&xn>0);
 int ac=1;char *avb[]={argv[0],NULL},**av=avb;CHECK(!lis_initialize(&ac,&av));init=1;
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,(LIS_INT)n));
 for(size_t k=0;k<nz;k++){size_t i=a[k].i,j=a[k].j;double q=(a[k].v/d[i])/d[j];CHECK(isfinite(q));CHECK(!lis_matrix_set_value(LIS_ADD_VALUE,(LIS_INT)i,(LIS_INT)j,q,A));if(i!=j)CHECK(!lis_matrix_set_value(LIS_ADD_VALUE,(LIS_INT)j,(LIS_INT)i,q,A));}
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 printf("STRONG_PRECON matrix=%s n=%zu triangle_nnz=%zu scaled=diagonal CG=SSOR others=ILUT rtol=1e-12 maxiter=2000 restart=50 forward_gate=1e-5 original_residual_gate=1e-10\n",argv[1],n,nz);
 for(int m=1;m<=4;m++){
  lis_rci_options o={n,2000,50,1e-12,0,(lis_rci_method)m};
  memset(y,0,n*sizeof(*y));
  lis_rci_info info={0};printf("START method=%d\n",m);fflush(stdout);
  double t=lis_wtime();int s=strong_solve(A,a,nz,d,n,&o,b,y,&info);double secs=lis_wtime()-t;
  double rr=NAN,er=NAN;if(s==LIS_RCI_DONE||s==LIS_RCI_MAXITER||s==LIS_RCI_BREAKDOWN){for(size_t i=0;i<n;i++)x[i]=y[i]/d[i];product(a,nz,x,ax,n);long double r=0,e=0;for(size_t i=0;i<n;i++){long double dr=orig[i]-ax[i],de=(long double)x[i]-exact[i];r+=dr*dr;e+=de*de;}rr=(double)sqrtl(r/bn);er=(double)sqrtl(e/xn);}
  int ok=s==LIS_RCI_DONE&&isfinite(rr)&&isfinite(er)&&rr<=1e-10&&er<=1e-5;
  printf("INDUSTRIAL_RESULT method=%d status=%d iter=%zu residual=%.9e forward_error=%.9e seconds=%.6f %s\n",m,s,info.iterations,rr,er,secs,ok?"PASS":"FAIL");fflush(stdout);
  failed+=!ok;
 }
 printf("INDUSTRIAL_SUMMARY %d/4 PASS\n",4-failed);rc=failed?1:0;
done:
 if(f)fclose(f);if(A)lis_matrix_destroy(A);if(init)lis_finalize();free(w);free(v);free(a);return rc;
}
