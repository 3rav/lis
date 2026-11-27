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
int main(int argc,char **argv){
 if(argc!=4){fprintf(stderr,"usage: %s matrix.mtx n triangle_entries\n",argv[0]);return 2;}
 char *endp=NULL;size_t n=(size_t)strtoull(argv[2],&endp,10);if(!endp||*endp||!n||n>20000)return 2;
 size_t nz=(size_t)strtoull(argv[3],&endp,10);if(!endp||*endp||!nz||nz>1000000)return 2;
 if(LIS_RCI_CHECK_CLIENT_ABI()!=LIS_RCI_DONE)return 2;
 int rc=2,init=0,failed=0;FILE *f=NULL;entry *a=NULL;double *v=NULL;long double *w=NULL;LIS_MATRIX A=NULL;lis_rci_bridge *p=NULL;
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
 printf("INDUSTRIAL matrix=%s n=%zu triangle_nnz=%zu scaled=diagonal rtol=1e-12 maxiter=2000 restart=50 forward_gate=1e-5 original_residual_gate=1e-10\n",argv[1],n,nz);
 for(int m=1;m<=4;m++){
  lis_rci_options o={n,2000,50,1e-12,0,(lis_rci_method)m};
  CHECK(!lis_rci_bridge_create(A,&o,1,&p));memset(y,0,n*sizeof(*y));
  lis_rci_info info={0};printf("START method=%d\n",m);fflush(stdout);
  double t=lis_wtime();int s=lis_rci_bridge_solve(p,b,y,&info);double secs=lis_wtime()-t;
  double rr=NAN,er=NAN;if(s==LIS_RCI_DONE){for(size_t i=0;i<n;i++)x[i]=y[i]/d[i];product(a,nz,x,ax,n);long double r=0,e=0;for(size_t i=0;i<n;i++){long double dr=orig[i]-ax[i],de=(long double)x[i]-exact[i];r+=dr*dr;e+=de*de;}rr=(double)sqrtl(r/bn);er=(double)sqrtl(e/xn);}
  int ok=s==LIS_RCI_DONE&&isfinite(rr)&&isfinite(er)&&rr<=1e-10&&er<=1e-5;
  printf("INDUSTRIAL_RESULT method=%d status=%d iter=%zu residual=%.9e forward_error=%.9e seconds=%.6f %s\n",m,s,info.iterations,rr,er,secs,ok?"PASS":"FAIL");fflush(stdout);
  failed+=!ok;lis_rci_bridge_destroy(p);p=NULL;
 }
 printf("INDUSTRIAL_SUMMARY %d/4 PASS\n",4-failed);rc=failed?1:0;
done:
 if(f)fclose(f);lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);if(init)lis_finalize();free(w);free(v);free(a);return rc;
}
