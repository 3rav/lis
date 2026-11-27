/* A small real structural matrix (HB/bcsstk01); native LIS matvec through RCI bridge.
 * Independently compute b and residual from the validated Matrix Market entries. */
#include "lis_rci_bridge.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#define N 48
#define CHECK(c) do {if(!(c)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c);rc=1;goto cleanup;}}while(0)
static double a[N][N];
static void product(const double*x,double*y){for(int i=0;i<N;i++){long double t=0;for(int j=0;j<N;j++)t+=(long double)a[i][j]*x[j];y[i]=(double)t;}}
int main(int argc,char**argv){
 if(argc!=2){fprintf(stderr,"usage: %s bcsstk01.mtx\n",argv[0]);return 2;}
 FILE*f=fopen(argv[1],"r");if(!f)return 2;
 LIS_MATRIX A=NULL;lis_rci_bridge*p=NULL;int rc=0,initialized=0;char line[128];int nr,nc,nz;
 CHECK(fgets(line,sizeof line,f)&&strcmp(line,"%%MatrixMarket matrix coordinate real symmetric\n")==0);
 CHECK(fscanf(f,"%d%d%d",&nr,&nc,&nz)==3&&nr==N&&nc==N&&nz==224);
 for(int k=0;k<nz;k++){int i,j;double v;CHECK(fscanf(f,"%d%d%lf",&i,&j,&v)==3&&i>=1&&i<=N&&j>=1&&j<=i&&isfinite(v)&&v!=0);CHECK(a[i-1][j-1]==0);a[i-1][j-1]=a[j-1][i-1]=v;}
 fclose(f);f=NULL;
 /* Do not send the matrix filename into LIS's option parser. */
 int ac=1;char *av_buf[]={argv[0],NULL};char **av=av_buf;
 CHECK(!lis_initialize(&ac,&av));initialized=1;
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,N));
 int entries=0;for(int i=0;i<N;i++)for(int j=0;j<N;j++)if(a[i][j]!=0){entries++;CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,j,a[i][j],A));}CHECK(entries==400);
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 for(int m=1;m<=4;m++){
  lis_rci_options o={N,1500,30,1e-12,0,(lis_rci_method)m};CHECK(!lis_rci_bridge_create(A,&o,1,&p));
  for(int rep=0;rep<2;rep++){
   double exact[N],b[N],x[N]={0},ax[N];for(int i=0;i<N;i++)exact[i]=rep?sin(.07*i)+1:sin(.37*(i+1))+.3*cos(.19*(i+3));
   product(exact,b);lis_rci_info info;int s=lis_rci_bridge_solve(p,b,x,&info);product(x,ax);
   long double e=0,d=0,r=0,bn=0;for(int i=0;i<N;i++){long double t=(long double)x[i]-exact[i],q=(long double)b[i]-ax[i];e+=t*t;d+=(long double)exact[i]*exact[i];r+=q*q;bn+=(long double)b[i]*b[i];}
   double er=(double)sqrtl(e/d),rr=(double)sqrtl(r/bn);
   int ok=s==LIS_RCI_DONE&&rr<=1.05e-12&&er<=5e-7&&isfinite(er)&&isfinite(rr);
   printf("BCSSTK01 method=%d rhs=%d status=%d iter=%zu residual=%.9e error=%.9e %s\n",m,rep,s,info.iterations,rr,er,ok?"PASS":"FAIL");CHECK(ok);
  }
  lis_rci_bridge_destroy(p);p=NULL;
 }
cleanup:
 if(f)fclose(f);lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);if(initialized)lis_finalize();return rc;
}
