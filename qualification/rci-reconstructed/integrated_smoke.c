#include "lis_rci_bridge.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define CHECK(c) do {if(!(c)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c);rc=1;goto cleanup;}}while(0)
int main(int argc,char **argv){
 LIS_MATRIX A=NULL; lis_rci_bridge *p=NULL; int rc=0; const LIS_INT n=64; double b[64],x[64];
 if(lis_initialize(&argc,&argv)!=0)return 1;
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,n));
 for(LIS_INT i=0;i<n;i++){CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));if(i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A));if(i+1<n)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A));}
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 for(int method=LIS_RCI_CG;method<=LIS_RCI_FGMRES;method++)for(int jac=0;jac<=1;jac++){
  lis_rci_options o={(size_t)n,200,12,1e-10,0,(lis_rci_method)method};CHECK(!lis_rci_bridge_create(A,&o,jac,&p));
  for(int rep=0;rep<2;rep++){
   for(LIS_INT i=0;i<n;i++){b[i]=(4-(i>0)-(i+1<n))*(rep+1.0);x[i]=0;}
   lis_rci_info info;CHECK(lis_rci_bridge_solve(p,b,x,&info)==LIS_RCI_DONE);
   for(LIS_INT i=0;i<n;i++)CHECK(fabs(x[i]-(rep+1.0))<1e-8);
   printf("method=%d jacobi=%d reuse=%d iterations=%zu PASS\n",method,jac,rep,info.iterations);
  }
  lis_rci_bridge_destroy(p);p=NULL;
 }
 cleanup:lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);lis_finalize();return rc;
}
