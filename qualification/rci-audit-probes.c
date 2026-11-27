#include "lis_rci.h"
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
static int run(void)
{
 lis_rci_options o={2,10,2,1,0,LIS_RCI_CG};
 lis_rci_context *c=lis_rci_create(&o);if(!c)return 1;
 double b[]={DBL_MAX,DBL_MAX};lis_rci_request q;
 int st=lis_rci_start(c,b,NULL,&q);
 if(st!=LIS_RCI_MATVEC)return 1;
 for(size_t i=0;i<2;i++)q.out[i]=q.in[i];
 st=lis_rci_advance(c,0,&q);lis_rci_info info;
 if(lis_rci_get_info(c,&info))return 1;
 printf("AUDIT RANGE LDBL_MAX_EXP=%d status=%d residual=%g isfinite=%d\n",LDBL_MAX_EXP,st,info.residual,!!isfinite(info.residual));
 lis_rci_destroy(c);
 o.rtol=1e-12;c=lis_rci_create(&o);if(!c)return 1;
 double bad[]={NAN,1};
 st=lis_rci_start(c,bad,NULL,&q);
 printf("AUDIT NAN_INPUT status=%d expected=%d\n",st,LIS_RCI_NONFINITE);
 if(st>0){lis_rci_get_info(c,&info);lis_rci_request_full f;if(lis_rci_get_pending(c,&f)>0)lis_rci_cancel(c,f.request.token);}
 double normal[]={1,2};st=lis_rci_start(c,normal,NULL,&q);if(st!=LIS_RCI_MATVEC)return 1;
 q.out[0]=0; /* leave q.out[1] at the sentinel */
 st=lis_rci_advance(c,0,&q);printf("AUDIT INCOMPLETE_OUTPUT status=%d expected=%d\n",st,LIS_RCI_NONFINITE);
 lis_rci_destroy(c);return 0;
}
int main(void){return run();}
