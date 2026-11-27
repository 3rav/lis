#include "lis_rci.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
void *__real_calloc(size_t,size_t);
void __real_free(void*);
static size_t calls,live,fail_at;
void *__wrap_calloc(size_t n,size_t z){calls++;if(calls==fail_at)return NULL;void*p=__real_calloc(n,z);if(p)live++;return p;}
void __wrap_free(void*p){if(p){if(!live)abort();live--;}__real_free(p);}
#define C(e) do{if(!(e)){fprintf(stderr,"OOM FAIL line=%d: %s\n",__LINE__,#e);return 1;}}while(0)
int main(void)
{
 for(int method=1;method<=4;method++){
  lis_rci_options o={16,100,8,1e-12,0,(lis_rci_method)method};
  for(size_t i=1;i<=2;i++){calls=0;fail_at=i;C(lis_rci_create(&o)==NULL&&live==0);}
  calls=0;fail_at=0;lis_rci_context*c=lis_rci_create(&o);C(c&&live==2&&calls==2);
  double b[16];for(int i=0;i<16;i++)b[i]=i+1;
  for(int rep=0;rep<3;rep++){
   lis_rci_request q;int st=lis_rci_start(c,b,NULL,&q);
   int limit=0;
   while(st>0&&++limit<1000){for(size_t i=0;i<q.n;i++)q.out[i]=(st==LIS_RCI_MATVEC?2:1)*q.in[i];st=lis_rci_advance(c,0,&q);}
   C(st==LIS_RCI_DONE&&calls==2&&live==2);
  }
  lis_rci_destroy(c);C(live==0);
 }
 puts("AUDIT CORE OOM 8 failure paths, 12 allocation-free standard solves PASS");return 0;
}
