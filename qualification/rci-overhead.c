#define _POSIX_C_SOURCE 200809L
#include "lis_rci.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
static double now(void){struct timespec t;if(clock_gettime(CLOCK_MONOTONIC,&t))exit(2);return t.tv_sec+1e-9*t.tv_nsec;}
static int run(size_t n,int method,int repeats)
{
 double *b=malloc(n*sizeof(double));if(!b)return 1;
 for(size_t i=0;i<n;i++)b[i]=1+(double)(i%17)/17;
 lis_rci_options opt={n,100,12,1e-10,0,(lis_rci_method)method};
 lis_rci_context *c=lis_rci_create(&opt);if(!c)return 1;
 size_t iterations=0,matvecs=0,precons=0;
 double start=now();
 for(int rep=-2;rep<repeats;rep++){
  if(rep==0)start=now();
  lis_rci_request q;int s=lis_rci_start(c,b,NULL,&q);
  while(s>0){
   if(s==LIS_RCI_MATVEC){for(size_t i=0;i<n;i++)q.out[i]=(4+(double)(i%13)/13)*q.in[i]-(i?q.in[i-1]:0)-(i+1<n?q.in[i+1]:0);}
   else if(s==LIS_RCI_PRECON){for(size_t i=0;i<n;i++)q.out[i]=q.in[i]/(4+(double)(i%13)/13);}
   else return 1;
   s=lis_rci_advance(c,0,&q);
  }
  if(s!=LIS_RCI_DONE)return 1;
  lis_rci_info info;lis_rci_get_info(c,&info);
  const double*x=lis_rci_solution(c);long double rn=0,bn=0;
  for(size_t i=0;i<n;i++){long double r=b[i]-((4+(double)(i%13)/13)*x[i]-(i?x[i-1]:0)-(i+1<n?x[i+1]:0));rn+=r*r;bn+=(long double)b[i]*b[i];}
  if(sqrtl(rn/bn)>1.05e-10L)return 1;
  if(rep>=0){iterations+=info.iterations;matvecs+=info.matvecs;precons+=info.precons;}
 }
 double elapsed=now()-start;
 printf("%zu,%d,%d,%zu,%zu,%zu,%.9f\n",n,method,repeats,iterations,matvecs,precons,elapsed);
 lis_rci_destroy(c);free(b);return 0;
}
int main(void){size_t sizes[]={16,256,4095,4096,8192};int reps[]={2000,500,40,40,20};for(int k=0;k<5;k++)for(int m=1;m<=4;m++)if(run(sizes[k],m,reps[k]))return 1;return 0;}
