/*
 * Memory allocator failure-injection regression.
 * Uses GNU ld --wrap=malloc to fail payload and bookkeeping allocations.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "lis.h"

extern void *__real_malloc(size_t size);

static int fail_on=0;
static int calls=0;

void *__wrap_malloc(size_t size)
{
    if(fail_on && ++calls==fail_on) return NULL;
    return __real_malloc(size);
}

static void arm_failure(int nth)
{
    calls=0;
    fail_on=nth;
}

static void disarm(void)
{
    fail_on=0;
    calls=0;
}
#define CHECK(cond) do { \
    if(!(cond)) { \
        disarm(); \
        fprintf(stderr,"Allocator assertion failed: %s (line %d)\n", \
                #cond,__LINE__); \
        return 1; \
    } \
} while(0)

int main(void)
{
    unsigned char *p=NULL,*q=NULL,*z=NULL,*raw=NULL;
    size_t i;
    void *out;

    /* Both payload and metadata failures must return NULL cleanly. */
    for(i=1;i<=2;i++) {
        arm_failure((int)i);
        out=lis_malloc(64,"allocfail");
        disarm();
        CHECK(out==NULL);
        arm_failure((int)i);
        out=lis_calloc(64,"callocfail");
        disarm();
        CHECK(out==NULL);
    }

    CHECK(lis_malloc((size_t)-1,"overflow")==NULL);
    CHECK(lis_calloc((size_t)-1,"overflow")==NULL);

    p=(unsigned char *)lis_malloc(64,"realloc");
    CHECK(p!=NULL && lis_is_malloc(p)==LIS_TRUE);
    for(i=0;i<64;i++) p[i]=(unsigned char)(i+1);

    arm_failure(1);
    out=lis_realloc(p,128);
    disarm();
    CHECK(out==NULL);
    CHECK(lis_is_malloc(p)==LIS_TRUE);
    for(i=0;i<64;i++) CHECK(p[i]==(unsigned char)(i+1));
    CHECK(lis_realloc(p,(size_t)-1)==NULL);
    CHECK(lis_is_malloc(p)==LIS_TRUE);

    /* Catch copy-of-old-size bugs when shrinking 64 bytes to 8 bytes. */
    q=(unsigned char *)lis_realloc(p,8);
    CHECK(q!=NULL && lis_is_malloc(q)==LIS_TRUE);
    for(i=0;i<8;i++) CHECK(q[i]==(unsigned char)(i+1));
    p=NULL; /* Realloc invalidated the original pointer. */

    p=(unsigned char *)lis_realloc(q,128);
    CHECK(p!=NULL && lis_is_malloc(p)==LIS_TRUE);
    for(i=0;i<8;i++) CHECK(p[i]==(unsigned char)(i+1));
    q=NULL;

    z=(unsigned char *)lis_calloc(80,"zeros");
    CHECK(z!=NULL && lis_is_malloc(z)==LIS_TRUE);
    for(i=0;i<80;i++) CHECK(z[i]==0);

    /* Untracked pointers keep the system realloc semantics. */
    raw=(unsigned char *)malloc(16);
    CHECK(raw!=NULL);
    raw[0]=99;
    q=(unsigned char *)lis_realloc(raw,32);
    CHECK(q!=NULL && q[0]==99);
    raw=NULL;
    lis_free(q);
    q=NULL;

    lis_free(z);
    z=NULL;
    lis_free(p);
    p=NULL;

    puts("MEMORY_ALLOCATOR_FAILURES_PASSED");
    return 0;
}
