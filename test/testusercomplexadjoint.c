#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H) && !defined(USE_MPI)
#include <complex.h>
#define N 3
typedef struct { LIS_SCALAR a[N][N]; LIS_INT bad_row; } CTX;
static LIS_INT row(void *raw,LIS_INT i,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals)
{
    CTX *c=(CTX *)raw;
    LIS_INT k;
    if(!nnz || i<0 || i>=N)return LIS_ERR_ILL_ARG;
    if(i==c->bad_row)return LIS_ERR_ILL_ARG;
    if(!cols && !vals){*nnz=N;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<N)return LIS_ERR_ILL_ARG;
    for(k=0;k<N;k++){cols[k]=k;vals[k]=c->a[i][k];}
    *nnz=N;
    return LIS_SUCCESS;
}
static LIS_INT ax(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)raw;
    LIS_INT i,j;
    for(i=0;i<N;i++){
        y[i]=0;
        for(j=0;j<N;j++)y[i]+=c->a[i][j]*x[j];
    }
    return LIS_SUCCESS;
}
static LIS_INT ahx(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)raw;
    LIS_INT i,j;
    for(i=0;i<N;i++){
        y[i]=0;
        for(j=0;j<N;j++)y[i]+=conj(c->a[j][i])*x[j];
    }
    return LIS_SUCCESS;
}
static int close_scalar(LIS_SCALAR a,LIS_SCALAR b)
{
    return cabs(a-b)<=1e-11*(1+cabs(b));
}
static int check(LIS_MATRIX A,LIS_VECTOR x,LIS_VECTOR z,LIS_VECTOR y,LIS_VECTOR h)
{
    LIS_SCALAR ex[N]={0},eh[N]={0},lhs=0,rhs=0;
    LIS_INT i,j;
    /* Reference uses literal asymmetric coefficients, not USER callbacks. */
    const LIS_SCALAR m[N][N]={
        {4+I, 1+2*I, -I},
        {2-I, 5-2*I, 3+I},
        {1+I, -2*I, 6+3*I}};
    if(lis_matvec(A,x,y) || lis_matvech(A,z,h))return 1;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            ex[i]+=m[i][j]*x->value[j];
            eh[i]+=conj(m[j][i])*z->value[j];
        }
        if(!close_scalar(y->value[i],ex[i]) ||
           !close_scalar(h->value[i],eh[i]))return 1;
        lhs+=conj(z->value[i])*y->value[i];
        rhs+=conj(h->value[i])*x->value[i];
    }
    return !close_scalar(lhs,rhs);
}
static int exercise(void)
{
    CTX c={{{4+I,1+2*I,-I},{2-I,5-2*I,3+I},{1+I,-2*I,6+3*I}},-1};
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,z=NULL,y=NULL,h=NULL;
    LIS_INT i,threadsafe;
    int failed=1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,N) ||
       lis_matrix_set_user(A,&c,ax,ahx) ||
       lis_vector_duplicate(A,&x) || lis_vector_duplicate(A,&z) ||
       lis_vector_duplicate(A,&y) || lis_vector_duplicate(A,&h))goto done;
    for(i=0;i<N;i++){
        x->value[i]=(1+i)+(2-i)*I;
        z->value[i]=(3-i)-(1+i)*I;
    }
    if(check(A,x,z,y,h))goto done;
    /* USER descriptors are initialized once; create a separate row-only matrix. */
    lis_matrix_destroy(A); A=NULL;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,N) ||
       lis_matrix_set_user(A,&c,NULL,NULL) || lis_matrix_set_user_row(A,row))goto done;
    for(threadsafe=0;threadsafe<=1;threadsafe++){
        if(lis_matrix_set_user_row_thread_safe(A,threadsafe) ||
           check(A,x,z,y,h))goto done;
        c.bad_row=1;
        if(lis_matvec(A,x,y)!=LIS_ERR_ILL_ARG ||
           lis_matvech(A,z,h)!=LIS_ERR_ILL_ARG)goto done;
        c.bad_row=-1;
        if(check(A,x,z,y,h))goto done;
    }
    failed=0;
done:
    if(h)lis_vector_destroy(h);
    if(y)lis_vector_destroy(y);
    if(z)lis_vector_destroy(z);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}
#endif
int main(int argc,char **argv)
{
#if defined(_COMPLEX) && defined(HAVE_COMPLEX_H) && !defined(USE_MPI)
    int failed;
    if(lis_initialize(&argc,&argv))return 1;
    failed=exercise();
    lis_finalize();
    printf("USER_COMPLEX_ADJOINT %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#else
    (void)argc;(void)argv;
    puts("USER_COMPLEX_ADJOINT SKIP: requires serial complex build");
    return 77;
#endif
}
