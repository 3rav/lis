/*
 * Assembly benchmark: three equivalent 5-point Laplacian matrices.
 * Time creation + assembly only, then verify SpMV numerically.
 * No timing-based assertions: CI hardware is not a benchmark reference.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"

static const char *names[]={"rco","rco_preallocated","csr"};

static LIS_INT insert_entry(LIS_MATRIX A,LIS_INT mode,LIS_INT row,
                            LIS_INT col,LIS_SCALAR v,LIS_INT *index,
                            LIS_SCALAR *value,LIS_INT *count)
{
    if(mode==2) {
        index[*count]=col;
        value[*count]=v;
        ++*count;
        return LIS_SUCCESS;
    }
    return lis_matrix_set_value(LIS_INS_VALUE,row,col,v,A);
}

static LIS_INT build(LIS_INT side,LIS_INT mode,LIS_MATRIX *out,double *time_s)
{
    LIS_MATRIX A=NULL;
    LIS_INT err,n,nnz,row,x,y,count=0;
    LIS_INT *ptr=NULL,*index=NULL;
    LIS_SCALAR *value=NULL;
    double t0=lis_wtime();

    *out=NULL;
    n=side*side;
    nnz=5*n-4*side;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) return err;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto fail;

    if(mode==1) {
        err=lis_matrix_malloc(A,5,NULL);
        if(err) goto fail;
    } else if(mode==2) {
        err=lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
        if(err) goto fail;
    }
    for(row=0;row<n;row++) {
        x=row%side; y=row/side;
        if(mode==2) ptr[row]=count;
        err=insert_entry(A,mode,row,row,4.0,index,value,&count);
        if(err) goto fail;
        if(x>0) {
            err=insert_entry(A,mode,row,row-1,-1.0,index,value,&count);
            if(err) goto fail;
        }
        if(x+1<side) {
            err=insert_entry(A,mode,row,row+1,-1.0,index,value,&count);
            if(err) goto fail;
        }
        if(y>0) {
            err=insert_entry(A,mode,row,row-side,-1.0,index,value,&count);
            if(err) goto fail;
        }
        if(y+1<side) {
            err=insert_entry(A,mode,row,row+side,-1.0,index,value,&count);
            if(err) goto fail;
        }
        if(mode==2) ptr[row+1]=count;
    }
    if(mode==2) {
        if(count!=nnz) { err=LIS_ERR_ILL_ARG; goto fail; }
        err=lis_matrix_set_csr(nnz,ptr,index,value,A);
        if(err) goto fail;
        ptr=NULL; index=NULL; value=NULL; /* LIS now owns the arrays. */
    }
    err=lis_matrix_assemble(A);
    if(err) goto fail;
    *time_s=lis_wtime()-t0;
    *out=A;
    return LIS_SUCCESS;

fail:
    if(ptr) lis_free(ptr);
    if(index) lis_free(index);
    if(value) lis_free(value);
    if(A) lis_matrix_destroy(A);
    return err;
}

static double probe_value(LIS_INT i)
{
    /* Exactly representable values, varying across columns and rows. */
    return 1.0+((i*17)%233)*0.03125;
}

static int verify(LIS_MATRIX A,LIS_INT side)
{
    LIS_VECTOR ones=NULL,y=NULL;
    LIS_INT err,row,x,z;
    LIS_SCALAR observed;
    int bad=0;

    err=lis_vector_duplicate(A,&ones);
    if(err) return 1;
    err=lis_vector_duplicate(A,&y);
    if(err) { bad=1; goto end; }
    for(row=0;row<side*side;row++) {
        err=lis_vector_set_value(
            LIS_INS_VALUE,row,(LIS_SCALAR)probe_value(row),ones);
        if(err) { bad=1; goto end; }
    }
    err=lis_matvec(A,ones,y);
    if(err) { bad=1; goto end; }
    for(row=0;row<side*side;row++) {
        double expected=4.0*probe_value(row);
        x=row%side; z=row/side;
        if(x>0) expected-=probe_value(row-1);
        if(x+1<side) expected-=probe_value(row+1);
        if(z>0) expected-=probe_value(row-side);
        if(z+1<side) expected-=probe_value(row+side);
        err=lis_vector_get_value(y,row,&observed);
        if(err || fabs((double)observed-expected)>1.0e-12) {
            fprintf(stderr,"SpMV mismatch: row=%ld\n",(long)row);
            bad=1;
            break;
        }
    }
end:
    if(y) lis_vector_destroy(y);
    if(ones) lis_vector_destroy(ones);
    return bad;
}

int main(int argc,char **argv)
{
    LIS_INT side=argc>1?(LIS_INT)atoi(argv[1]):64;
    LIS_INT reps=argc>2?(LIS_INT)atoi(argv[2]):3;
    LIS_INT rep,mode,err;
    LIS_MATRIX A=NULL;
    double best[]={1.e100,1.e100,1.e100},seconds;

    if(side<2 || side>512 || reps<1 || reps>100) {
        fprintf(stderr,"usage: %s [side 2..512] [repeats 1..100]\n",argv[0]);
        return 2;
    }
    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    for(rep=0;rep<reps;rep++) {
        for(mode=0;mode<3;mode++) {
            err=build(side,mode,&A,&seconds);
            if(err) {
                fprintf(stderr,"%s: assembly failed (%ld)\n",names[mode],(long)err);
                goto fail;
            }
            if(verify(A,side)) {
                fprintf(stderr,"%s: verification failed\n",names[mode]);
                goto fail;
            }
            if(seconds<best[mode]) best[mode]=seconds;
            printf("mode=%s rep=%ld n=%ld nnz=%ld assembly_s=%.9f verified=1\n",
                   names[mode],(long)rep,(long)(side*side),
                   (long)(5*side*side-4*side),seconds);
            lis_matrix_destroy(A);
            A=NULL;
        }
    }
    printf("best_seconds,rco=%.9f,rco_preallocated=%.9f,csr=%.9f\n",
           best[0],best[1],best[2]);
    puts("BENCH_RCO_CSR_PASSED");
    lis_finalize();
    return 0;

fail:
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return 1;
}