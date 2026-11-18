/* Optional batched USER row integration and cache lifecycle regression. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    LIS_INT n,scalar_queries,batch_queries,fail_mismatch;
} CONTEXT;

static LIS_INT count_row(CONTEXT *ctx,LIS_INT row)
{
    return (row>0)+(LIS_INT)1+(row+1<ctx->n);
}

static LIS_INT get_row(void *opaque,LIS_INT row,LIS_INT *nnz,
                       LIS_INT *cols,LIS_SCALAR *vals)
{
    CONTEXT *ctx=(CONTEXT *)opaque;
    LIS_INT k=0,need;
    if(row<0 || row>=ctx->n || nnz==NULL) return LIS_ERR_ILL_ARG;
    need=count_row(ctx,row);
#ifdef _OPENMP
#pragma omp atomic update
#endif
    ctx->scalar_queries++;
    if(cols==NULL && vals==NULL){*nnz=need;return LIS_SUCCESS;}
    if(cols==NULL || vals==NULL || *nnz<need)return LIS_ERR_ILL_ARG;
    if(row>0){cols[k]=row-1;vals[k++]=-1.0;}
    cols[k]=row;vals[k++]=4.0;
    if(row+1<ctx->n){cols[k]=row+1;vals[k++]=-1.0;}
    *nnz=need;
    return LIS_SUCCESS;
}

static LIS_INT get_batch(void *opaque,LIS_INT first,LIS_INT last,
                         LIS_INT *offsets,LIS_INT capacity,
                         LIS_INT *cols,LIS_SCALAR *vals)
{
    CONTEXT *ctx=(CONTEXT *)opaque;
    LIS_INT row,k=0,i;
    if(!offsets || first<0 || last>ctx->n || first>=last)
        return LIS_ERR_ILL_ARG;
#ifdef _OPENMP
#pragma omp atomic update
#endif
    ctx->batch_queries++;
    offsets[0]=0;
    for(row=first;row<last;row++)
        offsets[row-first+1]=offsets[row-first]+count_row(ctx,row);
    if(!cols && !vals)return LIS_SUCCESS;
    if(!cols || !vals || capacity<offsets[last-first])
        return LIS_ERR_ILL_ARG;
    for(row=first;row<last;row++)
    {
        if(row>0){cols[k]=row-1;vals[k++]=-1.0;}
        cols[k]=row;vals[k++]=4.0;
        if(row+1<ctx->n){cols[k]=row+1;vals[k++]=-1.0;}
    }
    if(ctx->fail_mismatch)
    {
        i=last-first;
        offsets[i]++; /* Deliberate contract violation on fill. */
    }
    return LIS_SUCCESS;
}

static int run_mode(int use_batch)
{
    enum { N=256 };
    CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_PRECON_ROW_SOURCE src;
    LIS_INT err=0,i,j,k,nnz;
    int failed=1,initialized=0;
    memset(&ctx,0,sizeof(ctx));
    memset(&src,0,sizeof(src));
    ctx.n=N;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,get_row);if(err)goto done;
    if(use_batch)
    {
        err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    }
    err=lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&y);if(err)goto done;
    for(i=0;i<N;i++) x->value[i]=1.0+0.01*i;

    err=lis_matvec(A,x,y);if(err)goto done;
    for(i=0;i<N;i++)
    {
        double expected=4.0*x->value[i];
        if(i>0)expected-=x->value[i-1];
        if(i+1<N)expected-=x->value[i+1];
        if(fabs((double)y->value[i]-expected)>1e-11)goto done;
    }
    err=lis_matvech(A,x,y);if(err)goto done;
    for(i=0;i<N;i++)
    {
        double expected=4.0*x->value[i];
        if(i>0)expected-=x->value[i-1];
        if(i+1<N)expected-=x->value[i+1];
        if(fabs((double)y->value[i]-expected)>1e-11)goto done;
    }
    err=lis_precon_row_source_init(A,&src);if(err)goto done;
    initialized=1;
    /* Non-monotone order must trigger safe cache refresh. */
    for(k=0;k<3;k++)
    {
        LIS_INT row=k==0?97:(k==1?3:99);
        err=lis_precon_row_source_get(&src,row);if(err)goto done;
        nnz=count_row(&ctx,row);
        if(src.nnz!=nnz)goto done;
        for(j=0;j<src.nnz;j++)
        {
            if(src.index[j]<0 || src.index[j]>=N)goto done;
        }
    }

    if(use_batch && ctx.scalar_queries!=0)
    {
        fprintf(stderr,"unexpected scalar callback count %d\n",
                (int)ctx.scalar_queries);
        goto done;
    }
    if(use_batch && ctx.batch_queries>=2*N)
    {
        fprintf(stderr,"batch callback count too high %d\n",
                (int)ctx.batch_queries);
        goto done;
    }
    if(!use_batch && ctx.scalar_queries<2*N)goto done;

    if(use_batch)
    {
        /* Explicit mismatch must be rejected, not produce stale row data. */
        ctx.fail_mismatch=1;
        err=lis_precon_row_source_get(&src,210);
        if(err!=LIS_ERR_ILL_ARG)goto done;
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"USER_BATCH FAIL mode=%d err=%d scalar=%d batch=%d\n",
                use_batch,(int)err,(int)ctx.scalar_queries,(int)ctx.batch_queries);
    else
        printf("USER_BATCH PASS mode=%d scalar_calls=%d batch_calls=%d\n",
               use_batch,(int)ctx.scalar_queries,(int)ctx.batch_queries);
    if(initialized)lis_precon_row_source_destroy(&src);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int failed;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 2;
    failed=run_mode(0)|run_mode(1);
    lis_finalize();
    return failed;
}
