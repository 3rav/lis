/* Microbenchmark of row-backed USER -> cached CSR on serial and OpenMP paths.
 * Validation-only. Each measurement uses a fresh cache for build and a reused
 * cache for values-only refresh. Neither affects a feature PR diff.
 */
#include "lis.h"
#include "lis_matrix.h"
#include "lis_matrix_materialize.h"
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

typedef struct {
    LIS_INT n, radius;
    LIS_REAL scale;
} BENCH_CTX;

static LIS_INT bench_row(void *data,LIS_INT row,LIS_INT *nnz,
                         LIS_INT *cols,LIS_SCALAR *vals)
{
    BENCH_CTX *c=(BENCH_CTX *)data;
    LIS_INT i,first,last,need,k=0;
    if(!c || !nnz || row<0 || row>=c->n) return LIS_ERR_ILL_ARG;
    first=row>c->radius ? row-c->radius : 0;
    last=row+c->radius<c->n ? row+c->radius : c->n-1;
    need=last-first+1;
    if(!cols && !vals){*nnz=need;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<need){*nnz=need;return LIS_ERR_ILL_ARG;}
    for(i=first;i<=last;i++)
    {
        cols[k]=i;
        vals[k]=(LIS_SCALAR)((i==row ? 8.0 : -0.125)*c->scale);
        k++;
    }
    *nnz=need;
    return LIS_SUCCESS;
}

static LIS_INT bench_matvec(void *data,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    BENCH_CTX *c=(BENCH_CTX *)data;
    LIS_INT row,k,first,last;
    if(!c || !x || !y)return LIS_ERR_ILL_ARG;
    for(row=0;row<c->n;row++)
    {
        LIS_SCALAR sum=(LIS_SCALAR)0.0;
        first=row>c->radius ? row-c->radius : 0;
        last=row+c->radius<c->n ? row+c->radius : c->n-1;
        for(k=first;k<=last;k++)
            sum+=(LIS_SCALAR)((k==row ? 8.0 : -0.125)*c->scale)*x[k];
        y[row]=sum;
    }
    return LIS_SUCCESS;
}

static double clock_seconds(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
}

static LIS_INT benchmark(LIS_INT n,LIS_INT radius,int threads,int safe,
                         int repeats,double *build_ms,double *values_ms)
{
    BENCH_CTX c;
    LIS_MATRIX A=NULL,csr=NULL;
    LIS_MATRIX_CSR_CACHE cache=NULL;
    LIS_INT err,action,rep;
    double tb=0.0,tv=0.0,t0;
    int refreshed=0;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#endif
    c.n=n;c.radius=radius;c.scale=(LIS_REAL)1.0;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err)return err;
    err=lis_matrix_set_size(A,0,n);
    if(err)goto cleanup;
    err=lis_matrix_set_user(A,&c,bench_matvec,bench_matvec);
    if(err)goto cleanup;
    err=lis_matrix_set_user_row(A,bench_row);
    if(err)goto cleanup;
    if(safe)
    {
        err=lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);
        if(err)goto cleanup;
    }
    for(rep=0;rep<repeats;rep++)
    {
        t0=clock_seconds();
        err=lis_matrix_materialize_csr_cached(A,&cache,&csr,&action);
        tb+=clock_seconds()-t0;
        if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD)goto cleanup;
        c.scale+=(LIS_REAL)0.125;
        err=lis_matrix_mark_values_changed(A);
        if(err)goto cleanup;
        t0=clock_seconds();
        err=lis_matrix_materialize_csr_cached(A,&cache,&csr,&action);
        tv+=clock_seconds()-t0;
        if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES)goto cleanup;
        refreshed++;
        lis_matrix_csr_cache_destroy(cache);
        cache=NULL;
        csr=NULL;
    }
    *build_ms=1000.0*tb/repeats;
    *values_ms=1000.0*tv/refreshed;
cleanup:
    if(cache)lis_matrix_csr_cache_destroy(cache);
    if(A)lis_matrix_destroy(A);
    return err;
}

int main(int argc,char **argv)
{
    const LIS_INT sizes[]={256,2048,16384};
    const LIS_INT radii[]={1,8};
    const int counts[]={1,2,4};
    int s,r,t,safe;
    LIS_INT err;
    err=lis_initialize(&argc,&argv);
    if(err)return 1;
    puts("size,radius,threads,callback_safe,build_ms,values_ms");
    for(s=0;s<3;s++)
      for(r=0;r<2;r++)
        for(t=0;t<3;t++)
          for(safe=0;safe<=1;safe++)
          {
            double build=0,values=0;
            err=benchmark(sizes[s],radii[r],counts[t],safe,3,
                          &build,&values);
            if(err)
            {
                fprintf(stderr,"bench failed n=%d radius=%d threads=%d safe=%d code=%d\n",
                        (int)sizes[s],(int)radii[r],counts[t],safe,(int)err);
                lis_finalize();
                return 1;
            }
            printf("%d,%d,%d,%d,%.6f,%.6f\n",
                   (int)sizes[s],(int)radii[r],counts[t],safe,build,values);
          }
    lis_finalize();
    return 0;
}
