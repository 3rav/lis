/* Research-only fresh preconditioner and direct solver comparison.
 * All modes load identical CSR/RHS files and use independent physical audits.
 * UMFPACK symbolic cache is online, keyed by exact CSC pattern, no lookahead. */
#include "lislib.h"
#include <suitesparse/umfpack.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef char check_native_int[(sizeof(LIS_INT)==sizeof(int))?1:-1];
static LIS_INT load_matrix(char *name,LIS_MATRIX *out)
{
    LIS_INT rc=lis_matrix_create(LIS_COMM_WORLD,out);
    if(!rc)rc=lis_matrix_set_type(*out,LIS_MATRIX_CSR);
    if(!rc)rc=lis_input_matrix(*out,name);
    if(!rc)rc=lis_matrix_sort_csr(*out);
    return rc;
}
/* Standard Matrix Market real general array, one column. */
static LIS_INT load_rhs(const char *name,LIS_VECTOR b)
{
    FILE *f=fopen(name,"r");char line[512];long long n,m;LIS_INT i;
    if(!f)return LIS_ERR_FILE_IO;
    if(!fgets(line,sizeof(line),f) || strncmp(line,"%%MatrixMarket matrix array real general",40)!=0)goto bad;
    do{if(!fgets(line,sizeof(line),f))goto bad;}while(line[0]=='%');
    if(sscanf(line,"%lld%lld",&n,&m)!=2 || n!=b->n || m!=1)goto bad;
    for(i=0;i<b->n;i++)if(fscanf(f,"%lf",&b->value[i])!=1 || !isfinite(b->value[i]))goto bad;
    fclose(f);return LIS_SUCCESS;
bad:
    fclose(f);return LIS_ERR_FILE_IO;
}
/* Opt-in audit only: exact round-trip decimal doubles, no timing claims. */
static LIS_INT save_solution(const char *dir,int step,LIS_VECTOR x)
{
    char path[4096]; FILE *f; LIS_INT i; int len,failed=0;
    len=snprintf(path,sizeof(path),"%s/x-%03d.mtx",dir,step);
    if(len<0 || (size_t)len>=sizeof(path))return LIS_ERR_FILE_IO;
    f=fopen(path,"wx");
    if(!f)return LIS_ERR_FILE_IO;
    if(fprintf(f,"%%%%MatrixMarket matrix array real general\n%lld 1\n",(long long)x->n)<0)failed=1;
    for(i=0;i<x->n && !failed;i++){
        if(!isfinite(x->value[i]) || fprintf(f,"%.17g\n",x->value[i])<0)failed=1;
    }
    if(fclose(f)!=0)failed=1;
    return failed?LIS_ERR_FILE_IO:LIS_SUCCESS;
}

typedef struct { int n,nnz,*ap,*ai; void *sym; } entry;
typedef struct { entry *items; int used,capacity; } cache;
typedef struct {
    int accepted, factors, symbolic, status, api, iterations;
    long long fill;
    double setup,total,residual;
} report;
static void monitor(LIS_SOLVER s,LIS_INT iter,LIS_REAL res,LIS_REAL tol,LIS_INT kind,void *data)
{
    report *r=data;(void)s;(void)res;(void)tol;
    if(kind==LIS_SOLVER_RESIDUAL_PROJECTED)r->iterations=(int)iter;
}
static long long factor_nnz(LIS_PRECON pc)
{
    LIS_INT i;
    long long count=0;
    if(!pc || !pc->L || !pc->U || !pc->D || !pc->L->nnz || !pc->U->nnz)return -1;
    count=pc->D->n;
    for(i=0;i<pc->L->n;i++)count+=pc->L->nnz[i];
    for(i=0;i<pc->U->n;i++)count+=pc->U->nnz[i];
    return count;
}
static int physical(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR scratch,report *r)
{
    LIS_REAL nr,nb;
    LIS_INT rc;
    if((rc=lis_matvec(A,x,scratch)) || (rc=lis_vector_axpy(-1,b,scratch)) ||
       (rc=lis_vector_nrm2(scratch,&nr)) || (rc=lis_vector_nrm2(b,&nb)))return (int)rc;
    r->residual=nb>0?nr/nb:nr;
    r->accepted=isfinite(nr)&&isfinite(nb)&&isfinite(r->residual)&&r->residual<=1e-8;
    return 0;
}
static int native(const char *mode,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,report *r)
{
    char opts[512];
    const char *profile=NULL;
    LIS_SOLVER s=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT rc=0,status;
    LIS_REAL reported;
    double start;
    if(!strcmp(mode,"fresh-ilut"))profile="-p ilut -ilut_drop 0.005 -ilut_rate 8";
    else if(!strcmp(mode,"fresh-ilutp-incomplete"))profile="-p ilut -ilut_drop 0.005 -ilut_rate 8 -ilu_pivot_tol 1e-9 -ilut_rcm 1 -ilut_row_scale 1 -ilutp 1 -ilutp_threshold 0.5 -ilutp_full 0";
    else if(!strcmp(mode,"fresh-ilu0"))profile="-p ilu -ilu_fill 0";
    else if(!strcmp(mode,"fresh-ilu1"))profile="-p ilu -ilu_fill 1";
    else if(!strcmp(mode,"fresh-mfamg"))profile="-p mfamg";
    if(!profile)return LIS_ERR_ILL_ARG;
    snprintf(opts,sizeof(opts),"-i fgmres -scale none -print none -restart 60 -maxiter 500 -tol 1e-8 -initx_zeros true %s",profile);
    if((rc=lis_solver_create(&s)) || (rc=lis_solver_set_option(opts,s)) ||
       (rc=lis_solver_set_num_threads(s,1)) || (rc=lis_solver_set_monitor(s,monitor,r)))goto done;
    start=lis_wtime();
    rc=lis_solver_sequence_prepare_precon(A,s,&pc);
    r->setup=lis_wtime()-start;
    if(rc)goto done;
    r->factors=1;r->fill=factor_nnz(pc);
    rc=lis_solve_kernel(A,b,x,s,pc);
    if(rc)goto done;
    if((rc=lis_solver_get_status(s,&status)) ||
       (rc=lis_solver_get_residualnorm(s,&reported)))goto done;
    r->status=(int)status;
    if(!isfinite(reported))rc=LIS_BREAKDOWN;
    if(status!=LIS_SUCCESS && status!=LIS_MAXITER && status!=LIS_BREAKDOWN)rc=status;
done:
    if(s)lis_solver_destroy(s);
    return (int)rc;
}
static void free_cache(cache *c)
{
    int i;
    for(i=0;i<c->used;i++){
        umfpack_di_free_symbolic(&c->items[i].sym);
        free(c->items[i].ap);free(c->items[i].ai);
    }
    free(c->items);
}
static int direct(int reuse,cache *c,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,report *r)
{
    int n=(int)A->n,nnz=(int)A->nnz,i,rc=UMFPACK_OK,owns=1;
    int *ap=NULL,*ai=NULL;
    double *ax=NULL,control[UMFPACK_CONTROL],info[UMFPACK_INFO];
    void *sym=NULL,*num=NULL;
    double start=lis_wtime();
    ap=malloc(((size_t)n+1)*sizeof(*ap));ai=malloc((size_t)nnz*sizeof(*ai));ax=malloc((size_t)nnz*sizeof(*ax));
    if(!ap||!ai||!ax){rc=UMFPACK_ERROR_out_of_memory;goto done;}
    /* CSR of A is CSC of A^T; transpose to canonical CSC of A. */
    rc=umfpack_di_transpose(n,n,A->ptr,A->index,A->value,NULL,NULL,ap,ai,ax);
    if(rc!=UMFPACK_OK)goto done;
    umfpack_di_defaults(control);
    if(reuse)for(i=0;i<c->used;i++){
        entry *e=&c->items[i];
        if(e->n==n&&e->nnz==nnz&&!memcmp(e->ap,ap,((size_t)n+1)*sizeof(*ap))&&!memcmp(e->ai,ai,(size_t)nnz*sizeof(*ai))){
            sym=e->sym;owns=0;break;
        }
    }
    if(!sym){
        rc=umfpack_di_symbolic(n,n,ap,ai,ax,&sym,control,info);
        if(rc!=UMFPACK_OK)goto done;
        r->symbolic=1;
        if(reuse){
            entry *e;
            if(c->used>=c->capacity){rc=UMFPACK_ERROR_out_of_memory;goto done;}
            e=&c->items[c->used];
            e->ap=malloc(((size_t)n+1)*sizeof(*ap));e->ai=malloc((size_t)nnz*sizeof(*ai));
            if(!e->ap||!e->ai){free(e->ap);free(e->ai);rc=UMFPACK_ERROR_out_of_memory;goto done;}
            memcpy(e->ap,ap,((size_t)n+1)*sizeof(*ap));memcpy(e->ai,ai,(size_t)nnz*sizeof(*ai));
            e->n=n;e->nnz=nnz;e->sym=sym;owns=0;c->used++;
        }
    }
    rc=umfpack_di_numeric(ap,ai,ax,sym,&num,control,info);
    r->setup=lis_wtime()-start;
    if(rc!=UMFPACK_OK)goto done;
    r->factors=1;r->fill=(long long)(info[UMFPACK_LNZ]+info[UMFPACK_UNZ])-n;
    rc=umfpack_di_solve(UMFPACK_A,ap,ai,ax,x->value,b->value,num,control,info);
done:
    if(r->setup==0)r->setup=lis_wtime()-start;
    if(num)umfpack_di_free_numeric(&num);
    if(sym&&owns)umfpack_di_free_symbolic(&sym);
    free(ap);free(ai);free(ax);
    r->status=rc;
    return rc;
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL,rvec=NULL;
    int i,steps,accepted=0,rc=0,hard=0,is_direct,reuse;
    cache c={0};
    const char *audit=getenv("LIS_ANCHOR_AUDIT_DIR");
    if(argc<4||argc%2)return 2;
    is_direct=!strcmp(argv[1],"umfpack-fresh")||!strcmp(argv[1],"umfpack-symbolic-cache");
    reuse=!strcmp(argv[1],"umfpack-symbolic-cache");
    steps=(argc-2)/2;
    c.capacity=steps;c.items=calloc((size_t)steps,sizeof(*c.items));
    if(!c.items)return 2;
    if(lis_initialize(&argc,&argv)){free_cache(&c);return 2;}
    for(i=2;i<argc;i+=2){
        report r={0};double start;
        r.fill=-1;r.residual=NAN;
        if(load_matrix(argv[i],&A)||lis_vector_duplicate(A,&b)||lis_vector_duplicate(A,&x)||
           lis_vector_duplicate(A,&rvec)||load_rhs(argv[i+1],b)){hard=1;goto done;}
        start=lis_wtime();
        r.api=is_direct?direct(reuse,&c,A,b,x,&r):native(argv[1],A,b,x,&r);
        if(!r.api)r.api=physical(A,b,x,rvec,&r);
        if(r.api)r.accepted=0;
        r.total=lis_wtime()-start;
        printf("BASE step=%d accepted=%d api=%d status=%d iterations=%d factors=%d symbolic=%d fill=%lld residual=%.17e setup_s=%.9f total_s=%.9f\n",
            i/2,r.accepted,r.api,r.status,r.iterations,r.factors,r.symbolic,r.fill,r.residual,r.setup,r.total);
        fflush(stdout);
        accepted+=r.accepted;
        if(audit&&r.accepted&&save_solution(audit,i/2,x)){hard=1;goto done;}
        lis_vector_destroy(rvec);rvec=NULL;lis_vector_destroy(x);x=NULL;
        lis_vector_destroy(b);b=NULL;lis_matrix_destroy(A);A=NULL;
    }
    printf("BASE_SUMMARY accepted=%d/%d\n",accepted,steps);
done:
    if(rvec)lis_vector_destroy(rvec);if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(A)lis_matrix_destroy(A);
    free_cache(&c);lis_finalize();
    rc=hard?2:(accepted==steps?0:1);return rc;
}
