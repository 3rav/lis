/* Research-only cost-aware rolling preconditioner experiment.
 * Timings are online observations, never read from future steps or other modes. */
#include "lislib.h"
#include "lis_solver_options.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static int valid_matrix(LIS_MATRIX A)
{
#ifdef _OPENMP
    if(omp_in_parallel() || omp_get_max_threads()!=1)return 0;
#endif
    return A && A->label==LIS_LABEL_MATRIX && A->matrix_type==LIS_MATRIX_CSR &&
        A->nprocs==1 && A->n>0 && A->n==A->gn && !A->is_scaled && !A->is_splited &&
        lis_matrix_is_assembled(A) && A->ptr && A->index && A->value;
}
static int valid_vector(LIS_MATRIX A,LIS_VECTOR v)
{
    return v && v->label==LIS_LABEL_VECTOR && v->precision==LIS_PRECISION_DEFAULT &&
        v->nprocs==1 && v->value && v->n==A->n && v->gn==A->gn && v->is==A->is && v->ie==A->ie;
}
static int valid_args(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x)
{
    return valid_matrix(A) && valid_vector(A,b) && valid_vector(A,x) && b!=x && b->value!=x->value;
}
static LIS_INT make_solver(LIS_SOLVER *s)
{
    /* Research-only selection: fresh and retained use exactly this profile. */
    const char *profile=getenv("LIS_ANCHOR_PROFILE"), *precon=NULL;
    char options[512];
    if(!profile)return LIS_ERR_ILL_ARG;
    if(!strcmp(profile,"full-ilutp"))
        precon="-p ilut -ilut_drop 0.005 -ilut_rate 8 -ilu_pivot_tol 1e-9 -ilut_rcm 1 -ilut_row_scale 1 -ilutp 1 -ilutp_threshold 0.5 -ilutp_full 1 -ilutp_full_limit 60000000";
    else if(!strcmp(profile,"incomplete-ilutp"))
        precon="-p ilut -ilut_drop 0.005 -ilut_rate 8 -ilu_pivot_tol 1e-9 -ilut_rcm 1 -ilut_row_scale 1 -ilutp 1 -ilutp_threshold 0.5 -ilutp_full 0";
    else if(!strcmp(profile,"ilut"))precon="-p ilut -ilut_drop 0.005 -ilut_rate 8";
    else if(!strcmp(profile,"ilu0"))precon="-p ilu -ilu_fill 0";
    else if(!strcmp(profile,"ilu1"))precon="-p ilu -ilu_fill 1";
    if(!precon)return LIS_ERR_ILL_ARG;
    if(snprintf(options,sizeof(options),"-i fgmres -scale none -print none -restart 60 -maxiter 500 -tol 1e-8 -initx_zeros true %s",precon)>=(int)sizeof(options))return LIS_ERR_ILL_ARG;
    LIS_INT rc=lis_solver_create(s);
    if(!rc)rc=lis_solver_set_option(options,*s);
    if(!rc)rc=lis_solver_set_num_threads(*s,1);
    return rc;
}
static void unbind(LIS_SOLVER s,LIS_MATRIX P)
{
    /* The caller can release A, b and x after the synchronous operation. */
    if(!s)return;
    s->A=P;s->b=NULL;s->x=NULL;s->xx=NULL;
    (void)lis_solver_set_monitor(s,NULL,NULL);
}
static LIS_INT factor(LIS_MATRIX P,LIS_SOLVER *s,LIS_PRECON *pc,lis_anchor_pilot_report *r)
{
    LIS_INT rc=make_solver(s);
    double start;
    if(rc)return rc;
    start=lis_wtime();
    rc=lis_solver_sequence_prepare_precon(P,*s,pc);
    r->setup_seconds+=lis_wtime()-start;
    if(!rc){
        LIS_INT i;
        long long stored=-1;
        r->factorizations++;
        if((*pc)->L && (*pc)->U && (*pc)->D && (*pc)->L->nnz && (*pc)->U->nnz){
            stored=(*pc)->D->n;
            for(i=0;i<(*pc)->L->n;i++)stored+=(*pc)->L->nnz[i];
            for(i=0;i<(*pc)->U->n;i++)stored+=(*pc)->U->nnz[i];
        }
        printf("FACTOR_NNZ %lld\n",stored);
    }
    return rc;
}
static int numerical(LIS_INT status)
{
    return status==LIS_SUCCESS || status==LIS_MAXITER || status==LIS_BREAKDOWN;
}
static LIS_INT attempt_impl(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR work,LIS_VECTOR residual,
    LIS_SOLVER s,LIS_PRECON pc,LIS_INT budget,int zero,LIS_INT *iterations,
    lis_anchor_pilot_report *r,LIS_FGMRES_STOP stop,void *data)
{
    LIS_INT rc;
    LIS_REAL nr,nb,reported;
    if((rc=lis_solver_set_maxiter(s,budget)))return rc;
    s->options[LIS_OPTIONS_INITGUESS_ZEROS]=zero?LIS_TRUE:LIS_FALSE;
    if(zero && (rc=lis_vector_set_all(0,work)))return rc;
    rc=stop?lis_solve_kernel_fgmres_stop(A,b,work,s,pc,stop,data):lis_solve_kernel(A,b,work,s,pc);
    if(rc)return rc;
    if((rc=lis_solver_get_status(s,&r->raw_status)) ||
       (rc=lis_solver_get_iter(s,iterations)) ||
       (rc=lis_solver_get_residualnorm(s,&reported)))return rc;
    if(!numerical(r->raw_status))return r->raw_status;
    if((rc=lis_matvec(A,work,residual)) || (rc=lis_vector_axpy(-1,b,residual)) ||
       (rc=lis_vector_nrm2(residual,&nr)) || (rc=lis_vector_nrm2(b,&nb)))return rc;
    r->relative_residual=nb>0?nr/nb:nr;
    r->accepted=isfinite((double)reported) && isfinite((double)nr) && isfinite((double)nb) &&
        isfinite(r->relative_residual) && r->relative_residual<=1e-8;
    return LIS_SUCCESS;
}

static LIS_INT attempt(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR work,LIS_VECTOR residual,
    LIS_SOLVER s,LIS_PRECON pc,LIS_INT budget,int zero,LIS_INT *iterations,
    lis_anchor_pilot_report *r)
{
    return attempt_impl(A,b,work,residual,s,pc,budget,zero,iterations,r,NULL,NULL);
}
typedef struct { double r10,r20; LIS_INT last; } pilot_history;
static LIS_INT checkpoint(LIS_INT iter,LIS_REAL res,void *data)
{
    pilot_history *h=data;
    h->last=iter;
    if(iter==10)h->r10=res;
    if(iter!=20)return LIS_FALSE;
    h->r20=res;
    return !(isfinite(h->r10) && h->r10>0 && isfinite(h->r20) &&
             h->r20>=0 && h->r20/h->r10<=0.25);
}
static void monitor(LIS_SOLVER s,LIS_INT iter,LIS_REAL res,LIS_REAL tol,LIS_INT kind,void *data)
{
    pilot_history *h=data;(void)s;(void)tol;
    if(kind!=LIS_SOLVER_RESIDUAL_PROJECTED)return;
    h->last=iter;
    if(iter==10)h->r10=res;
    if(iter==20)h->r20=res;
}

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

typedef struct { LIS_SOLVER s; LIS_PRECON pc; LIS_MATRIX P; } anchor;
typedef struct { pilot_history history; double start,budget,prev_res; LIS_INT prev_iter; int reason; } cost_history;
static void reset(lis_anchor_pilot_report *r)
{
    memset(r,0,sizeof(*r));r->raw_status=LIS_ERR_ILL_ARG;r->relative_residual=r->ratio=NAN;
}
static void release(anchor *a)
{
    if(a->s)lis_solver_destroy(a->s);
    memset(a,0,sizeof(*a));
}
static LIS_INT build(anchor *a,LIS_MATRIX A,lis_anchor_pilot_report *r)
{
    a->P=A;return factor(A,&a->s,&a->pc,r);
}
static LIS_INT fresh(anchor *a,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR residual,lis_anchor_pilot_report *r)
{
    pilot_history h={NAN,NAN,0};
    LIS_INT rc=lis_solver_set_monitor(a->s,monitor,&h);
    if(!rc)rc=attempt(A,b,x,residual,a->s,a->pc,500,1,&r->fresh_iterations,r);
    r->fresh_iterations=h.last;
    unbind(a->s,a->P);
    return rc;
}
/* Conservative projected time to the same residual target, checked every
 * five iterations after ten. At most 80% of last observed fresh total.
 * Actual elapsed pilot time includes overhead before the first callback. */
static LIS_INT cost_stop(LIS_INT iter,LIS_REAL res,void *data)
{
    cost_history *h=data;
    double elapsed=lis_wtime()-h->start;
    h->history.last=iter;
    if(iter==10)h->history.r10=res;
    if(iter==20)h->history.r20=res;
    if(!isfinite((double)res)||res<0){h->reason=1;return LIS_TRUE;}
    if(res<=1e-8)return LIS_FALSE;
    if(elapsed>=h->budget){h->reason=2;return LIS_TRUE;}
    if(iter%5==0){
        if(h->prev_iter>0){
            double slope=log(res/h->prev_res)/(iter-h->prev_iter);
            double remain=slope<0?log(1e-8/res)/slope:INFINITY;
            double predicted=elapsed+1.25*remain*elapsed/iter;
            if(!isfinite(predicted)||predicted>h->budget||remain>100-iter){
                h->reason=3;return LIS_TRUE;
            }
        }
        h->prev_iter=iter;h->prev_res=res;
    }
    return LIS_FALSE;
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL,residual=NULL;
    anchor kept={0},candidate={0};
    int i,steps,accepted=0,mode,hard=0,streak=0,cooldown=0;
    LIS_INT rc=0;
    double last_fresh=0,last_iter_cost=0,total=0,setup=0;
    int creates=0;
    const char *audit=getenv("LIS_ANCHOR_AUDIT_DIR");
    if(argc<4||argc%2)return 2;
    if(!strcmp(argv[1],"fresh"))mode=0;
    else if(!strcmp(argv[1],"fixed"))mode=1;
    else if(!strcmp(argv[1],"rolling"))mode=2;
    else if(!strcmp(argv[1],"cost-aware"))mode=3;
    else return 2;
    steps=(argc-2)/2;
    if(lis_initialize(&argc,&argv))return 2;
    for(i=2;i<argc;i+=2){
        lis_anchor_pilot_report r;
        pilot_history h={NAN,NAN,0};
        cost_history ch={{NAN,NAN,0},0,0,0,0,0};
        int reason=0,skip=0,promote=0,attempted=0;
        double start,elapsed,fs=0,fc=0;
        LIS_INT reported=0;
        reset(&r);
        if(load_matrix(argv[i],&A)||lis_vector_duplicate(A,&b)||lis_vector_duplicate(A,&x)||
           lis_vector_duplicate(A,&residual)||load_rhs(argv[i+1],b)||!valid_args(A,b,x)){
            hard=1;goto done;
        }
        start=lis_wtime();
        if(mode==1&&!kept.s){
            rc=build(&kept,A,&r);
            if(rc){hard=1;goto done;}
        }
        if(mode==3&&kept.s){
            if(cooldown>0){cooldown--;skip=1;reason=4;}
            else if(last_fresh>0&&5*last_iter_cost>.8*last_fresh){skip=1;reason=5;}
        }
        if(mode&&kept.s&&!skip){
            attempted=1;
            if(mode==3){
                ch.start=lis_wtime();ch.budget=.8*last_fresh;
                rc=attempt_impl(A,b,x,residual,kept.s,kept.pc,100,1,&reported,&r,cost_stop,&ch);
                h=ch.history;reason=ch.reason;
            }else{
                rc=attempt_impl(A,b,x,residual,kept.s,kept.pc,100,1,&reported,&r,checkpoint,&h);
                if(!r.accepted)reason=6;
            }
            unbind(kept.s,kept.P);
            if(rc){hard=1;goto done;}
            r.pilot_iterations=h.last<20?h.last:20;
            r.continuation_iterations=h.last>20?h.last-20:0;
            r.continued=h.last>20;
            if(isfinite(h.r10)&&h.r10>0&&isfinite(h.r20)&&h.r20>=0)r.ratio=h.r20/h.r10;
            if(mode==3){
                if(r.accepted)streak=0;
                else if(++streak>=2){cooldown=2;streak=0;}
            }
        }
        if(!r.accepted){
            r.fallback=attempted;
            fs=lis_wtime();
            rc=build(&candidate,A,&r);
            if(rc){hard=1;goto done;}
            fc=lis_wtime();
            rc=fresh(&candidate,A,b,x,residual,&r);
            if(rc){hard=1;goto done;}
            if(mode>=2){
                last_fresh=lis_wtime()-fs;
                if(r.fresh_iterations>0)last_iter_cost=(lis_wtime()-fc)/r.fresh_iterations;
                if(r.accepted)promote=1;
            }
            if(promote){
                LIS_MATRIX old=kept.P;
                release(&kept);
                kept=candidate;memset(&candidate,0,sizeof(candidate));
                if(old&&old!=A)lis_matrix_destroy(old);
            }else release(&candidate);
        }
        elapsed=lis_wtime()-start;
        accepted+=r.accepted;creates+=r.factorizations;total+=elapsed;setup+=r.setup_seconds;
        printf("COST step=%d accepted=%d api=%lld raw_status=%lld pilot_iter=%lld continuation_iter=%lld fresh_iter=%lld continued=%d fallback=%d factors=%d ratio=%.17e residual=%.17e setup_s=%.9f total_s=%.9f skip=%d reason=%d promote=%d attempted=%d\n",
            i/2,r.accepted,(long long)rc,(long long)r.raw_status,(long long)r.pilot_iterations,(long long)r.continuation_iterations,(long long)r.fresh_iterations,r.continued,r.fallback,r.factorizations,r.ratio,r.relative_residual,r.setup_seconds,elapsed,skip,reason,promote,attempted);
        fflush(stdout);
        if(audit&&r.accepted&&save_solution(audit,i/2,x)){hard=1;goto done;}
        lis_vector_destroy(x);x=NULL;lis_vector_destroy(b);b=NULL;lis_vector_destroy(residual);residual=NULL;
        if(A!=kept.P)lis_matrix_destroy(A);
        A=NULL;
    }
    printf("COST_SUMMARY accepted=%d/%d factors=%d setup_s=%.9f total_s=%.9f\n",accepted,steps,creates,setup,total);
done:
    release(&candidate);
    if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(residual)lis_vector_destroy(residual);
    {LIS_MATRIX P=kept.P;release(&kept);if(A&&A!=P)lis_matrix_destroy(A);if(P)lis_matrix_destroy(P);}
    lis_finalize();
    return hard?2:accepted!=steps;
}
