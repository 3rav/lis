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

/* Only the research benchmark owns these objects. A is sorted CSR;
 * the anchor remains alive and immutable until its factor is released. */
typedef struct { LIS_SOLVER s; LIS_PRECON pc; LIS_MATRIX P; } anchor;
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
/* Max relative row-L1 change on the union of the sorted CSR patterns.
 * Row scaling prevents large rows hiding a change in a small row.
 * This is a heuristic, not a bound on preconditioned convergence. */
static double distance(LIS_MATRIX A,LIS_MATRIX P)
{
    LIS_INT i,j,k,je,ke;
    double worst=0;
    if(!A||!P||A->n!=P->n)return INFINITY;
    for(i=0;i<A->n;i++){
        double delta=0,norm=0;
        j=A->ptr[i];je=A->ptr[i+1];k=P->ptr[i];ke=P->ptr[i+1];
        while(j<je||k<ke){
            if(k==ke||(j<je&&A->index[j]<P->index[k])){delta+=fabs(A->value[j++]);}
            else if(j==je||P->index[k]<A->index[j]){
                double p=P->value[k++];delta+=fabs(p);norm+=fabs(p);
            }else{
                double p=P->value[k++];delta+=fabs(A->value[j++]-p);norm+=fabs(p);
            }
        }
        if(!isfinite(delta)||!isfinite(norm))return INFINITY;
        if(norm==0){if(delta>0)return INFINITY;}
        else if(delta/norm>worst)worst=delta/norm;
    }
    return worst;
}
/* A warm guess changes NRM2_R's normalization. Convert the solver tolerance
 * from the physical RHS-relative target to the initial-residual target.
 * Reinitialize the guess before every attempt, including fallback. */
static LIS_INT solve_attempt(anchor *a,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,
    LIS_VECTOR residual,LIS_VECTOR prev,int retained,lis_anchor_pilot_report *r,pilot_history *h)
{
    LIS_INT rc,it=0;
    LIS_REAL nr,nb;
    double tol=1e-8;
    r->accepted=0;r->relative_residual=NAN;r->raw_status=LIS_ERR_ILL_ARG;
    if(prev){
        if((rc=lis_vector_copy(prev,x))||(rc=lis_matvec(A,x,residual))||
           (rc=lis_vector_axpy(-1,b,residual))||(rc=lis_vector_nrm2(residual,&nr))||
           (rc=lis_vector_nrm2(b,&nb)))return rc;
        r->relative_residual=nb>0?nr/nb:nr;
        if(isfinite(nr)&&isfinite(nb)&&isfinite(r->relative_residual)){
            if(r->relative_residual<=1e-8){
                r->raw_status=LIS_SUCCESS;r->accepted=1;return LIS_SUCCESS;
            }
            tol=1e-8/r->relative_residual;
        }else return LIS_ERR_ILL_ARG;
    }
    if((rc=lis_solver_set_tolerance(a->s,tol)))return rc;
    if(!retained&&(rc=lis_solver_set_monitor(a->s,monitor,h)))return rc;
    rc=attempt_impl(A,b,x,residual,a->s,a->pc,retained?100:500,prev==NULL,&it,r,
                    retained?checkpoint:NULL,retained?h:NULL);
    unbind(a->s,a->P);
    return rc;
}
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL,residual=NULL,prev=NULL;
    anchor kept={0},candidate={0};
    int i,steps,accepted=0,mode,warm,similar,hard=0,creates=0;
    LIS_INT rc=0;
    double total=0,setup=0;
    const char *audit=getenv("LIS_ANCHOR_AUDIT_DIR");
    if(argc<4||argc%2)return 2;
    if(!strcmp(argv[1],"fresh"))mode=0;
    else if(!strcmp(argv[1],"fresh-warm"))mode=1;
    else if(!strcmp(argv[1],"rolling"))mode=2;
    else if(!strcmp(argv[1],"rolling-warm"))mode=3;
    else if(!strcmp(argv[1],"similarity"))mode=4;
    else if(!strcmp(argv[1],"similarity-warm"))mode=5;
    else return 2;
    warm=mode%2;similar=mode>=4;steps=(argc-2)/2;
    if(lis_initialize(&argc,&argv))return 2;
    for(i=2;i<argc;i+=2){
        lis_anchor_pilot_report r;
        pilot_history h={NAN,NAN,0};
        int reason=0,skip=0,promote=0,attempted=0;
        double start,elapsed,change=NAN;
        reset(&r);
        if(load_matrix(argv[i],&A)||lis_vector_duplicate(A,&b)||lis_vector_duplicate(A,&x)||
           lis_vector_duplicate(A,&residual)||load_rhs(argv[i+1],b)||!valid_args(A,b,x)){
            hard=1;goto done;
        }
        start=lis_wtime();
        if(mode>=2&&kept.s){
            if(similar){
                change=distance(A,kept.P);
                if(!isfinite(change)||change>.05){skip=1;reason=7;}
            }
            if(!skip){
                attempted=1;
                rc=solve_attempt(&kept,A,b,x,residual,prev,1,&r,&h);
                if(rc){hard=1;goto done;}
                r.pilot_iterations=h.last<20?h.last:20;
                r.continuation_iterations=h.last>20?h.last-20:0;
                r.continued=h.last>20;
                if(isfinite(h.r10)&&h.r10>0&&isfinite(h.r20)&&h.r20>=0)r.ratio=h.r20/h.r10;
                if(!r.accepted)reason=6;
            }
        }
        if(!r.accepted){
            pilot_history f={NAN,NAN,0};
            r.fallback=attempted;
            /* Retire the old factor before fresh setup for all rolling modes;
             * no simultaneous live factors in this research comparison. */
            if(kept.s){
                LIS_MATRIX old=kept.P;
                release(&kept);if(old!=A)lis_matrix_destroy(old);
            }
            rc=build(&candidate,A,&r);
            if(rc){hard=1;goto done;}
            rc=solve_attempt(&candidate,A,b,x,residual,prev,0,&r,&f);
            if(rc){hard=1;goto done;}
            r.fresh_iterations=f.last;
            if(mode>=2&&r.accepted){
                kept=candidate;memset(&candidate,0,sizeof(candidate));promote=1;
            }else release(&candidate);
        }
        if(warm&&r.accepted){
            if(!prev&&(rc=lis_vector_duplicate(A,&prev))){hard=1;goto done;}
            if((rc=lis_vector_copy(x,prev))){hard=1;goto done;}
        }
        elapsed=lis_wtime()-start;
        accepted+=r.accepted;creates+=r.factorizations;total+=elapsed;setup+=r.setup_seconds;
        printf("COST step=%d accepted=%d api=%lld raw_status=%lld pilot_iter=%lld continuation_iter=%lld fresh_iter=%lld continued=%d fallback=%d factors=%d ratio=%.17e residual=%.17e setup_s=%.9f total_s=%.9f skip=%d reason=%d promote=%d attempted=%d change=%.17e\n",
            i/2,r.accepted,(long long)rc,(long long)r.raw_status,(long long)r.pilot_iterations,(long long)r.continuation_iterations,(long long)r.fresh_iterations,r.continued,r.fallback,r.factorizations,r.ratio,r.relative_residual,r.setup_seconds,elapsed,skip,reason,promote,attempted,change);
        fflush(stdout);
        if(audit&&r.accepted&&save_solution(audit,i/2,x)){hard=1;goto done;}
        lis_vector_destroy(x);x=NULL;lis_vector_destroy(b);b=NULL;lis_vector_destroy(residual);residual=NULL;
        if(A!=kept.P)lis_matrix_destroy(A);
        A=NULL;
    }
    printf("COST_SUMMARY accepted=%d/%d factors=%d setup_s=%.9f total_s=%.9f\n",accepted,steps,creates,setup,total);
done:
    release(&candidate);
    if(prev)lis_vector_destroy(prev);
    if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(residual)lis_vector_destroy(residual);
    {LIS_MATRIX P=kept.P;release(&kept);if(A&&A!=P)lis_matrix_destroy(A);if(P)lis_matrix_destroy(P);}
    lis_finalize();
    return hard?2:accepted!=steps;
}
