/* Diagnostic: compare FGMRES, GCRO-DR and FGCRO-DR on value-changing
 * nonsymmetric sparse systems with a persistent solver per sequence.
 *
 * A is CSR with fixed sparsity but deliberately changing values, including
 * intermittent high-stiffness rows modelling an active-set transition.
 * These are synthetic profiles, NOT actual FEA/CalculiX contact matrices.
 * The change is recorded via lis_matrix_mark_values_changed before solve.
 *
 * Deliberately outside TESTS; timing never gates CI. The accompanying Python
 * checker gates convergence and parses paired warmup + measured trials.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#ifdef _COMPLEX
int main(void) { puts("SKIP: diagnostic is real-only"); return 77; }
#else
#define NSTEP 6
typedef struct {
    const char *name;
    int k,stride,intra;
} METHOD;
static const METHOD METHODS[]={
    {"fgmres",0,1,0},
    {"gcrodr",4,1,0},{"gcrodr",4,4,0},
    {"gcrodr",8,1,0},{"gcrodr",8,4,0},
    {"fgcrodr",0,1,0},
    {"fgcrodr",4,1,0},{"fgcrodr",4,4,0},
    {"fgcrodr",8,1,0},{"fgcrodr",8,4,0},
    {"fgcrodr",4,4,1}
};
#define N_METHODS ((int)(sizeof(METHODS)/sizeof(METHODS[0])))

static double wall_now(void)
{
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)t.tv_sec+(double)t.tv_usec*1e-6;
}

static int create_matrix(int side,LIS_MATRIX *out)
{
    LIS_MATRIX a=NULL;
    LIS_INT i,rc,is,ie;
    LIS_INT n=(LIS_INT)side*(LIS_INT)side;
    rc=lis_matrix_create(LIS_COMM_WORLD,&a);
    if(rc) return rc;
    if((rc=lis_matrix_set_size(a,0,n)))goto bad;
    if((rc=lis_matrix_get_range(a,&is,&ie)))goto bad;
    for(i=is;i<ie;i++)
    {
        int row=(int)(i/side),col=(int)(i%side);
        LIS_REAL d=(LIS_REAL)(4.26+0.006*(i%13));
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,a)))goto bad;
        if(row>0 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-side,-1.0,a)))goto bad;
        if(row+1<side &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+side,-1.0,a)))goto bad;
        if(col>0 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.13,a)))goto bad;
        if(col+1<side &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.87,a)))goto bad;
        /* Fixed long-range sparsity with nonsymmetric coupling. */
        if(i+side+3<n &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i+side+3,-0.095,a)))
            goto bad;
        if(i>=side+3 &&
           (rc=lis_matrix_set_value(LIS_INS_VALUE,i,i-side-3,0.060,a)))
            goto bad;
    }
    if((rc=lis_matrix_set_type(a,LIS_MATRIX_CSR)))goto bad;
    if((rc=lis_matrix_assemble(a)))goto bad;
    *out=a;
    return LIS_SUCCESS;
bad:
    if(a)lis_matrix_destroy(a);
    return rc;
}

/* Value-only changes: contact-like diagonal stiffening/release as an
 * active-set proxy. No sparsity changes or private matrix snapshots.
 * Min diagonal >=4.23, max absolute off-diagonal row sum <=4.155.
 */
static int matrix_step(LIS_MATRIX a,int step)
{
    LIS_INT row,p,is,ie;
    LIS_INT rc;
    if(a==NULL || a->matrix_type!=LIS_MATRIX_CSR ||
       a->ptr==NULL || a->index==NULL || a->value==NULL)
        return LIS_ERR_ILL_ARG;
    if((rc=lis_matrix_get_range(a,&is,&ie)))return rc;
    for(row=0;row<a->n;row++)
    {
        LIS_INT global=is+row;
        int found=0;
        int active=((global*37+step*13)%19 < 3+(step%4));
        double d=4.26+0.006*(global%13) +
                 0.06*sin((double)global*0.015+step*0.27);
        if(step>0 && active)
            d+=0.28+0.17*(step%3); /* selective stiffening */
        for(p=a->ptr[row];p<a->ptr[row+1];p++)
        {
            if(a->index[p]==row)
            {
                a->value[p]=(LIS_SCALAR)d;
                found=1;
                break;
            }
        }
        if(!found)return LIS_ERR_ILL_ARG;
    }
    return lis_matrix_mark_values_changed(a);
}

static void fill_truth(LIS_VECTOR u,int step)
{
    LIS_INT i;
    for(i=0;i<u->n;i++)
    {
        double global=(double)(u->is+i);
        u->value[i]=(LIS_SCALAR)(
           1.0+0.12*sin(0.013*global+step*0.32)
              +0.07*cos(0.029*global-step*0.19));
    }
}

static int run_sequence(LIS_MATRIX a,int side,const char *precon,
                        int restart,int trial,const METHOD *method)
{
    LIS_SOLVER solver=NULL;
    LIS_VECTOR u=NULL,b=NULL,x=NULL,ax=NULL;
    LIS_INT rc=0,status=-1,step;
    LIS_REAL nr=0.,nb=0.;
    char opts[650],detail[260],extra[170];
    int failed=1;
    if(!strcmp(method->name,"gcrodr"))
        snprintf(detail,sizeof(detail),
            " -gcrodr_keep %d -gcrodr_keep_quality 0.8 "
            "-gcrodr_reject_quality 0.2 -gcrodr_ritz_stride %d",
            method->k,method->stride);
    else if(!strcmp(method->name,"fgcrodr"))
        snprintf(detail,sizeof(detail),
            " -fgcrodr_k %d -fgcrodr_intra %d -fgcrodr_ritz_stride %d",
            method->k,method->intra,method->stride);
    else detail[0]='\0';

    if(!strcmp(precon,"ilutp"))
        snprintf(extra,sizeof(extra),
           " -ilutp 1 -ilut_rate 2 -ilut_drop 0 -ilutp_threshold 0.5");
    else extra[0]='\0';

    snprintf(opts,sizeof(opts),
        "-i %s -p %s -scale none -print none -restart %d "
        "-maxiter 500 -tol 1e-8%s%s",
        method->name,!strcmp(precon,"ilutp")?"ilut":"jacobi",
        restart,detail,extra);

    if((rc=lis_solver_create(&solver)))goto done;
    if((rc=lis_solver_set_option(opts,solver)))goto done;
    if((rc=lis_vector_duplicate(a,&u)))goto done;
    if((rc=lis_vector_duplicate(a,&b)))goto done;
    if((rc=lis_vector_duplicate(a,&x)))goto done;
    if((rc=lis_vector_duplicate(a,&ax)))goto done;

    failed=0;
    for(step=0;step<NSTEP;step++)
    {
        double w0,wall,cpu,residual=INFINITY;
        clock_t c0;
        LIS_UNSIGNED_INT builds=0,updates=0,uses=0;
        if((rc=matrix_step(a,step))) {failed=1;break;}
        fill_truth(u,step);
        if((rc=lis_matvec(a,u,b))) {failed=1;break;}
        if((rc=lis_vector_set_all(0.0,x))) {failed=1;break;}

        w0=wall_now();c0=clock();
        rc=lis_solve_sequence(a,b,x,solver);
        cpu=(double)(clock()-c0)/CLOCKS_PER_SEC;
        wall=wall_now()-w0;
        status=solver->retcode;
        if(!rc && !lis_matvec(a,x,ax) &&
           !lis_vector_xpay(b,-1.0,ax) &&
           !lis_vector_nrm2(ax,&nr) &&
           !lis_vector_nrm2(b,&nb))
            residual=nb>0.0 ? (double)(nr/nb) : (double)nr;
        if(!strcmp(method->name,"gcrodr"))
        {
            builds=(LIS_UNSIGNED_INT)solver->recycle_builds;
            updates=(LIS_UNSIGNED_INT)solver->recycle_updates;
            uses=(LIS_UNSIGNED_INT)solver->recycle_cross_uses;
        }
        else if(!strcmp(method->name,"fgcrodr"))
        {
            (void)lis_fgcrodr_get_stats(
                solver,&builds,&updates,NULL,NULL);
            (void)lis_fgcrodr_get_policy_stats(
                solver,&uses,NULL,NULL,NULL);
        }
        printf("%s,%d,%d,%d,%d,%d,%s,%d,%d,%d,"
               "%.9f,%.9f,%lld,%lld,%lld,%.12g,%llu,%llu,%llu,%.9f\n",
            precon,side,restart,trial,step,NSTEP,
            method->name,method->k,method->stride,method->intra,
            wall,cpu,(long long)solver->iter,(long long)rc,
            (long long)status,residual,
            (unsigned long long)builds,
            (unsigned long long)updates,
            (unsigned long long)uses,solver->ptime);
        fflush(stdout);
        if(rc || status!=LIS_SUCCESS ||
           !isfinite(residual) || residual>1e-8)
        {
            failed=1;
            break;
        }
    }
done:
    if(failed)fprintf(stderr,"DYNAMIC_BENCH_FAIL method=%s k=%d stride=%d "
        "intra=%d trial=%d step=%lld err=%lld status=%lld\n",
        method->name,method->k,method->stride,method->intra,
        trial,(long long)step,(long long)rc,(long long)status);
    if(ax)lis_vector_destroy(ax);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(u)lis_vector_destroy(u);
    if(solver)lis_solver_destroy(solver);
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX a=NULL;
    LIS_INT rc;
    int side=(argc>1)?atoi(argv[1]):28;
    int trials=(argc>2)?atoi(argv[2]):4;
    const char *precon=(argc>3)?argv[3]:"jacobi";
    int restart=(argc>4)?atoi(argv[4]):12;
    int trial,i,k,failed=0;
    if(side<12 || side>64 || trials<2 || trials>8 ||
       (restart!=12 && restart!=24) ||
       (strcmp(precon,"jacobi") && strcmp(precon,"ilutp")))
    {
        fprintf(stderr,"Usage: %s [side 12..64] [trials 2..8] "
                "[jacobi|ilutp] [restart 12|24]\n",argv[0]);
        return 2;
    }
    if((rc=lis_initialize(&argc,&argv)))return 3;
    if((rc=create_matrix(side,&a)))
    {
        lis_finalize();
        return 4;
    }
    puts("precon,side,restart,trial,step,nstep,method,k,stride,intra,"
         "wall_seconds,cpu_seconds,iterations,api_error,solver_status,"
         "true_relres,recycle_builds,recycle_updates,recycle_cross_uses,"
         "precon_apply_wall");
    for(trial=0;trial<trials;trial++)
    {
        for(i=0;i<N_METHODS;i++)
        {
            k=(trial%2) ? N_METHODS-i-1 : i;
            failed |= run_sequence(a,side,precon,restart,trial,&METHODS[k]);
            if(failed)break;
        }
        if(failed)break;
    }
    lis_matrix_destroy(a);
    lis_finalize();
    return failed;
}
#endif
