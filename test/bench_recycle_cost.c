/* Diagnostic-only CPU and wall-time comparison of restarted Krylov solvers.
 * Run FGMRES/GCRO-DR/FGCRO-DR on identical sparse systems, with k=0
 * controls and recycled subsequent solves. Setup input is outside timing.
 * This is not a proof of numerical correctness or a production benchmark.
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
int main(void) { puts("SKIP: real-only diagnostic"); return 77; }
#else
typedef struct {
    const char *method;
    int k,intra;
} METHOD;

static double now_seconds(void)
{
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)t.tv_sec+1e-6*(double)t.tv_usec;
}
static LIS_INT build_matrix(int side,int extended,LIS_MATRIX *out)
{
    LIS_MATRIX A=NULL;
    LIS_INT i,n=(LIS_INT)side*(LIS_INT)side,err;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err)return err;
    err=lis_matrix_set_size(A,0,n);
    if(err)goto fail;
    for(i=0;i<n;i++)
    {
        int row=(int)(i/side),col=(int)(i%side);
        LIS_REAL d=extended ? (LIS_REAL)(4.215+0.003*(i%13)) :
                              (LIS_REAL)(4.03+0.007*(i%13));
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)d,A);
        if(err)goto fail;
        if(row>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-side,-1.0,A);
            if(err)goto fail;
        }
        if(row+1<side)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+side,-1.0,A);
            if(err)goto fail;
        }
        if(col>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.18,A);
            if(err)goto fail;
        }
        if(col+1<side)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-0.82,A);
            if(err)goto fail;
        }
        if(extended)
        {
            /* Deterministic nonsymmetric long-range links.
             * Strict diagonal dominance remains: min diag 4.215,
             * sum off-diagonal magnitudes <= 4.17. */
            if(i+side+3<n)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i+side+3,-0.10,A);
                if(err)goto fail;
            }
            if(i>=side+3)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i-side-3,0.07,A);
                if(err)goto fail;
            }
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err)goto fail;
    err=lis_matrix_assemble(A);
    if(err)goto fail;
    *out=A;
    return LIS_SUCCESS;
fail:
    lis_matrix_destroy(A);
    return err;
}
static int run_method(LIS_MATRIX A,int side,int restart,int repeats,
                      const METHOD *method,const char *precon)
{
    LIS_SOLVER solver=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL,ax=NULL;
    LIS_INT err=0,status,i,n=A->n;
    char opt[512],extra[256];
    int solve;
    if(strcmp(method->method,"gcrodr")==0)
        snprintf(extra,sizeof(extra),
            " -gcrodr_keep %d -gcrodr_keep_quality 0.80 -gcrodr_reject_quality 0.20",
            method->k);
    else if(strcmp(method->method,"fgcrodr")==0)
        snprintf(extra,sizeof(extra),
            " -fgcrodr_k %d -fgcrodr_intra %d -fgcrodr_update true",
            method->k,method->intra);
    else extra[0]='\0';
    snprintf(opt,sizeof(opt),
        "-i %s -p %s -scale none -print none -restart %d "
        "-maxiter 360 -tol 1e-8%s%s",
        method->method,strcmp(precon,"ilutp")==0?"ilut":"jacobi",
        restart,extra,
        strcmp(precon,"ilutp")==0?
        " -ilutp 1 -ilut_rate 2 -ilut_drop 0 -ilutp_threshold 0.5":"");
    err=lis_solver_create(&solver);
    if(err)goto finish;
    err=lis_solver_set_option(opt,solver);
    if(err)goto finish;
    err=lis_vector_duplicate(A,&truth);
    if(err)goto finish;
    err=lis_vector_duplicate(A,&b);
    if(err)goto finish;
    err=lis_vector_duplicate(A,&x);
    if(err)goto finish;
    err=lis_vector_duplicate(A,&ax);
    if(err)goto finish;
    for(solve=0;solve<repeats;solve++)
    {
        double wall0,wall,cpu,resid=0.0;
        LIS_REAL nr=0.,nb=0.;
        clock_t cpu0;
        for(i=0;i<n;i++)
            truth->value[i]=(LIS_SCALAR)(1.0+
                0.1*sin((double)i*0.013+0.3*(double)solve)+
                0.05*cos((double)i*0.071-0.2*(double)solve));
        err=lis_matvec(A,truth,b);
        if(err)goto finish;
        err=lis_vector_set_all(0.0,x);
        if(err)goto finish;

        wall0=now_seconds();
        cpu0=clock();
        err=lis_solve(A,b,x,solver);
        cpu=(double)(clock()-cpu0)/(double)CLOCKS_PER_SEC;
        wall=now_seconds()-wall0;
        status=-1;
        if(!err) (void)lis_solver_get_status(solver,&status);

        /* Independently recomputed true residual; not part of time. */
        if(!lis_matvec(A,x,ax))
        {
            (void)lis_vector_xpay(b,-1.0,ax);
            (void)lis_vector_nrm2(ax,&nr);
            (void)lis_vector_nrm2(b,&nb);
            resid=(nb>0.) ? (double)(nr/nb) : (double)nr;
        }
        printf("%s,%d,%d,%d,%d,%d,%.9f,%.9f,%d,%d,%d,%.10g,%.9f\n",
            method->method,method->k,method->intra,restart,solve,side,
            wall,cpu,(int)solver->iter,(int)err,(int)status,resid,
            solver->ptime);
        fflush(stdout);
        if(err) {
            /* Keep the benchmark artifact even on solver failure.
             * Report the failure explicitly, then continue methods. */
            break;
        }
    }
finish:
    if(ax)lis_vector_destroy(ax);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(truth)lis_vector_destroy(truth);
    if(solver)lis_solver_destroy(solver);
    if(err)fprintf(stderr,
        "bench_recycle_cost: method=%s k=%d restart=%d last_error=%d\n",
        method->method,method->k,restart,(int)err);
    return 0;
}
int main(int argc,char **argv)
{
    const METHOD methods[]={
        {"fgmres",0,0},
        {"gcrodr",0,0},{"gcrodr",4,0},{"gcrodr",8,0},
        {"fgcrodr",0,0},{"fgcrodr",4,0},{"fgcrodr",8,0},
        {"fgcrodr",4,1}
    };
    const int restarts[]={12,24};
    LIS_MATRIX A=NULL;
    LIS_INT err;
    int side=(argc>1)?atoi(argv[1]):36;
    int repeats=(argc>2)?atoi(argv[2]):4;
    const char *precon=(argc>3)?argv[3]:"jacobi";
    const char *profile=(argc>4)?argv[4]:"standard";
    size_t j,k;
    if(side<8 || side>80 || repeats<2 || repeats>10 ||
       (strcmp(precon,"jacobi") && strcmp(precon,"ilutp")) ||
       (strcmp(profile,"standard") && strcmp(profile,"extended")))
    {
        fprintf(stderr,"Usage: %s [side 8..80] [repeats 2..10] [jacobi|ilutp] [standard|extended]\n",argv[0]);
        return 2;
    }
    err=lis_initialize(&argc,&argv);
    if(err)return 3;
    err=build_matrix(side,strcmp(profile,"extended")==0,&A);
    if(err){lis_finalize();return 4;}
    puts("method,k,intra,restart,solve_index,side,wall_seconds,cpu_seconds,iterations,api_error,solver_status,true_relres,precon_apply_wall");
    for(j=0;j<sizeof(restarts)/sizeof(restarts[0]);j++)
        for(k=0;k<sizeof(methods)/sizeof(methods[0]);k++)
            (void)run_method(A,side,restarts[j],repeats,&methods[k],precon);
    lis_matrix_destroy(A);
    lis_finalize();
    return 0;
}
#endif
