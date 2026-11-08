/* Regression: exact MFAMG near-nullspace modes are captured without A.
 * Only the serial real-valued path is exercised.
 */
#ifndef _WIN32
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#endif
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lislib.h"
#include "lis_snapshot.h"

#define N 240
#define PREFIX "mfamg-v7-regression"

static int toggle_export(int enabled)
{
#ifdef _WIN32
    return _putenv(enabled ?
        "LIS_MFAMG_MODES_EXPORT_PREFIX=" PREFIX :
        "LIS_MFAMG_MODES_EXPORT_PREFIX=");
#else
    return enabled ?
        setenv("LIS_MFAMG_MODES_EXPORT_PREFIX",PREFIX,1) :
        unsetenv("LIS_MFAMG_MODES_EXPORT_PREFIX");
#endif
}

/*
 * Match the already-established v6 multivector hierarchy regression:
 * the default pairwise path on a tiny 1-D matrix is not guaranteed to
 * coarsen a two-mode tentative space.  Force the opt-in fine target=4.
 */
static int toggle_fine_aggregation(int enabled)
{
#ifdef _WIN32
    if(_putenv(enabled ?
        "LIS_MFAMG_EXPERIMENT_FINE_AGG=1" :
        "LIS_MFAMG_EXPERIMENT_FINE_AGG=")!=0)
        return 1;
    return _putenv(enabled ?
        "LIS_MFAMG_FINE_AGG_FACTOR=2" :
        "LIS_MFAMG_FINE_AGG_FACTOR=");
#else
    if(enabled)
        return setenv("LIS_MFAMG_EXPERIMENT_FINE_AGG","1",1) ||
               setenv("LIS_MFAMG_FINE_AGG_FACTOR","2",1);
    return unsetenv("LIS_MFAMG_EXPERIMENT_FINE_AGG") ||
           unsetenv("LIS_MFAMG_FINE_AGG_FACTOR");
#endif
}

static int inspect_file(const char *path)
{
    FILE *f=fopen(path,"r");
    char line[256];
    long long n=0,m=0;
    int i,j,failed=0;
    double value,want;
    if(f==NULL) return 1;
    if(fgets(line,sizeof(line),f)==NULL ||
       strstr(line,"MatrixMarket matrix array real general")==NULL ||
       fscanf(f,"%lld %lld",&n,&m)!=2 || n!=N || m!=2)
        failed=1;
    if(!failed)
        for(j=0;j<2;j++)
            for(i=0;i<N;i++)
            {
                want=(j==0) ? 1.0 :
                    ((double)i-((double)N-1.0)/2.0)/(double)N;
                if(fscanf(f,"%lf",&value)!=1 ||
                   fabs(value-want)>1.0e-12)
                {
                    failed=1;
                    break;
                }
            }
    if(fclose(f)!=0) failed=1;
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR modes[2]={NULL,NULL};
    LIS_SOLVER solver=NULL;
    LIS_PRECON precon=NULL;
    LIS_INT err=0,i;
    char base[256],path[288],other[288];
    FILE *f=NULL;
    int failed=1,export_enabled=0,fine_enabled=0;

    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#ifdef _COMPLEX
    lis_finalize();
    return 0;
#endif
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1) { lis_finalize(); return 0; }
    }
#endif

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,N);
    if(err) goto done;
    for(i=0;i<N;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
            if(err) goto done;
        }
        if(i+1<N)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_vector_duplicate(A,&modes[0]);
    if(err) goto done;
    err=lis_vector_duplicate(A,&modes[1]);
    if(err) goto done;
    for(i=0;i<N;i++)
    {
        modes[0]->value[i]=1.0;
        modes[1]->value[i]=
            ((double)i-((double)N-1.0)/2.0)/(double)N;
    }

    if(snprintf(base,sizeof(base),PREFIX "_000001_n%d_nnz%lld",
                N,(long long)A->nnz)>=sizeof(base))
        goto done;
    snprintf(path,sizeof(path),"%s_modes.mtx",base);
    snprintf(other,sizeof(other),"%s_A.mtx",base);
    /* Clean only this narrowly scoped test's previously interrupted run. */
    remove(path);
    if(toggle_export(1)!=0) goto done;
    export_enabled=1;
    if(toggle_fine_aggregation(1)!=0) goto done;
    fine_enabled=1;

    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i gmres -p mfamg -scale none -print none",solver);
    if(err) goto done;
    err=lis_solver_set_near_nullspace(solver,2,modes);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;

    err=lis_precon_create(solver,&precon);
    if(err) { precon=NULL; goto done; }

    if(inspect_file(path)!=0)
    {
        fprintf(stderr,"MFAMG V7: exported basis mismatch\n");
        goto done;
    }
    f=fopen(other,"r");
    if(f!=NULL)
    {
        fclose(f); f=NULL;
        fprintf(stderr,"MFAMG V7: unexpected matrix export\n");
        goto done;
    }

    /* Reuse of destination must never overwrite existing data. */
    err=lis_snapshot_export_modes_mm(2,modes,base);
    if(err!=LIS_ERR_FILE_IO || inspect_file(path)!=0)
    {
        fprintf(stderr,"MFAMG V7: destination collision not rejected\n");
        goto done;
    }
    err=0;
    printf("MFAMG V7 modes-only export and no-overwrite PASS\n");
    failed=0;

done:
    if(err) fprintf(stderr,"MFAMG V7 test LIS error=%d\n",(int)err);
    if(precon) lis_precon_destroy(precon);
    if(solver) lis_solver_destroy(solver);
    if(modes[1]) lis_vector_destroy(modes[1]);
    if(modes[0]) lis_vector_destroy(modes[0]);
    if(A) lis_matrix_destroy(A);
    if(fine_enabled && toggle_fine_aggregation(0)!=0) failed=1;
    if(export_enabled && toggle_export(0)!=0) failed=1;
    if(export_enabled) remove(path);
    lis_finalize();
    return failed;
}
