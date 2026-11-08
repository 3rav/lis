/*
 * Offline serial real-valued LIS MFAMG replay from CCX snapshots.
 * Usage: mfamgreplay A.mtx b.mtx modes.mtx [x0.mtx]
 *
 * MFAMG environment knobs are inherited from the caller. The default
 * options match the 2026-11-08 contact experiments. No output is written.
 */
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

static int read_modes(const char *path,LIS_MATRIX A,
                      LIS_INT *count_out,LIS_VECTOR **vectors_out)
{
    FILE *f=NULL;
    LIS_VECTOR *vecs=NULL;
    LIS_INT m,i,columns=0,n=A->n,err;
    long long file_n=0,file_m=0;
    char line[256];
    double value;
    int failed=1;

    *count_out=0;
    *vectors_out=NULL;
    f=fopen(path,"r");
    if(f==NULL) { perror(path); return 1; }

    if(fgets(line,sizeof(line),f)==NULL ||
       strncmp(line,"%%MatrixMarket matrix array real general",40)!=0)
    {
        fprintf(stderr,"REPLAY: modes must be a real Matrix Market array\n");
        goto done;
    }

    do
    {
        if(fgets(line,sizeof(line),f)==NULL) goto done;
    } while(line[0]=='%');

    if(sscanf(line,"%lld %lld",&file_n,&file_m)!=2 ||
       file_n!=(long long)n || file_m<=0 || file_m>256)
    {
        fprintf(stderr,
            "REPLAY: modes dimensions mismatch, matrix n=%lld modes=%lldx%lld\n",
            (long long)n,file_n,file_m);
        goto done;
    }

    columns=(LIS_INT)file_m;
    vecs=(LIS_VECTOR *)calloc((size_t)columns,sizeof(LIS_VECTOR));
    if(vecs==NULL) goto done;
    for(m=0;m<columns;m++)
    {
        err=lis_vector_duplicate(A,&vecs[m]);
        if(err) goto done;
        for(i=0;i<n;i++)
        {
            if(fscanf(f,"%lf",&value)!=1 || !isfinite(value))
            {
                fprintf(stderr,
                    "REPLAY: invalid modes value at column=%lld row=%lld\n",
                    (long long)m,(long long)i);
                goto done;
            }
            vecs[m]->value[i]=(LIS_SCALAR)value;
        }
    }
    *vectors_out=vecs;
    *count_out=columns;
    vecs=NULL;
    failed=0;
done:
    if(f) fclose(f);
    if(vecs)
    {
        for(m=0;m<columns;m++)
            if(vecs[m]) lis_vector_destroy(vecs[m]);
        free(vecs);
    }
    return failed;
}

/*
 * The generic LIS vector importer expects a LIS-specific Matrix Market
 * object: "vector" followed by indexed values. The opt-in CCX snapshot
 * exporter deliberately writes standard Matrix Market n-by-1 arrays:
 * "matrix array real general" followed by n values.
 *
 * Read that exact snapshot format directly. Do not modify source files
 * or change the semantics of the generic lis_input_vector() API.
 */
static int read_snapshot_column(const char *path,LIS_VECTOR v)
{
    FILE *f=NULL;
    char line[256];
    long long rows=0,columns=0;
    LIS_INT i;
    double value;
    int failed=1;

    if(path==NULL || v==NULL || v->n<=0 || v->value==NULL)
        return 1;
    f=fopen(path,"r");
    if(f==NULL)
    {
        perror(path);
        return 1;
    }
    if(fgets(line,sizeof(line),f)==NULL ||
       strncmp(line,"%%MatrixMarket matrix array real general",40)!=0)
    {
        fprintf(stderr,
            "REPLAY: %s: expected Matrix Market matrix array real general\n",
            path);
        goto done;
    }
    do
    {
        if(fgets(line,sizeof(line),f)==NULL)
        {
            fprintf(stderr,"REPLAY: %s: missing dimensions\n",path);
            goto done;
        }
    } while(line[0]=='%' || line[0]=='\n' || line[0]=='\r');

    if(sscanf(line,"%lld %lld",&rows,&columns)!=2 ||
       rows!=(long long)v->n || columns!=1)
    {
        fprintf(stderr,
            "REPLAY: %s: RHS dimensions %lldx%lld, expected %lldx1\n",
            path,rows,columns,(long long)v->n);
        goto done;
    }
    for(i=0;i<v->n;i++)
    {
        if(fscanf(f,"%lf",&value)!=1 || !isfinite(value))
        {
            fprintf(stderr,
                "REPLAY: %s: invalid value at row %lld\n",
                path,(long long)i+1);
            goto done;
        }
        v->value[i]=(LIS_SCALAR)value;
    }
    if(fscanf(f,"%lf",&value)==1 || ferror(f))
    {
        fprintf(stderr,"REPLAY: %s: excess or invalid trailing data\n",path);
        goto done;
    }
    failed=0;
done:
    if(fclose(f)!=0) failed=1;
    return failed;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,ax=NULL,*modes=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,nmodes=0,i,iter=0,status=-1;
    LIS_REAL reported=0.0;
    double b2=0.0,r2=0.0,t0=0.0,t1=0.0;
    int rc=1;

    if(argc!=4 && argc!=5)
    {
        fprintf(stderr,
            "Usage: %s A.mtx b.mtx modes.mtx [x0.mtx]\n",argv[0]);
        return 2;
    }
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
#if defined(_COMPLEX)
    fprintf(stderr,"REPLAY: real-valued LIS build required\n");
    goto done;
#endif
#ifdef USE_MPI
    {
        int ranks=1;
        MPI_Comm_size(LIS_COMM_WORLD,&ranks);
        if(ranks!=1)
        {
            fprintf(stderr,"REPLAY: serial LIS only\n");
            goto done;
        }
    }
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_input_matrix(A,argv[1]);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&ax);
    if(err) goto done;
    if(read_snapshot_column(argv[2],b)!=0) goto done;
    if(read_modes(argv[3],A,&nmodes,&modes)!=0) goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);
    if(err) goto done;
    if(argc==5)
    {
        if(read_snapshot_column(argv[4],x)!=0) goto done;
    }

    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i fgmres -p mfamg -scale none -tol 1.0e-5 "
        "-maxiter 1000 -restart 180 -print none",solver);
    if(err) goto done;
    err=lis_solver_set_near_nullspace(solver,nmodes,modes);
    if(err) goto done;

    printf("REPLAY n=%lld nnz=%lld modes=%lld restart=180 maxiter=1000 tol=1e-5\n",
           (long long)A->n,(long long)A->nnz,(long long)nmodes);
    t0=lis_wtime();
    err=lis_solve(A,b,x,solver);
    t1=lis_wtime();
    lis_solver_get_iter(solver,&iter);
    lis_solver_get_residualnorm(solver,&reported);
    lis_solver_get_status(solver,&status);

    if(lis_matvec(A,x,ax)==LIS_SUCCESS)
    {
        for(i=0;i<A->n;i++)
        {
            double bv=(double)b->value[i];
            double diff=bv-(double)ax->value[i];
            b2+=bv*bv;
            r2+=diff*diff;
        }
    }
    printf("REPLAY_RESULT solve_rc=%lld status=%lld iter=%lld "
           "reported_residual=%.12e true_relative=%.12e "
           "wall_seconds=%.3f\n",
           (long long)err,(long long)status,(long long)iter,
           (double)reported,b2>0 ? sqrt(r2/b2) : sqrt(r2),t1-t0);
    rc=(err==LIS_SUCCESS && status==0) ? 0 : 1;
    err=0;
done:
    if(err) fprintf(stderr,"REPLAY LIS error=%lld\n",(long long)err);
    if(solver) lis_solver_destroy(solver);
    if(modes)
    {
        for(i=0;i<nmodes;i++)
            if(modes[i]) lis_vector_destroy(modes[i]);
        free(modes);
    }
    if(ax) lis_vector_destroy(ax);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return rc;
}
