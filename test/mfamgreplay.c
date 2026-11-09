/*
 * Offline serial real-valued LIS MFAMG replay from CCX snapshots.
 * Usage: mfamgreplay A.mtx b.mtx modes.mtx [x0.mtx]
 *
 * MFAMG environment knobs are inherited from the caller. The default
 * FGMRES options match the 2026-11-08 contact experiments.
 *
 * Opt-in offline recycling experiment:
 *   LIS_MFAMG_REPLAY_METHOD=fgcrodr
 *   LIS_MFAMG_REPLAY_RECYCLE_DIM=8  (allowed 0..64, less than restart)
 *
 * FGCRO-DR maintains an A*Z=C recycle relation while allowing flexible
 * right preconditioning. No solver or MFAMG production path is modified.
 * No output files are written by this program.
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
#include <errno.h>
#include "lislib.h"
#include "lis_solver.h"

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

/*
 * Offline-only tuning. The original restart=180 and maxiter=1000
 * remain unchanged unless the caller explicitly selects an override.
 * The restart bound limits vector allocation; the maxiter bound limits
 * runaway work on difficult contact snapshots.
 */
static int replay_int_option(const char *name,int fallback,
                             int lower,int upper,int *result)
{
    const char *raw=getenv(name);
    char *end=NULL;
    long value;
    if(raw==NULL)
    {
        *result=fallback;
        return 0;
    }
    errno=0;
    value=strtol(raw,&end,10);
    if(errno!=0 || end==raw || *end!='\0' ||
       value<lower || value>upper)
    {
        fprintf(stderr,
            "REPLAY: invalid %s=%s (expected integer %d..%d)\n",
            name,raw,lower,upper);
        return 1;
    }
    *result=(int)value;
    return 0;
}


/*
 * Same-structure, different-values replay ONLY; never alter a sparsity
 * graph or blindly assume two CCX contact assemblies have the same graph.
 *
 * Preflight every row pointer, every column index and every value before
 * changing the destination.  This function does not allocate or expose
 * contact data. It accepts only owned, assembled, serial real CSR arrays.
 * The second imported matrix remains an INDEPENDENT A2 for the residual.
 */
static int replay_check_or_replace_csr_values(LIS_MATRIX target,LIS_MATRIX source,
                                              LIS_INT *changed_out,int apply)
{
    LIS_INT i,changed=0;

    *changed_out=0;
    if(target==NULL || source==NULL ||
       target->matrix_type!=LIS_MATRIX_CSR ||
       source->matrix_type!=LIS_MATRIX_CSR ||
       target->is_splited || source->is_splited ||
       target->is_scaled || source->is_scaled ||
       target->csr_external || source->csr_external ||
       target->n<=0 || source->n!=target->n ||
       target->np!=source->np || target->gn!=source->gn ||
       target->nnz<=0 || source->nnz!=target->nnz ||
       target->ptr==NULL || source->ptr==NULL ||
       target->index==NULL || source->index==NULL ||
       target->value==NULL || source->value==NULL)
    {
        fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values requires identical owned unsplit real CSR shape\n");
        return 1;
    }
    if(target->ptr[0]!=0 || source->ptr[0]!=0 ||
       target->ptr[target->n]!=target->nnz ||
       source->ptr[source->n]!=source->nnz)
    {
        fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values invalid CSR row pointers\n");
        return 1;
    }
    for(i=0;i<target->n;i++)
    {
        if(target->ptr[i]!=source->ptr[i] ||
           target->ptr[i+1]!=source->ptr[i+1] ||
           target->ptr[i]<0 || target->ptr[i+1]<target->ptr[i] ||
           target->ptr[i+1]>target->nnz)
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values CSR row graph mismatch row=%lld\n",
                    (long long)(i+1));
            return 1;
        }
    }
    for(i=0;i<target->nnz;i++)
    {
        if(target->index[i]!=source->index[i] ||
           target->index[i]<0 || target->index[i]>=target->gn)
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values CSR column graph mismatch offset=%lld\n",
                    (long long)i);
            return 1;
        }
        if(!isfinite((double)target->value[i]) ||
           !isfinite((double)source->value[i]))
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values nonfinite CSR value offset=%lld\n",
                    (long long)i);
            return 1;
        }
        if(target->value[i]!=source->value[i]) changed++;
    }
    if(changed==0)
    {
        fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values requires a genuinely changed operator\n");
        return 1;
    }
    /* Dry-run option deliberately never modifies the destination. */
    if(apply)
    {
        for(i=0;i<target->nnz;i++)
            target->value[i]=source->value[i];
    }
    *changed_out=changed;
    return 0;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,A2=NULL,A3=NULL;
    LIS_VECTOR b=NULL,x=NULL,ax=NULL,b2=NULL,x2=NULL,ax2=NULL;
    LIS_VECTOR b3=NULL,x3=NULL,ax3=NULL,*modes=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,nmodes=0,i,iter=0,status=-1;
    LIS_REAL reported=0.0;
    LIS_REAL recycle_orth=0.0,recycle_relation=0.0;
    LIS_UNSIGNED_INT recycle_builds=0,recycle_updates=0;
    LIS_UNSIGNED_INT recycle_uses=0,recycle_rejects=0;
    double b_norm=0.0,r_norm=0.0,true_relative=0.0,t0=0.0,t1=0.0;
    int rc=1,restart=180,maxiter=1000,recycle_k=8,use_recycling=0;
    int recycle_validation_failed=0;
    char *next_matrix=getenv("LIS_MFAMG_REPLAY_NEXT_A");
    const char *next_rhs=getenv("LIS_MFAMG_REPLAY_NEXT_B");
    const char *next_diag_shift=getenv("LIS_MFAMG_REPLAY_NEXT_DIAG_SHIFT");
    const char *next_replace_values=getenv("LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES");
    char *third_matrix=getenv("LIS_MFAMG_REPLAY_THIRD_A");
    const char *third_rhs=getenv("LIS_MFAMG_REPLAY_THIRD_B");
    const char *third_replace_values=getenv("LIS_MFAMG_REPLAY_THIRD_REPLACE_VALUES");
    const char *method="fgmres";
    const char *method_env;
    LIS_INT recycle_dim=0,recycle_action=0;
    LIS_REAL recycle_quality=0.0,recycle_quality_min=0.0;
    int optlen;
    char solver_options[256];

    /*
     * Exact dry-run preflight, no modes/RHS/solver needed. Compare the
     * actual LIS-imported CSR graph and finite values using precisely
     * the same gate as the opt-in in-place value replacement. Do not
     * pass --compare-csr through the LIS option parser.
     */
    if(argc>1 && strcmp(argv[1],"--compare-csr")==0)
    {
        int lis_argc=1;
        char *lis_args_storage[2]={argv[0],NULL};
        char **lis_args=lis_args_storage;
        LIS_INT changed=0;

        if(argc!=4)
        {
            fprintf(stderr,"Usage: %s --compare-csr A1.mtx A2.mtx\n",argv[0]);
            return 2;
        }
        if(lis_initialize(&lis_argc,&lis_args)!=LIS_SUCCESS) return 1;
#if defined(_COMPLEX)
        fprintf(stderr,"REPLAY_CSR_COMPARE_ERROR real-valued LIS required\n");
        goto done;
#endif
#ifdef USE_MPI
        {
            int ranks=1;
            MPI_Comm_size(LIS_COMM_WORLD,&ranks);
            if(ranks!=1)
            {
                fprintf(stderr,"REPLAY_CSR_COMPARE_ERROR serial LIS only\n");
                goto done;
            }
        }
#endif
        err=lis_matrix_create(LIS_COMM_WORLD,&A);
        if(err) goto done;
        err=lis_input_matrix(A,argv[2]);
        if(err) goto done;
        err=lis_matrix_create(LIS_COMM_WORLD,&A2);
        if(err) goto done;
        err=lis_input_matrix(A2,argv[3]);
        if(err) goto done;
        if(replay_check_or_replace_csr_values(A,A2,&changed,0)!=0)
        {
            fprintf(stderr,"REPLAY_CSR_COMPARE compatible=0\n");
            goto done;
        }
        printf("REPLAY_CSR_COMPARE compatible=1 n=%lld nnz=%lld changed_nnz=%lld\n",
               (long long)A->n,(long long)A->nnz,(long long)changed);
        rc=0;
        err=0;
        goto done;
    }

    if((third_matrix==NULL)!=(third_rhs==NULL) ||
       (third_matrix!=NULL && next_matrix==NULL) ||
       (third_replace_values!=NULL &&
        (strcmp(third_replace_values,"1")!=0 || third_matrix==NULL)))
    {
        fprintf(stderr,"REPLAY_SEQUENCE_ERROR third step requires NEXT_A/B and THIRD_A/B, replace flag=1\n");
        return 2;
    }
    if(next_replace_values!=NULL &&
       (strcmp(next_replace_values,"1")!=0 ||
        next_matrix==NULL || next_rhs==NULL || next_diag_shift!=NULL))
    {
        fprintf(stderr,"REPLAY_NEWTON_ERROR replace-values requires NEXT_A and NEXT_B, no diagonal shift, flag=1\n");
        return 2;
    }
    if(next_diag_shift!=NULL && (next_matrix==NULL || next_rhs==NULL))
    {
        fprintf(stderr,"REPLAY_NEWTON_ERROR diagonal shift requires NEXT_A and NEXT_B\n");
        return 2;
    }
    if((next_matrix==NULL)!=(next_rhs==NULL))
    {
        fprintf(stderr,"REPLAY_NEWTON_ERROR specify both NEXT_A and NEXT_B\n");
        return 2;
    }
    /* Preserve every completed step even if a later step crashes. */
    if(next_matrix!=NULL)
        setvbuf(stdout,NULL,_IONBF,0);
    if(argc!=4 && argc!=5)
    {
        fprintf(stderr,
            "Usage: %s A.mtx b.mtx modes.mtx [x0.mtx]\n",argv[0]);
        return 2;
    }
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 1;
    if(replay_int_option("LIS_MFAMG_REPLAY_RESTART",180,2,600,
                         &restart)!=0 ||
       replay_int_option("LIS_MFAMG_REPLAY_MAXITER",1000,100,6000,
                         &maxiter)!=0)
        goto done;

    method_env=getenv("LIS_MFAMG_REPLAY_METHOD");
    if(method_env!=NULL)
    {
        if(strcmp(method_env,"fgmres")==0)
            method="fgmres";
        else if(strcmp(method_env,"fgcrodr")==0)
        {
            method="fgcrodr";
            use_recycling=1;
        }
        else
        {
            fprintf(stderr,
                "REPLAY: LIS_MFAMG_REPLAY_METHOD must be fgmres or fgcrodr\n");
            goto done;
        }
    }
    if(use_recycling)
    {
        if(replay_int_option("LIS_MFAMG_REPLAY_RECYCLE_DIM",8,0,64,
                             &recycle_k)!=0)
            goto done;
        if(recycle_k>=restart)
        {
            fprintf(stderr,
                "REPLAY: recycle dimension must be less than restart\n");
            goto done;
        }
    }
    else if(getenv("LIS_MFAMG_REPLAY_RECYCLE_DIM")!=NULL)
    {
        fprintf(stderr,
            "REPLAY: recycle dimension requires LIS_MFAMG_REPLAY_METHOD=fgcrodr\n");
        goto done;
    }
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
    optlen=snprintf(solver_options,sizeof(solver_options),
        "-i %s -p mfamg -scale none -tol 1.0e-5 "
        "-maxiter %d -restart %d -print none",
        method,maxiter,restart);
    if(optlen<0 || (size_t)optlen>=sizeof(solver_options))
    {
        fprintf(stderr,"REPLAY: solver options too long\n");
        goto done;
    }
    err=lis_solver_set_option(solver_options,solver);
    if(err) goto done;
    if(use_recycling)
    {
        optlen=snprintf(solver_options,sizeof(solver_options),
            "-fgcrodr_k %d",recycle_k);
        if(optlen<0 || (size_t)optlen>=sizeof(solver_options))
            goto done;
        err=lis_solver_set_option(solver_options,solver);
        if(err) goto done;
    }
    err=lis_solver_set_near_nullspace(solver,nmodes,modes);
    if(err) goto done;

    printf("REPLAY n=%lld nnz=%lld modes=%lld restart=%d maxiter=%d tol=1e-5 method=%s recycle_k=%d\n",
           (long long)A->n,(long long)A->nnz,(long long)nmodes,
           restart,maxiter,method,use_recycling?recycle_k:0);
    t0=lis_wtime();
    err=lis_solve(A,b,x,solver);
    t1=lis_wtime();
    lis_solver_get_iter(solver,&iter);
    lis_solver_get_residualnorm(solver,&reported);
    lis_solver_get_status(solver,&status);
    if(use_recycling)
    {
        LIS_INT recycle_rc=lis_solver_get_fgcrodr_recycle(
            solver,&recycle_dim,&recycle_action,
            &recycle_quality,&recycle_quality_min);
        if(recycle_rc==LIS_SUCCESS)
            printf("REPLAY_RECYCLE dim=%lld action=%lld quality=%.9e min_quality=%.9e\n",
                   (long long)recycle_dim,(long long)recycle_action,
                   (double)recycle_quality,(double)recycle_quality_min);
        else
            fprintf(stderr,"REPLAY: FGCRO-DR recycle diagnostic rc=%lld\n",
                    (long long)recycle_rc);

        recycle_rc=lis_fgcrodr_get_stats(solver,
            &recycle_builds,&recycle_updates,
            &recycle_uses,&recycle_rejects);
        if(recycle_rc==LIS_SUCCESS)
            printf("REPLAY_RECYCLE_STATS builds=%llu updates=%llu "
                   "uses=%llu rejects=%llu\n",
                   (unsigned long long)recycle_builds,
                   (unsigned long long)recycle_updates,
                   (unsigned long long)recycle_uses,
                   (unsigned long long)recycle_rejects);

        /* Expensive A*Z=C check is opt-in, never alters the result. */
        {
            const char *validate=getenv("LIS_MFAMG_REPLAY_VALIDATE_RECYCLE");
            if(validate!=NULL && strcmp(validate,"1")==0)
            {
                LIS_INT validation_rc=lis_fgcrodr_validate_recycle(
                    solver,&recycle_orth,&recycle_relation);
                printf("REPLAY_RECYCLE_CHECK rc=%lld orth_error=%.12e "
                       "relation_error=%.12e\n",
                       (long long)validation_rc,
                       (double)recycle_orth,(double)recycle_relation);
                if(validation_rc!=LIS_SUCCESS ||
                   !isfinite((double)recycle_orth) ||
                   !isfinite((double)recycle_relation) ||
                   recycle_orth>1.0e-8 || recycle_relation>1.0e-8)
                    recycle_validation_failed=1;
            }
        }
    }

    /*
     * Fail closed: a failed explicit A*x must NEVER turn an
     * uninitialized residual accumulator into apparent convergence.
     * hypot() avoids overflow when squaring large contact RHS entries.
     */
    {
        LIS_INT matvec_rc=lis_matvec(A,x,ax);
        if(matvec_rc!=LIS_SUCCESS)
        {
            fprintf(stderr,"REPLAY_TRUE_RESIDUAL_ERROR matvec_rc=%lld\n",
                    (long long)matvec_rc);
            err=matvec_rc;
            goto done;
        }
    }
    for(i=0;i<A->n;i++)
    {
        const double bv=(double)b->value[i];
        const double axv=(double)ax->value[i];
        const double diff=bv-axv;
        if(!isfinite(bv) || !isfinite(axv) || !isfinite(diff))
        {
            fprintf(stderr,"REPLAY_TRUE_RESIDUAL_ERROR nonfinite row=%lld\n",
                    (long long)(i+1));
            rc=1;
            goto done;
        }
        b_norm=hypot(b_norm,bv);
        r_norm=hypot(r_norm,diff);
    }
    true_relative=b_norm>0.0 ? r_norm/b_norm : r_norm;
    if(!isfinite(b_norm) || !isfinite(r_norm) ||
       !isfinite(true_relative) || !isfinite((double)reported))
    {
        fprintf(stderr,"REPLAY_TRUE_RESIDUAL_ERROR nonfinite norm\n");
        rc=1;
        goto done;
    }
    printf("REPLAY_RESULT solve_rc=%lld status=%lld iter=%lld "
           "reported_residual=%.12e true_relative=%.12e "
           "wall_seconds=%.3f\n",
           (long long)err,(long long)status,(long long)iter,
           (double)reported,true_relative,t1-t0);
    rc=(err==LIS_SUCCESS && status==0 &&
        true_relative<=1.0e-5 &&
        !recycle_validation_failed) ? 0 : 1;
    /*
     * Optional two-operator in-process experiment: retain exactly the
     * same LIS_SOLVER and MFAMG baseline modes while changing A and b.
     * Each step starts from zero, like a Newton linear correction.  This
     * is not a complete nonlinear Newton/contact trajectory.
     */
    if(next_matrix!=NULL)
    {
        LIS_INT err2,iter2=0,status2=-1;
        LIS_MATRIX solve_matrix=NULL;
        LIS_UNSIGNED_INT before_cross=0,before_keep=0;
        LIS_UNSIGNED_INT before_transport=0,before_reject=0;
        LIS_UNSIGNED_INT after_cross=0,after_keep=0;
        LIS_UNSIGNED_INT after_transport=0,after_reject=0;
        LIS_INT policies_available=LIS_FALSE;
        LIS_REAL reported2=0.0,orth2=0.0,relation2=0.0;
        LIS_UNSIGNED_INT builds2=0,updates2=0,uses2=0,rejects2=0;
        double bnorm2=0.0,rnorm2=0.0,true2=0.0;
        double start2,end2;
        int second_ok=1;
        int first_ok=(rc==0);

        /* Default to failure before any second-step allocation or I/O:
         * early goto done must never reuse a successful step-1 rc. */
        rc=1;
        err2=lis_matrix_create(LIS_COMM_WORLD,&A2);
        if(err2) { err=err2; goto done; }
        err2=lis_input_matrix(A2,next_matrix);
        if(err2) { err=err2; goto done; }
        /* A2 was NULL at declaration; select it only once initialized. */
        solve_matrix=A2;
        if(A2->n!=A->n || A2->np!=A->np)
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR operator dimensions/distribution changed\n");
            goto done;
        }
        err2=lis_vector_duplicate(A2,&b2);
        if(err2) { err=err2; goto done; }
        err2=lis_vector_duplicate(A2,&x2);
        if(err2) { err=err2; goto done; }
        err2=lis_vector_duplicate(A2,&ax2);
        if(err2) { err=err2; goto done; }
        if(read_snapshot_column(next_rhs,b2)!=0) goto done;
        err2=lis_vector_set_all((LIS_SCALAR)0.0,x2);
        if(err2) { err=err2; goto done; }

        if(use_recycling)
        {
            err2=lis_fgcrodr_get_policy_stats(
                solver,&before_cross,&before_keep,
                &before_transport,&before_reject);
            if(err2!=LIS_SUCCESS) { err=err2; goto done; }
            policies_available=LIS_TRUE;
        }
        if(next_replace_values!=NULL)
        {
            LIS_INT changed=0;
            if(replay_check_or_replace_csr_values(A,A2,&changed,1)!=0)
                goto done;
            err2=lis_matrix_mark_values_changed(A);
            if(err2) { err=err2; goto done; }
            solve_matrix=A;
            printf("REPLAY_NEWTON_MATRIX mode=inplace_replace_values changed_nnz=%lld\n",
                   (long long)changed);
        }
        else if(next_diag_shift!=NULL)
        {
            char *end=NULL;
            double shift;
            errno=0;
            shift=strtod(next_diag_shift,&end);
            if(errno!=0 || end==next_diag_shift ||
               *end!='\0' || !isfinite(shift))
            {
                fprintf(stderr,"REPLAY_NEWTON_ERROR invalid diagonal shift\n");
                goto done;
            }
            /*
             * Controlled same-object / same-structure experiment:
             * the caller must supply the independently generated A2
             * for true-residual verification.  A has NOT changed
             * structure: notify the recycle policy of changed values.
             */
            err2=lis_matrix_shift_diagonal(A,(LIS_SCALAR)shift);
            if(err2) { err=err2; goto done; }
            err2=lis_matrix_mark_values_changed(A);
            if(err2) { err=err2; goto done; }
            solve_matrix=A;
            printf("REPLAY_NEWTON_MATRIX mode=inplace_diag_shift shift=%.12g\n",
                   shift);
        }
        else
        {
            printf("REPLAY_NEWTON_MATRIX mode=separate_object\n");
        }

        printf("REPLAY_NEWTON_STEP step=2 n=%lld nnz=%lld modes=%lld method=%s recycle_k=%d\n",
               (long long)A2->n,(long long)A2->nnz,(long long)nmodes,
               method,use_recycling?recycle_k:0);
        start2=lis_wtime();
        err2=lis_solve(solve_matrix,b2,x2,solver);
        end2=lis_wtime();
        if(lis_solver_get_iter(solver,&iter2)!=LIS_SUCCESS ||
           lis_solver_get_residualnorm(solver,&reported2)!=LIS_SUCCESS ||
           lis_solver_get_status(solver,&status2)!=LIS_SUCCESS)
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR missing second-solve diagnostics\n");
            goto done;
        }
        if(use_recycling)
        {
            LIS_INT vrc=lis_fgcrodr_get_stats(
                solver,&builds2,&updates2,&uses2,&rejects2);
            if(vrc!=LIS_SUCCESS)
            {
                fprintf(stderr,"REPLAY_NEWTON_ERROR missing recycle statistics\n");
                second_ok=0;
            }
            printf("REPLAY_NEWTON_RECYCLE_STATS builds=%llu updates=%llu uses=%llu rejects=%llu\n",
                   (unsigned long long)builds2,(unsigned long long)updates2,
                   (unsigned long long)uses2,(unsigned long long)rejects2);
            if(policies_available)
            {
                vrc=lis_fgcrodr_get_policy_stats(
                    solver,&after_cross,&after_keep,
                    &after_transport,&after_reject);
                if(vrc!=LIS_SUCCESS ||
                   after_cross<before_cross ||
                   after_keep<before_keep ||
                   after_transport<before_transport ||
                   after_reject<before_reject)
                {
                    fprintf(stderr,"REPLAY_NEWTON_ERROR invalid policy counters\n");
                    second_ok=0;
                }
                else
                {
                    printf("REPLAY_NEWTON_POLICY cross_uses_delta=%llu keeps_delta=%llu transports_delta=%llu rejects_delta=%llu\n",
                           (unsigned long long)(after_cross-before_cross),
                           (unsigned long long)(after_keep-before_keep),
                           (unsigned long long)(after_transport-before_transport),
                           (unsigned long long)(after_reject-before_reject));
                }
            }
            if(getenv("LIS_MFAMG_REPLAY_VALIDATE_RECYCLE")!=NULL &&
               strcmp(getenv("LIS_MFAMG_REPLAY_VALIDATE_RECYCLE"),"1")==0)
            {
                vrc=lis_fgcrodr_validate_recycle(solver,&orth2,&relation2);
                printf("REPLAY_NEWTON_RECYCLE_CHECK rc=%lld orth_error=%.12e relation_error=%.12e\n",
                       (long long)vrc,(double)orth2,(double)relation2);
                if(vrc!=LIS_SUCCESS || !isfinite((double)orth2) ||
                   !isfinite((double)relation2) || orth2>1.0e-8 ||
                   relation2>1.0e-8)
                    second_ok=0;
            }
        }
        if(lis_matvec(A2,x2,ax2)!=LIS_SUCCESS)
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR explicit A2*x2 failed\n");
            goto done;
        }
        /*
         * A2 is an independent externally supplied control operator.
         * For the in-place mutation mode, explicitly verify the solution
         * against A2; success must not rely on a caller's shift claim.
         */

        for(i=0;i<A2->n;i++)
        {
            const double bv=(double)b2->value[i];
            const double av=(double)ax2->value[i];
            const double diff=bv-av;
            if(!isfinite(bv) || !isfinite(av) || !isfinite(diff))
            {
                fprintf(stderr,"REPLAY_NEWTON_ERROR nonfinite residual row=%lld\n",
                        (long long)(i+1));
                goto done;
            }
            bnorm2=hypot(bnorm2,bv);
            rnorm2=hypot(rnorm2,diff);
        }
        true2=bnorm2>0.0 ? rnorm2/bnorm2 : rnorm2;
        if(!isfinite(true2) || !isfinite((double)reported2))
        {
            fprintf(stderr,"REPLAY_NEWTON_ERROR nonfinite second residual norm\n");
            goto done;
        }
        printf("REPLAY_NEWTON_RESULT step=2 solve_rc=%lld status=%lld iter=%lld "
               "reported_residual=%.12e true_relative=%.12e wall_seconds=%.3f\n",
               (long long)err2,(long long)status2,(long long)iter2,
               (double)reported2,true2,end2-start2);
        if(first_ok && err2==LIS_SUCCESS && status2==0 &&
           true2<=1.0e-5 && second_ok)
            rc=0;

        /*
         * Optional third consecutive linear system, using the SAME solver
         * object retained by steps 1 and 2. No new preconditioner policy:
         * the solver handles KEEP/TRANSPORT/REJECT upon each lis_solve.
         *
         * For fixed CSR, mutate the matrix object actually used by step 2.
         * For changing graph, pass a distinct A3 to lis_solve and allow
         * the solver to reject stale recycle directions.
         */
        if(third_matrix!=NULL)
        {
            LIS_MATRIX solve_matrix3=NULL;
            LIS_INT err3,iter3=0,status3=-1,changed3=0;
            LIS_REAL reported3=0.0,orth3=0.0,relation3=0.0;
            LIS_UNSIGNED_INT cross0=0,keep0=0,transport0=0,reject0=0;
            LIS_UNSIGNED_INT cross1=0,keep1=0,transport1=0,reject1=0;
            LIS_UNSIGNED_INT builds3=0,updates3=0,uses3=0,rejects3=0;
            double bnorm3=0.0,rnorm3=0.0,true3=0.0;
            double start3,end3;
            int previous_ok=(rc==0),last_ok=1;

            rc=1; /* Any early third-step failure must fail closed. */
            err3=lis_matrix_create(LIS_COMM_WORLD,&A3);
            if(err3) { err=err3; goto done; }
            err3=lis_input_matrix(A3,third_matrix);
            if(err3) { err=err3; goto done; }
            if(A3->n!=A->n || A3->np!=A->np)
            {
                fprintf(stderr,"REPLAY_SEQUENCE_ERROR third operator dimensions changed\n");
                goto done;
            }
            err3=lis_vector_duplicate(A3,&b3);
            if(err3) { err=err3; goto done; }
            err3=lis_vector_duplicate(A3,&x3);
            if(err3) { err=err3; goto done; }
            err3=lis_vector_duplicate(A3,&ax3);
            if(err3) { err=err3; goto done; }
            if(read_snapshot_column(third_rhs,b3)!=0) goto done;
            err3=lis_vector_set_all((LIS_SCALAR)0.0,x3);
            if(err3) { err=err3; goto done; }

            if(use_recycling)
            {
                err3=lis_fgcrodr_get_policy_stats(
                    solver,&cross0,&keep0,&transport0,&reject0);
                if(err3!=LIS_SUCCESS) { err=err3; goto done; }
            }
            if(third_replace_values!=NULL)
            {
                if(replay_check_or_replace_csr_values(
                       solve_matrix,A3,&changed3,1)!=0)
                {
                    fprintf(stderr,"REPLAY_SEQUENCE_ERROR third CSR preflight rejected\n");
                    goto done;
                }
                err3=lis_matrix_mark_values_changed(solve_matrix);
                if(err3) { err=err3; goto done; }
                solve_matrix3=solve_matrix;
                printf("REPLAY_SEQUENCE_MATRIX step=3 mode=inplace_replace_values changed_nnz=%lld\n",
                       (long long)changed3);
            }
            else
            {
                solve_matrix3=A3;
                printf("REPLAY_SEQUENCE_MATRIX step=3 mode=separate_object\n");
            }
            printf("REPLAY_SEQUENCE_STEP step=3 n=%lld nnz=%lld modes=%lld method=%s recycle_k=%d\n",
                   (long long)A3->n,(long long)A3->nnz,
                   (long long)nmodes,method,use_recycling?recycle_k:0);
            start3=lis_wtime();
            err3=lis_solve(solve_matrix3,b3,x3,solver);
            end3=lis_wtime();
            if(lis_solver_get_iter(solver,&iter3)!=LIS_SUCCESS ||
               lis_solver_get_residualnorm(solver,&reported3)!=LIS_SUCCESS ||
               lis_solver_get_status(solver,&status3)!=LIS_SUCCESS)
            {
                fprintf(stderr,"REPLAY_SEQUENCE_ERROR missing third-step diagnostics\n");
                goto done;
            }
            if(use_recycling)
            {
                LIS_INT vrc=lis_fgcrodr_get_stats(
                    solver,&builds3,&updates3,&uses3,&rejects3);
                if(vrc!=LIS_SUCCESS) last_ok=0;
                printf("REPLAY_SEQUENCE_RECYCLE_STATS step=3 builds=%llu updates=%llu uses=%llu rejects=%llu\n",
                       (unsigned long long)builds3,
                       (unsigned long long)updates3,
                       (unsigned long long)uses3,
                       (unsigned long long)rejects3);
                vrc=lis_fgcrodr_get_policy_stats(
                    solver,&cross1,&keep1,&transport1,&reject1);
                if(vrc!=LIS_SUCCESS || cross1<cross0 || keep1<keep0 ||
                   transport1<transport0 || reject1<reject0)
                {
                    fprintf(stderr,"REPLAY_SEQUENCE_ERROR invalid third policy counters\n");
                    last_ok=0;
                }
                else
                {
                    printf("REPLAY_SEQUENCE_POLICY step=3 cross_uses_delta=%llu keeps_delta=%llu transports_delta=%llu rejects_delta=%llu\n",
                           (unsigned long long)(cross1-cross0),
                           (unsigned long long)(keep1-keep0),
                           (unsigned long long)(transport1-transport0),
                           (unsigned long long)(reject1-reject0));
                }
                if(getenv("LIS_MFAMG_REPLAY_VALIDATE_RECYCLE")!=NULL &&
                   strcmp(getenv("LIS_MFAMG_REPLAY_VALIDATE_RECYCLE"),"1")==0)
                {
                    vrc=lis_fgcrodr_validate_recycle(
                        solver,&orth3,&relation3);
                    printf("REPLAY_SEQUENCE_RECYCLE_CHECK step=3 rc=%lld orth_error=%.12e relation_error=%.12e\n",
                           (long long)vrc,(double)orth3,(double)relation3);
                    if(vrc!=LIS_SUCCESS ||
                       !isfinite((double)orth3) ||
                       !isfinite((double)relation3) ||
                       orth3>1.0e-8 || relation3>1.0e-8)
                        last_ok=0;
                }
            }
            /* Independently imported A3 is the source of truth,
             * including when solve_matrix3 is an in-place mutation. */
            if(lis_matvec(A3,x3,ax3)!=LIS_SUCCESS)
            {
                fprintf(stderr,"REPLAY_SEQUENCE_ERROR explicit third A3*x failed\n");
                goto done;
            }
            for(i=0;i<A3->n;i++)
            {
                const double bv=(double)b3->value[i];
                const double av=(double)ax3->value[i];
                const double diff=bv-av;
                if(!isfinite(bv) || !isfinite(av) || !isfinite(diff))
                {
                    fprintf(stderr,"REPLAY_SEQUENCE_ERROR third nonfinite row=%lld\n",
                            (long long)(i+1));
                    goto done;
                }
                bnorm3=hypot(bnorm3,bv);
                rnorm3=hypot(rnorm3,diff);
            }
            true3=bnorm3>0.0 ? rnorm3/bnorm3 : rnorm3;
            if(!isfinite(true3) || !isfinite((double)reported3))
            {
                fprintf(stderr,"REPLAY_SEQUENCE_ERROR nonfinite third residual\n");
                goto done;
            }
            printf("REPLAY_SEQUENCE_RESULT step=3 solve_rc=%lld status=%lld iter=%lld reported_residual=%.12e true_relative=%.12e wall_seconds=%.3f\n",
                   (long long)err3,(long long)status3,(long long)iter3,
                   (double)reported3,true3,end3-start3);
            if(previous_ok && err3==LIS_SUCCESS && status3==0 &&
               true3<=1.0e-5 && last_ok)
                rc=0;
        }
    }
    err=0;
done:
    if(err) fprintf(stderr,"REPLAY LIS error=%lld\n",(long long)err);
    if(solver) lis_solver_destroy(solver);
    if(ax3) lis_vector_destroy(ax3);
    if(x3) lis_vector_destroy(x3);
    if(b3) lis_vector_destroy(b3);
    if(A3) lis_matrix_destroy(A3);
    if(ax2) lis_vector_destroy(ax2);
    if(x2) lis_vector_destroy(x2);
    if(b2) lis_vector_destroy(b2);
    if(A2) lis_matrix_destroy(A2);
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
