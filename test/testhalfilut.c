#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "lislib.h"


static LIS_INT
create_matrix(
	LIS_MATRIX *A)
{
	const LIS_INT n = 64;

	LIS_INT i;
	LIS_INT err;

	err =
		lis_matrix_create(
			LIS_COMM_WORLD,
			A);

	if( err ) return err;

	err =
		lis_matrix_set_size(
			*A,
			0,
			n);

	if( err ) return err;

	for(i=0;i<n;i++)
	{
		if( i>1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i-2,
					(LIS_SCALAR)-0.08,
					*A);

			if( err ) return err;
		}

		if( i>0 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i-1,
					(LIS_SCALAR)-0.71,
					*A);

			if( err ) return err;
		}

		err =
			lis_matrix_set_value(
				LIS_INS_VALUE,
				i,
				i,
				(LIS_SCALAR)
				(4.5 +
				 0.01*(double)i),
				*A);

		if( err ) return err;

		if( i<n-1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+1,
					(LIS_SCALAR)-0.53,
					*A);

			if( err ) return err;
		}

		if( i<n-2 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+2,
					(LIS_SCALAR)-0.06,
					*A);

			if( err ) return err;
		}
	}

	err =
		lis_matrix_set_type(
			*A,
			LIS_MATRIX_CSR);

	if( err ) return err;

	return lis_matrix_assemble(*A);
}


static int
normal_values_present(
	LIS_PRECON precon)
{
	LIS_INT i;

	for(i=0;i<precon->L->n;i++)
	{
		if( precon->L->nnz[i]>0 &&
			precon->L->value[i]==NULL )
		{
			return 0;
		}

		if( precon->U->nnz[i]>0 &&
			precon->U->value[i]==NULL )
		{
			return 0;
		}
	}

	return 1;
}


static int
normal_values_released(
	LIS_PRECON precon)
{
	LIS_INT i;

	for(i=0;i<precon->L->n;i++)
	{
		if( precon->L->nnz[i]>0 &&
			precon->L->value[i]!=NULL )
		{
			return 0;
		}

		if( precon->U->nnz[i]>0 &&
			precon->U->value[i]!=NULL )
		{
			return 0;
		}
	}

	return 1;
}


static double
relative_difference(
	LIS_VECTOR A,
	const LIS_SCALAR *ref,
	LIS_INT n)
{
	LIS_INT i;

	double nr;
	double nd;
	double d;

	nr = 0.0;
	nd = 0.0;

	for(i=0;i<n;i++)
	{
		nr =
			hypot(
				nr,
				(double)ref[i]);

		d =
			(double)A->value[i] -
			(double)ref[i];

		nd =
			hypot(
				nd,
				d);
	}

	if( nr==0.0 )
	{
		return nd;
	}

	return nd/nr;
}


static int
run_mode(
	LIS_MATRIX A,
	LIS_INT mode,
	int strict,
	LIS_SCALAR *ref,
	LIS_SCALAR *refh)
{
	LIS_SOLVER solver;
	LIS_PRECON precon;

	LIS_VECTOR B;
	LIS_VECTOR X;
	LIS_VECTOR XH;

	LIS_INT i;
	LIS_INT err;

	double rel;
	double relh;

	char options[256];

	int failed;

	solver = NULL;
	precon = NULL;

	B = NULL;
	X = NULL;
	XH = NULL;

	failed = 0;

	err =
		lis_solver_create(
			&solver);

	if( err )
	{
		return 1;
	}

	snprintf(
		options,
		sizeof(options),
		"-i gmres "
		"-p ilut "
		"-scale none "
		"-print none "
		"-precon_half %d",
		(int)mode);

	err =
		lis_solver_set_option(
			options,
			solver);

	if( err )
	{
		fprintf(
			stderr,
			"option parse mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_set_matrix(
			A,
			solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_precon_create(
			solver,
			&precon);

	if( err )
	{
		fprintf(
			stderr,
			"ILUT create mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		precon = NULL;
		failed = 1;
		goto cleanup;
	}

	solver->precon = precon;

	if( mode==LIS_PRECON_HALF_OFF ||
		strict )
	{
		if( precon->half!=NULL )
		{
			fprintf(
				stderr,
				"unexpected half mode=%d strict=%d\n",
				(int)mode,
				strict);

			failed = 1;
			goto cleanup;
		}

		if( !normal_values_present(precon) )
		{
			fprintf(
				stderr,
				"normal factors missing mode=%d strict=%d\n",
				(int)mode,
				strict);

			failed = 1;
			goto cleanup;
		}
	}
	else
	{
		if( precon->half==NULL )
		{
			fprintf(
				stderr,
				"half factors missing mode=%d\n",
				(int)mode);

			failed = 1;
			goto cleanup;
		}

#ifdef USE_QUAD_PRECISION

		if( !normal_values_present(precon) )
		{
			fprintf(
				stderr,
				"quad normal factors missing\n");

			failed = 1;
			goto cleanup;
		}

#else

		if( !normal_values_released(precon) )
		{
			fprintf(
				stderr,
				"normal factor values not released\n");

			failed = 1;
			goto cleanup;
		}

#endif
	}

	err = lis_vector_duplicate(A,&B);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&X);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&XH);
	if( err ) goto vector_fail;

	for(i=0;i<A->n;i++)
	{
		B->value[i] =
			(LIS_SCALAR)
			(0.2 +
			 0.013*(double)i);
	}

	err =
		lis_psolve_ilut_csr(
			solver,
			B,
			X);

	if( err )
	{
		fprintf(
			stderr,
			"psolve mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		failed = 1;
		goto cleanup;
	}

	err =
		lis_psolveh_ilut_csr(
			solver,
			B,
			XH);

	if( err )
	{
		fprintf(
			stderr,
			"psolveh mode=%d rc=%d\n",
			(int)mode,
			(int)err);

		failed = 1;
		goto cleanup;
	}

	if( mode==LIS_PRECON_HALF_OFF )
	{
		for(i=0;i<A->n;i++)
		{
			ref[i] =
				X->value[i];

			refh[i] =
				XH->value[i];
		}

		rel = 0.0;
		relh = 0.0;
	}
	else
	{
		rel =
			relative_difference(
				X,
				ref,
				A->n);

		relh =
			relative_difference(
				XH,
				refh,
				A->n);

		if( strict )
		{
			if( rel>1.0e-14 ||
				relh>1.0e-14 )
			{
				fprintf(
					stderr,
					"fallback mismatch "
					"rel=%.17g relh=%.17g\n",
					rel,
					relh);

				failed = 1;
				goto cleanup;
			}
		}
		else if( rel>1.0e-3 ||
			relh>1.0e-3 )
		{
			fprintf(
				stderr,
				"half mismatch mode=%d "
				"rel=%.17g relh=%.17g\n",
				(int)mode,
				rel,
				relh);

			failed = 1;
			goto cleanup;
		}
	}

	printf(
		"HALF_ILUT_STAGE3 PASSED "
		"mode=%d strict=%d "
		"rel=%.17g relh=%.17g "
		"half=%d",
		(int)mode,
		strict,
		rel,
		relh,
		precon->half!=NULL);

#ifdef USE_QUAD_PRECISION

	printf(
		" normal=retained quad=1\n");

#else

	if( precon->half!=NULL )
	{
		printf(
			" normal=released quad=0\n");
	}
	else
	{
		printf(
			" normal=retained quad=0\n");
	}

#endif

	goto cleanup;


vector_fail:

	fprintf(
		stderr,
		"vector allocation mode=%d rc=%d\n",
		(int)mode,
		(int)err);

	failed = 1;


cleanup:

	if( XH!=NULL )
		lis_vector_destroy(XH);

	if( X!=NULL )
		lis_vector_destroy(X);

	if( B!=NULL )
		lis_vector_destroy(B);

	if( precon!=NULL )
	{
		solver->precon = NULL;
		lis_precon_destroy(precon);
	}

	if( solver!=NULL )
	{
		lis_solver_destroy(solver);
	}

	return failed;
}



/* Exercise the opt-in symmetric-graph RCM and row-equilibrated ILUT path.
 * The old ILUT path remains unchanged when neither option is set.
 */
static int
check_ilut_ordering(LIS_MATRIX A)
{
    const char *variants[] = {
        "-ilut_rcm 1 -ilut_row_scale 0",
        "-ilut_rcm 0 -ilut_row_scale 1",
        "-ilut_rcm 1 -ilut_row_scale 1"
    };
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL,h=NULL,r=NULL,e=NULL;
    LIS_INT i,j,err,seen[64],nnz_l=0,nnz_u=0;
    LIS_REAL rn=0.0,bn=0.0,min_piv=0.0,max_piv=0.0;
    int v,failed=0;
    char opts[320];

    if(A->n!=64 || A->nprocs!=1) return 0;
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&h) ||
       lis_vector_duplicate(A,&r) ||
       lis_vector_duplicate(A,&e))
    { failed=1; goto cleanup; }
    lis_vector_set_all((LIS_SCALAR)1.0,e);
    err=lis_matvec(A,e,b);
    if(err) { failed=1; goto cleanup; }
    lis_vector_nrm2(b,&bn);

    for(v=0;v<3;v++)
    {
        solver=NULL;pc=NULL;
        err=lis_solver_create(&solver);
        if(err) { failed=1; goto cleanup; }
        snprintf(opts,sizeof(opts),
                 "-i gmres -p ilut -scale none -print none "
                 "-ilut_drop 0 -ilut_rate 10 -ilu_pivot_tol 1e-12 %s",
                 variants[v]);
        err=lis_solver_set_option(opts,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_solver_set_matrix(A,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_precon_create(solver,&pc);
        if(err) { failed=1; goto cleanup; }
        solver->precon=pc;
        err=lis_precon_ilut_factor_stats(pc,&nnz_l,&nnz_u,
                                         &min_piv,&max_piv);
        if(err || nnz_l<=0 || nnz_u<=0 ||
           !(min_piv>0.0 && max_piv>=min_piv))
        { failed=1; goto cleanup; }
        printf("ILUT_RCM_FACTOR_STATS variant=%d nnz_l=%d nnz_u=%d min_pivot=%.6g max_pivot=%.6g\n",
               v,(int)nnz_l,(int)nnz_u,(double)min_piv,(double)max_piv);
        if(!pc->ilut_perm || !pc->ilut_row_scale)
        { failed=1; goto cleanup; }

        memset(seen,0,sizeof(seen));
        for(i=0;i<64;i++)
        {
            j=pc->ilut_perm[i];
            if(j<0 || j>=64 || seen[j]++ ||
               !isfinite((double)pc->ilut_row_scale[i]) ||
               pc->ilut_row_scale[i]<=0.0)
            { failed=1; goto cleanup; }
        }

        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err=lis_psolve(solver,b,x);
        if(err) { failed=1; goto cleanup; }
        err=lis_matvec(A,x,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) { failed=1; goto cleanup; }
        lis_vector_nrm2(r,&rn);
        if(!(rn/bn<0.10))
        {
            fprintf(stderr,"ILUT RCM variant %d poor residual %.3g\n",
                    v,(double)(rn/bn));
            failed=1; goto cleanup;
        }

        /* Check transpose mapping independently from the normal solve. */
        lis_vector_set_all((LIS_SCALAR)0.0,h);
        err=lis_psolveh(solver,b,h);
        if(err) { failed=1; goto cleanup; }
        err=lis_matvech(A,h,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) { failed=1; goto cleanup; }
        lis_vector_nrm2(r,&rn);
        if(!(rn/bn<0.10))
        {
            fprintf(stderr,"ILUT RCM adjoint variant %d residual %.3g\n",
                    v,(double)(rn/bn));
            failed=1; goto cleanup;
        }
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }

    printf("ILUT_RCM_ROW_SCALE PASSED\n");
cleanup:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(e) lis_vector_destroy(e);
    if(r) lis_vector_destroy(r);
    if(h) lis_vector_destroy(h);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    return failed;
}


/*
 * Deterministic stabilized contact-KKT-like system:
 * [ K  B^T ][u] = [f]
 * [ B -0.1I][λ]   [g]
 *
 * K is an SPD spring chain, B enforces displacement differences.
 * The negative multiplier block makes the matrix indefinite.
 * This is an algebraic regression, NOT a full nonlinear contact FEA run.
 */
static int
check_ilut_contact_kkt(void)
{
    const LIS_INT nu=36,nlambda=12,n=48;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,u=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,q,iters=0,variant;
    LIS_REAL bn=0.0,rn=0.0;
    char opts[320];
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<nu;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
        if(i+1<nu)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
    }
    for(q=0;q<nlambda;q++)
    {
        LIS_INT l=nu+q,a=3*q,c=3*q+1;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,l,
                                 (LIS_SCALAR)-0.1,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,a,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,c,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,a,l,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,c,l,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&u) ||
       lis_vector_duplicate(A,&r)) goto done;
    lis_vector_set_all((LIS_SCALAR)1.0,u);
    err=lis_matvec(A,u,b);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || bn<=0.0) goto done;
    for(variant=0;variant<2;variant++)
    {
        solver=NULL;
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
            "-i fgmres -p ilut -scale none -print none "
            "-tol 1e-8 -maxiter 120 -restart 30 "
            "-ilut_drop 0.005 -ilut_rate 8 -ilu_pivot_tol 1e-9 %s",
            variant ? "-ilut_rcm 1 -ilut_row_scale 1" : "");
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err=lis_solve(A,b,x,solver);
        if(err) goto done;
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&rn);
        if(err) goto done;
        err=lis_solver_get_iter(solver,&iters);
        if(err || rn/bn>1e-7 || iters<1)
        {
            fprintf(stderr,"ILUT_CONTACT_KKT variant=%d iter=%d residual=%.3e\n",
                    variant,(int)iters,(double)(rn/bn));
            goto done;
        }
        printf("ILUT_CONTACT_KKT variant=%d iter=%d residual=%.6e\n",
               variant,(int)iters,(double)(rn/bn));
        lis_solver_destroy(solver); solver=NULL;
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_CONTACT_KKT_FAILED rc=%d\n",(int)err);
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(u) lis_vector_destroy(u);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    return failed;
}

/* Opt-in serial partitioned RAS ILUT with diagonal-approximate Schur
 * correction must remain separate from the legacy ILUT default.
 */
static int
check_ilut_ras_schur(LIS_MATRIX A)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL,e=NULL,r=NULL;
    LIS_INT err,overlap,schur,parallel,i;
    LIS_SCALAR reference[64];
    LIS_REAL bn=0.0,rn=0.0;
    int failed=0;
    char opts[300];

    if(A->nprocs!=1 || A->n!=64) return 0;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&e) || lis_vector_duplicate(A,&r))
    { failed=1; goto cleanup; }
    lis_vector_set_all((LIS_SCALAR)1.0,e);
    err=lis_matvec(A,e,b);
    if(err) { failed=1; goto cleanup; }
    err=lis_vector_nrm2(b,&bn);
    if(err) { failed=1; goto cleanup; }

    for(overlap=0;overlap<=2;overlap++)
    for(schur=0;schur<=1;schur++)
    for(parallel=0;parallel<=1;parallel++)
    {
        solver=NULL;pc=NULL;
        err=lis_solver_create(&solver);
        if(err) { failed=1; goto cleanup; }
        snprintf(opts,sizeof(opts),
            "-i gmres -p ilut -scale none -print none "
            "-iluc_drop 0 -iluc_rate 8 -ilu_pivot_tol 1e-12 "
            "-ilut_ras_domains 4 -ilut_ras_overlap %d -ilut_schur %d "
            "-ilut_ras_parallel %d",
            (int)overlap,(int)schur,(int)parallel);
        err=lis_solver_set_option(opts,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_solver_set_matrix(A,solver);
        if(err) { failed=1; goto cleanup; }
        err=lis_precon_create(solver,&pc);
        if(err) { failed=1; goto cleanup; }
        solver->precon=pc;
        if(pc->ilut_ras==NULL)
        { failed=1; goto cleanup; }
        err=lis_psolve(solver,b,x);
        if(err) { failed=1; goto cleanup; }
        if(parallel==0)
            for(i=0;i<64;i++) reference[i]=x->value[i];
        else
            for(i=0;i<64;i++)
                if(!isfinite((double)x->value[i]) ||
                   fabs((double)(x->value[i]-reference[i])) >
                   1e-10*(1.0+fabs((double)reference[i])))
                {
                    fprintf(stderr,
                        "ILUT RAS parallel mismatch: row=%ld overlap=%ld schur=%ld\n",
                        (long)i,(long)overlap,(long)schur);
                    failed=1; goto cleanup;
                }
        err=lis_matvec(A,x,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) { failed=1; goto cleanup; }
        err=lis_vector_nrm2(r,&rn);
        if(err || !(rn/bn<0.80))
        {
            fprintf(stderr,"ILUT RAS overlap=%d schur=%d residual %.5e\n",
                    (int)overlap,(int)schur,(double)(rn/bn));
            failed=1; goto cleanup;
        }
        /* No mathematically incorrect adjoint implementation. */
        err=lis_psolveh(solver,b,x);
        if(err!=LIS_ERR_NOT_IMPLEMENTED)
        {
            fprintf(stderr,"ILUT RAS unexpectedly accepted adjoint\n");
            failed=1; goto cleanup;
        }
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    printf("ILUT_RAS_SCHUR PASSED\n");

cleanup:
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(e) lis_vector_destroy(e);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    return failed;
}

/*
 * Deterministic stabilized contact-KKT-like system:
 * [ K  B^T ][u] = [f]
 * [ B -0.1I][λ]   [g]
 *
 * K is an SPD spring chain, B enforces displacement differences.
 * The negative multiplier block makes the matrix indefinite.
 * This is an algebraic regression, NOT a full nonlinear contact FEA run.
 */
static int
check_ilut_ras_contact_kkt(void)
{
    const LIS_INT nu=36,nlambda=12,n=48;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,u=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,i,q,iters=0,variant;
    LIS_REAL bn=0.0,rn=0.0;
    char opts[320];
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<nu;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,(LIS_SCALAR)4.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
        if(i+1<nu)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.8,A);
            if(err) goto done;
        }
    }
    for(q=0;q<nlambda;q++)
    {
        LIS_INT l=nu+q,a=3*q,c=3*q+1;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,l,
                                 (LIS_SCALAR)-0.1,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,a,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,l,c,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,a,l,(LIS_SCALAR)1.0,A);
        if(err) goto done;
        err=lis_matrix_set_value(LIS_INS_VALUE,c,l,(LIS_SCALAR)-1.0,A);
        if(err) goto done;
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    if(lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&u) ||
       lis_vector_duplicate(A,&r)) goto done;
    lis_vector_set_all((LIS_SCALAR)1.0,u);
    err=lis_matvec(A,u,b);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || bn<=0.0) goto done;
    for(variant=0;variant<2;variant++)
    {
        solver=NULL;
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
            "-i fgmres -p ilut -scale none -print none "
            "-tol 1e-8 -maxiter 120 -restart 30 "
            "-iluc_drop 0.005 -iluc_rate 8 -ilu_pivot_tol 1e-9 %s",
            variant ? "-ilut_ras_domains 4 -ilut_ras_overlap 1 -ilut_schur 1" : "");
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err=lis_solve(A,b,x,solver);
        if(err) goto done;
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&rn);
        if(err) goto done;
        err=lis_solver_get_iter(solver,&iters);
        if(err || rn/bn>1e-7 || iters<1)
        {
            fprintf(stderr,"ILUT_RAS_CONTACT_KKT variant=%d iter=%d residual=%.3e\n",
                    variant,(int)iters,(double)(rn/bn));
            goto done;
        }
        printf("ILUT_RAS_CONTACT_KKT variant=%d iter=%d residual=%.6e\n",
               variant,(int)iters,(double)(rn/bn));
        lis_solver_destroy(solver); solver=NULL;
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_RAS_CONTACT_KKT_FAILED rc=%d\n",(int)err);
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(u) lis_vector_destroy(u);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    return failed;
}


/*
 * Negative regression: the Schur approximation MUST reject a weak
 * interior diagonal instead of quietly treating diag(A_II)^-1 as zero.
 * Row 1 is interior for the 2-domain partition of this 8x8 system.
 */
static int
check_ilut_ras_unsafe_interior(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT i,err=0;
    int failed=1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,8);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<8;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
             (i==1 || i==2) ? (LIS_SCALAR)1.0e-20 : (LIS_SCALAR)3.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.4,A);
            if(err) goto done;
        }
        if(i<7)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.4,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i fgmres -p ilut -scale none -print none "
        "-iluc_drop 0 -iluc_rate 8 -ilu_pivot_tol 1e-9 "
        "-ilut_ras_domains 2 -ilut_ras_overlap 1 -ilut_schur 1",
        solver);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_precon_create(solver,&pc);
    if(err==LIS_ERR_ILL_ARG)
    {
        pc=NULL; /* failing constructor releases its partial preconditioner */
        printf("ILUT_RAS_WEAK_DIAG_REJECTED\n");
        failed=0;
        goto done;
    }
    fprintf(stderr,"Schur accepted weak interior pivot; rc=%ld\n",(long)err);
done:
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(A) lis_matrix_destroy(A);
    return failed;
}



/*
 * A single weak interior diagonal is promoted to Gamma. This tests only
 * finite, successful setup and application: promotion does not guarantee
 * convergence for indefinite KKT systems or large 3D interfaces.
 */
static int
check_ilut_ras_promoted_interior(void)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_INT i,err=0;
    int failed=1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,8);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<8;i++)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
               i==1 ? (LIS_SCALAR)1.0e-20 : (LIS_SCALAR)3.0,A);
        if(err) goto done;
        if(i>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                                     (LIS_SCALAR)-0.4,A);
            if(err) goto done;
        }
        if(i<7)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                                     (LIS_SCALAR)-0.4,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(
        "-i fgmres -p ilut -scale none -print none "
        "-iluc_drop 0 -iluc_rate 8 -ilu_pivot_tol 1e-9 "
        "-ilut_ras_domains 2 -ilut_ras_overlap 1 -ilut_schur 1",
        solver);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_precon_create(solver,&pc);
    if(err || !pc || !pc->ilut_ras) goto done;
    solver->precon=pc;
    if(lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x))
        goto done;
    lis_vector_set_all((LIS_SCALAR)1.0,b);
    err=lis_psolve(solver,b,x);
    if(err) goto done;
    for(i=0;i<8;i++)
        if(!isfinite((double)x->value[i])) goto done;
    printf("ILUT_RAS_WEAK_DIAG_PROMOTED\n");
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_RAS_WEAK_DIAG_PROMOTION_FAILED rc=%ld\n",
                (long)err);
    if(solver) solver->precon=NULL;
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    return failed;
}

/*
 * High-unit-contrast 48x48 contact-KKT algebraic regression:
 *       [ 1e8 K   B^T ] [u] = [f]
 *       [ B     -1e-8 I] [l]   [g]
 *
 * Compare unscaled, row equilibration, RCM+row equilibration, and
 * a manually block-balanced congruence D A D, where
 * D_u=1e-4 and D_l=1e4. The latter is a diagnostic reference,
 * NOT an implemented automatic two-sided equilibrator.
 *
 * Always check both the original true residual and D(Ax-b); the
 * latter gives multiplier equations meaningful weight despite
 * displacement stiffness being orders of magnitude larger.
 * Only the block-balanced reference is a hard convergence gate:
 * the other variants intentionally record numerical failures.
 * This is not a production contact FEA performance benchmark.
 */
static int
check_ilut_kkt_unit_contrast(void)
{
    const LIS_INT nu=36,nlambda=12,n=48;
    const double stiffness=1.0e8;
    const double du=1.0e-4,dl=1.0e4;
    const char *labels[]={"baseline","row","rcm_row","block_balanced"};
    const char *flags[]={"","-ilut_row_scale 1",
                          "-ilut_rcm 1 -ilut_row_scale 1",
                          "-ilut_rcm 1"};
    LIS_MATRIX A=NULL,Ab=NULL,M;
    LIS_VECTOR truth=NULL,b=NULL,bb=NULL,x=NULL,xb=NULL,r=NULL,wr=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,pass,i,q,variant,iters=0;
    LIS_REAL bn=0.0,bbn=0.0,rn=0.0,wrn=0.0;
    char opts[320];
    int failed=1,balanced_ok=0;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    err=lis_matrix_create(LIS_COMM_WORLD,&Ab);
    if(err) goto done;
    err=lis_matrix_set_size(Ab,0,n);
    if(err) goto done;

    for(pass=0;pass<2;pass++)
    {
        M=pass ? Ab : A;
        for(i=0;i<nu;i++)
        {
            double factor=pass ? du*du : 1.0;
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
                                     (LIS_SCALAR)(4.0*stiffness*factor),M);
            if(err) goto done;
            if(i>0)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i-1,
                         (LIS_SCALAR)(-0.8*stiffness*factor),M);
                if(err) goto done;
            }
            if(i+1<nu)
            {
                err=lis_matrix_set_value(LIS_INS_VALUE,i,i+1,
                         (LIS_SCALAR)(-0.8*stiffness*factor),M);
                if(err) goto done;
            }
        }
        for(q=0;q<nlambda;q++)
        {
            LIS_INT l=nu+q,a=3*q,c=3*q+1;
            double cross=pass ? du*dl : 1.0;
            double cc=pass ? dl*dl : 1.0;
            err=lis_matrix_set_value(LIS_INS_VALUE,l,l,
                                     (LIS_SCALAR)(-1e-8*cc),M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,l,a,
                                     (LIS_SCALAR)cross,M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,l,c,
                                     (LIS_SCALAR)(-cross),M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,a,l,
                                     (LIS_SCALAR)cross,M);
            if(err) goto done;
            err=lis_matrix_set_value(LIS_INS_VALUE,c,l,
                                     (LIS_SCALAR)(-cross),M);
            if(err) goto done;
        }
        err=lis_matrix_set_type(M,LIS_MATRIX_CSR);
        if(err) goto done;
        err=lis_matrix_assemble(M);
        if(err) goto done;
    }

    if(lis_vector_duplicate(A,&truth) ||
       lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&bb) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&xb) ||
       lis_vector_duplicate(A,&r) ||
       lis_vector_duplicate(A,&wr)) goto done;

    for(i=0;i<n;i++)
    {
        double di=i<nu ? du : dl;
        double yi=(i%2 ? -1.0 : 1.0)*(1.0+0.01*i);
        truth->value[i]=(LIS_SCALAR)(di*yi);
    }
    err=lis_matvec(A,truth,b);
    if(err) goto done;
    for(i=0;i<n;i++)
        bb->value[i]=b->value[i]*(i<nu ? du : dl);
    if(lis_vector_nrm2(b,&bn) || lis_vector_nrm2(bb,&bbn) ||
       !(bn>0.0) || !(bbn>0.0)) goto done;

    for(variant=0;variant<4;variant++)
    {
        LIS_MATRIX target=variant==3 ? Ab : A;
        LIS_VECTOR rhs=variant==3 ? bb : b;
        LIS_VECTOR solution=variant==3 ? xb : x;
        LIS_INT nnzl=0,nnzu=0;
        LIS_REAL minpivot=0.0,maxpivot=0.0;
        int solve_rc;
        solver=NULL;
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(opts,sizeof(opts),
                 "-i fgmres -p ilut -scale none -print none "
                 "-tol 1e-9 -maxiter 250 -restart 40 "
                 "-ilut_drop 0.001 -ilut_rate 12 "
                 "-ilu_pivot_tol 1e-12 %s",flags[variant]);
        err=lis_solver_set_option(opts,solver);
        if(err) goto done;
        lis_vector_set_all((LIS_SCALAR)0.0,solution);
        solve_rc=lis_solve(target,rhs,solution,solver);
        if(solve_rc)
        {
            printf("ILUT_KKT_CONTRAST variant=%s solve_rc=%d\n",
                   labels[variant],(int)solve_rc);
            lis_solver_destroy(solver); solver=NULL;
            continue;
        }
        if(variant==3)
            for(i=0;i<n;i++)
                x->value[i]=xb->value[i]*(i<nu ? du : dl);
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        for(i=0;i<n;i++)
            wr->value[i]=r->value[i]*(i<nu ? du : dl);
        if(lis_vector_nrm2(r,&rn) || lis_vector_nrm2(wr,&wrn))
            goto done;
        err=lis_solver_get_iter(solver,&iters);
        if(err) goto done;
        if(solver->precon)
            (void)lis_precon_ilut_factor_stats(solver->precon,
                    &nnzl,&nnzu,&minpivot,&maxpivot);
        printf("ILUT_KKT_CONTRAST variant=%s iter=%d "
               "true_residual=%.6e balanced_residual=%.6e "
               "nnz_lu=%d min_pivot=%.6e max_pivot=%.6e\n",
               labels[variant],(int)iters,(double)(rn/bn),
               (double)(wrn/bbn),(int)(nnzl+nnzu),
               (double)minpivot,(double)maxpivot);
        if(variant==3 && isfinite((double)rn) &&
           isfinite((double)wrn) && rn/bn<1e-7 && wrn/bbn<1e-7 &&
           iters>=1)
            balanced_ok=1;
        lis_solver_destroy(solver); solver=NULL;
    }
    failed=!balanced_ok;
done:
    if(failed)
        fprintf(stderr,"ILUT_KKT_CONTRAST_FAILED rc=%d balanced_ok=%d\n",
                (int)err,balanced_ok);
    if(solver) lis_solver_destroy(solver);
    if(wr) lis_vector_destroy(wr);
    if(r) lis_vector_destroy(r);
    if(xb) lis_vector_destroy(xb);
    if(x) lis_vector_destroy(x);
    if(bb) lis_vector_destroy(bb);
    if(b) lis_vector_destroy(b);
    if(truth) lis_vector_destroy(truth);
    if(Ab) lis_matrix_destroy(Ab);
    if(A) lis_matrix_destroy(A);
    return failed;
}



/*
 * Explicit synthetic 3D scalar seven-point stencil, 8^3 = 512 unknowns.
 * Four contiguous serial domains create a Gamma interface of 384 DOF,
 * forcing the sparse CSR-ILUT Schur path (dense cutoff is 256).
 *
 * This is NOT 3D vector elasticity, contact FEA, or an application
 * performance benchmark. It only checks finite, consistent application
 * on a larger-than-dense-limit synthetic interface.
 */
static int
check_ilut_sparse_schur_smoke(int scrambled,int half_storage,int row_scale,int local_rcm)
{
    const LIS_INT side=8,n=512;
    LIS_MATRIX A=NULL;
    LIS_VECTOR truth=NULL,b=NULL,x=NULL,ref=NULL,r=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT err=0,i,ix,iy,iz,parallel,rcm,row,ngamma=0,baseline_gamma=-1;
    LIS_REAL bn=0.0,rn=0.0;
    char options[320];
    int failed=1;
#ifdef _OPENMP
    int tsan_mode=getenv("LIS_RAS_TSAN")!=NULL;
    /* Isolate the OpenMP RAS *apply* from unrelated matrix/setup
     * OpenMP work in the existing LIS library for targeted TSan.
     * Production tests continue to use the normal thread settings. */
    if(tsan_mode) omp_set_num_threads(1);
#endif
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,n);
    if(err) goto done;
    if(A->nprocs!=1) { failed=0; goto done; }
    for(i=0;i<n;i++)
    {
        row=scrambled ? ((257*i+13)%n) : i;
        ix=i%side;
        iy=(i/side)%side;
        iz=i/(side*side);
        err=lis_matrix_set_value(LIS_INS_VALUE,row,row,
                                 (LIS_SCALAR)7.0,A);
        if(err) goto done;
        if(ix>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,row,(scrambled ? ((257*(i-1)+13)%n) : (i-1)),
                                     (LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
        if(ix+1<side)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,row,(scrambled ? ((257*(i+1)+13)%n) : (i+1)),
                                     (LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
        if(iy>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,row,(scrambled ? ((257*(i-side)+13)%n) : (i-side)),
                                     (LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
        if(iy+1<side)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,row,(scrambled ? ((257*(i+side)+13)%n) : (i+side)),
                                     (LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
        if(iz>0)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,row,(scrambled ? ((257*(i-side*side)+13)%n) : (i-side*side)),
                                     (LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
        if(iz+1<side)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,row,(scrambled ? ((257*(i+side*side)+13)%n) : (i+side*side)),
                                     (LIS_SCALAR)-1.0,A);
            if(err) goto done;
        }
    }
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_assemble(A);
    if(err) goto done;
    if(lis_vector_duplicate(A,&truth) ||
       lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) ||
       lis_vector_duplicate(A,&ref) ||
       lis_vector_duplicate(A,&r)) goto done;
    lis_vector_set_all((LIS_SCALAR)1.0,truth);
    err=lis_matvec(A,truth,b);
    if(err) goto done;
    if(lis_vector_nrm2(b,&bn) || !(bn>0.0)) goto done;
    for(rcm=0;rcm<=1;rcm++)
    for(parallel=0;parallel<=1;parallel++)
    {
        err=lis_solver_create(&solver);
        if(err) goto done;
        snprintf(options,sizeof(options),
                 "-i fgmres -p ilut -scale none -print none "
                 "-ilut_drop 0.001 -ilut_rate 8 -ilu_pivot_tol 1e-12 "
                 "-ilut_ras_domains 4 -ilut_ras_overlap 1 -ilut_schur 1 "
                 "-ilut_ras_parallel %d -ilut_ras_rcm %d "
                 "-ilut_row_scale %d -ilut_rcm %d",
                 (int)parallel,(int)rcm,row_scale,local_rcm);
        err=lis_solver_set_option(options,solver);
        if(err) goto done;
        if(half_storage)
        {
            err=lis_solver_set_precon_half(solver,LIS_PRECON_HALF_ON);
            if(err) goto done;
        }
        err=lis_solver_set_matrix(A,solver);
        if(err) goto done;
        err=lis_precon_create(solver,&pc);
        if(err || !pc || !pc->ilut_ras) goto done;
        err=lis_precon_ilut_ras_interface_size(pc,&ngamma);
        if(err || ngamma<=0 || ngamma>n) goto done;
        if(scrambled && parallel==0)
        {
            if(rcm==0) baseline_gamma=ngamma;
            else if(baseline_gamma<=0 || ngamma>=baseline_gamma)
            {
                fprintf(stderr,
                    "RCM did not reduce the interface on the scrambled fixture: raw=%ld rcm=%ld\n",
                    (long)baseline_gamma,(long)ngamma);
                goto done;
            }
        }
        solver->precon=pc;
#ifdef _OPENMP
        if(tsan_mode && parallel) omp_set_num_threads(4);
#endif
        err=lis_psolve(solver,b,x);
#ifdef _OPENMP
        if(tsan_mode) omp_set_num_threads(1);
#endif
        if(err) goto done;
        err=lis_matvec(A,x,r);
        if(err) goto done;
        err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
        if(err) goto done;
        err=lis_vector_nrm2(r,&rn);
        if(err || !isfinite((double)rn) ||
           !(rn/bn < 1.0)) goto done;
        for(i=0;i<n;i++)
        {
            if(!isfinite((double)x->value[i])) goto done;
            if(parallel==0) ref->value[i]=x->value[i];
            else if(fabs((double)(x->value[i]-ref->value[i])) >
                    1e-10*(1.0+fabs((double)ref->value[i]))) goto done;
        }
        printf("ILUT_RAS_SPARSE_SCHUR_CHECK n=%ld "
               "scrambled=%d rcm=%ld ngamma=%ld parallel=%ld true_residual=%.6e\n",
               (long)n,scrambled,(long)rcm,(long)ngamma,
               (long)parallel,(double)(rn/bn));
        solver->precon=NULL;
        lis_precon_destroy(pc);pc=NULL;
        lis_solver_destroy(solver);solver=NULL;
    }
    printf("ILUT_RAS_SPARSE_SCHUR_PASSED\n");
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_RAS_SPARSE_SCHUR_FAILED rc=%ld\n",(long)err);
    if(solver) solver->precon=NULL;
    if(pc) lis_precon_destroy(pc);
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(ref) lis_vector_destroy(ref);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(truth) lis_vector_destroy(truth);
    if(A) lis_matrix_destroy(A);
    return failed;
}


/*
 * Manually-invoked external Matrix Market FEA review harness.
 *
 * ./test/testhalfilut --fea-matrix A.mtx BASE_RATE RCM_RATE [DROP]
 *
 * Both modes solve the *original* system A x = A x_true, not the reordered
 * or equilibrated surrogate. Adjust the rates until factor nnz match to
 * within 5%; the output states explicitly if they do not.
 *
 * Intentionally NOT part of make check: no private FEA data are shipped.
 */
typedef struct {
    LIS_INT iterations,nnz_l,nnz_u;
    LIS_REAL residual,min_pivot,max_pivot;
    double setup_cpu_s,solve_cpu_s;
} LIS_ILUT_FEA_METRICS;

static int
run_ilut_fea_variant(LIS_MATRIX A,LIS_VECTOR b,LIS_INT reordered,
                     double rate,double drop,LIS_ILUT_FEA_METRICS *m)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR x=NULL,r=NULL;
    LIS_INT err=0;
    LIS_REAL bn=0.0,rn=0.0;
    clock_t t0,t1,t2;
    char opts[360];
    int failed=1;
    if(!m || !isfinite(rate) || !isfinite(drop) ||
       rate<=0.0 || drop<0.0 || drop>=1.0) return 1;
    memset(m,0,sizeof(*m));
    snprintf(opts,sizeof(opts),
             "-i fgmres -p ilut -scale none -print none -tol 1e-8 "
             "-maxiter 500 -restart 40 -ilu_pivot_tol 1e-10 "
             "-ilut_rate %.12g -ilut_drop %.12g %s",
             rate,drop,
             reordered ? "-ilut_rcm 1 -ilut_row_scale 1" : "");
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(opts,solver);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&r);
    if(err) goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);
    if(err) goto done;
    t0=clock();
    err=lis_precon_create(solver,&pc);
    t1=clock();
    if(err) { pc=NULL; goto done; }
    solver->precon=pc;
    err=lis_precon_ilut_factor_stats(pc,&m->nnz_l,&m->nnz_u,
                                     &m->min_pivot,&m->max_pivot);
    if(err) goto done;
    err=lis_solve_kernel(A,b,x,solver,pc);
    t2=clock();
    if(err) goto done;
    err=lis_solver_get_iter(solver,&m->iterations);
    if(err) goto done;
    err=lis_matvec(A,x,r);
    if(err) goto done;
    err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
    if(err) goto done;
    err=lis_vector_nrm2(r,&rn);
    if(err) goto done;
    err=lis_vector_nrm2(b,&bn);
    if(err || !(bn>0.0)) goto done;
    m->residual=rn/bn;
    m->setup_cpu_s=(double)(t1-t0)/CLOCKS_PER_SEC;
    m->solve_cpu_s=(double)(t2-t1)/CLOCKS_PER_SEC;
    if(!isfinite((double)m->residual) || m->residual>1e-7)
        goto done;
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_FEA_BENCH_FAILED mode=%s rc=%ld\n",
                reordered?"rcm":"baseline",(long)err);
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(r) lis_vector_destroy(r);
    if(x) lis_vector_destroy(x);
    return failed;
}

static int
run_ilut_fea_matrix(const char *file,double baseline_rate,
                    double rcm_rate,double drop)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR exact=NULL,b=NULL;
    LIS_ILUT_FEA_METRICS m[2];
    LIS_INT i,err=0;
    double fill_gap;
    int failed=1;
    if(!file || !isfinite(baseline_rate) || !isfinite(rcm_rate) ||
       baseline_rate<=0.0 || rcm_rate<=0.0 ||
       !isfinite(drop) || drop<0.0 || drop>=1.0)
        return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_input_matrix(A,(char *)file);
    if(err) goto done;
    if(A->nprocs!=1 || A->matrix_type!=LIS_MATRIX_CSR ||
       A->n!=A->gn || A->n<2)
    {
        err=LIS_ERR_NOT_IMPLEMENTED;
        goto done;
    }
    err=lis_vector_duplicate(A,&exact);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    /* Avoid an all-ones null mode in elasticity stiffness matrices. */
    for(i=0;i<A->n;i++)
        exact->value[i]=(LIS_SCALAR)(1.0+0.017*(double)((i*19)%29));
    err=lis_matvec(A,exact,b);
    if(err) goto done;
    if(run_ilut_fea_variant(A,b,0,baseline_rate,drop,&m[0]) ||
       run_ilut_fea_variant(A,b,1,rcm_rate,drop,&m[1]))
        goto done;
    fill_gap=fabs((double)(m[0].nnz_l+m[0].nnz_u) -
                  (double)(m[1].nnz_l+m[1].nnz_u))/
                  (double)(m[0].nnz_l+m[0].nnz_u);
    printf("ILUT_FEA_MATRIX file=%s dof=%ld nnz_A=%ld drop=%.9g "
           "fill_gap=%.5f fair_fill=%s\n",
           file,(long)A->n,(long)A->nnz,drop,fill_gap,
           fill_gap<=0.05 ? "yes" : "no");
    printf("variant,rate,nnz_L,nnz_U,nnz_LU,iterations,"
           "true_relres,min_pivot,max_pivot,setup_cpu_s,solve_cpu_s\n");
    for(i=0;i<2;i++)
        printf("%s,%.9g,%ld,%ld,%ld,%ld,%.9g,%.9g,%.9g,%.6f,%.6f\n",
               i==0?"baseline":"rcm_scaled",
               i==0?baseline_rate:rcm_rate,
               (long)m[i].nnz_l,(long)m[i].nnz_u,
               (long)(m[i].nnz_l+m[i].nnz_u),
               (long)m[i].iterations,(double)m[i].residual,
               (double)m[i].min_pivot,(double)m[i].max_pivot,
               m[i].setup_cpu_s,m[i].solve_cpu_s);
    failed=0;
done:
    if(failed)
        fprintf(stderr,"ILUT_FEA_MATRIX_FAILED file=%s rc=%ld\n",
                file?file:"(null)",(long)err);
    if(b) lis_vector_destroy(b);
    if(exact) lis_vector_destroy(exact);
    if(A) lis_matrix_destroy(A);
    return failed;
}

/* Optional, externally supplied Matrix Market comparison.
 *
 * ./test/testhalfilut --fea-ras A.mtx 8 8 8 0.005 8 1 1 1
 *
 * Rates are independent: base, RAS, RAS+Schur. The last four numbers
 * specify RAS domains, overlap, RCM-ranked partition, OpenMP local
 * apply. The benchmark measures wall time, NOT process CPU time.
 * Fill entries are a storage proxy, not peak RSS; tune rates until
 * the proxy fill_gap is <=5% before inferring cost/performance.
 *
 * No synthetic matrix is presented as a real contact-Newton snapshot.
 */
typedef struct {
    LIS_INT iters,entries,ngamma,err;
    LIS_REAL true_relres,probe_apply_ratio,probe_true_relres;
    LIS_INT worst_domain;
    LIS_REAL worst_local_ratio,min_local_pivot;
    double setup_wall,solve_wall;
    int converged;
} LIS_RAS_BENCH_RESULT;

static void
run_fea_ras_case(LIS_MATRIX A,LIS_VECTOR b,int mode,
                 double rate,double drop,LIS_INT domains,
                 LIS_INT overlap,LIS_INT rcm,LIS_INT parallel,LIS_INT row_scale,
                 LIS_INT local_rcm,LIS_RAS_BENCH_RESULT *out)
{
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR x=NULL,r=NULL,probe=NULL;
    LIS_INT nl=0,nu=0,err=0;
    LIS_REAL minpivot=0.0,maxpivot=0.0,bn=0.0,rn=0.0,pn=0.0;
    double t0=0.0,t1=0.0,t2=0.0;
    char opts[512];
    const char *stage="initial";
    const char *baseline_threads=getenv("LIS_FEA_BASELINE_THREADS");
    const char *baseline_rcm_env=getenv("LIS_FEA_BASELINE_RCM");
    LIS_INT baseline_rcm=(baseline_rcm_env &&
                         strcmp(baseline_rcm_env,"1")==0) ? 1 : 0;
#ifdef _OPENMP
    int saved_threads=omp_get_max_threads();
    int override_threads=0;
    if(mode==0 && baseline_threads && baseline_threads[0])
    {
        char *last=NULL;
        long num=strtol(baseline_threads,&last,10);
        if(last==baseline_threads || *last || num<1 || num>64)
        {
            memset(out,0,sizeof(*out));
            out->err=LIS_ERR_ILL_ARG;
            return;
        }
        omp_set_num_threads((int)num);
        override_threads=1;
    }
#else
    if(mode==0 && baseline_threads && strcmp(baseline_threads,"1")!=0)
    {
        memset(out,0,sizeof(*out));
        out->err=LIS_ERR_NOT_IMPLEMENTED;
        return;
    }
#endif
    memset(out,0,sizeof(*out));
    out->iters=-1;
    out->true_relres=-1.0;
    out->probe_apply_ratio=-1.0;
    out->probe_true_relres=-1.0;
    out->worst_domain=-1;
    out->worst_local_ratio=-1.0;
    out->min_local_pivot=-1.0;
    if(mode==0)
        snprintf(opts,sizeof(opts),
                 "-i fgmres -p ilut -scale none -print none "
                 "-tol 1e-8 -maxiter 500 -restart 40 "
                 "-ilu_pivot_tol 1e-10 -ilut_rate %.12g "
                 "-ilut_drop %.12g -ilut_row_scale %ld -ilut_rcm %ld",
                 rate,drop,(long)row_scale,(long)baseline_rcm);
    else
        snprintf(opts,sizeof(opts),
                 "-i fgmres -p ilut -scale none -print none "
                 "-tol 1e-8 -maxiter 500 -restart 40 "
                 "-ilu_pivot_tol 1e-10 -ilut_rate %.12g "
                 "-ilut_drop %.12g -ilut_ras_domains %ld "
                 "-ilut_ras_overlap %ld -ilut_ras_rcm %ld "
                 "-ilut_ras_parallel %ld -ilut_schur %d "
                 "-ilut_row_scale %ld -ilut_rcm %ld",
                 rate,drop,(long)domains,(long)overlap,
                 (long)rcm,(long)parallel,mode==2,(long)row_scale,
                 (long)local_rcm);
    err=lis_solver_create(&solver);
    if(err) goto done;
    err=lis_solver_set_option(opts,solver);
    if(err) goto done;
    err=lis_solver_set_matrix(A,solver);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x);
    if(err) goto done;
    err=lis_vector_duplicate(A,&r);
    if(err) goto done;
    err=lis_vector_duplicate(A,&probe);
    if(err) goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);
    if(err) goto done;
    t0=lis_wtime();
    stage="precon_create";
    err=lis_precon_create(solver,&pc);
    t1=lis_wtime();
    out->setup_wall=t1-t0;
    if(err) { pc=NULL; goto done; }
    solver->precon=pc;
    if(mode==0)
    {
        stage="factor_stats";
        err=lis_precon_ilut_factor_stats(pc,&nl,&nu,
                                          &minpivot,&maxpivot);
        if(err) goto done;
        if(nl<0 || nu<0) { err=LIS_ERR_ILL_ARG; goto done; }
        out->entries=nl+nu;
    }
    else
    {
        stage="ras_stats";
        err=lis_precon_ilut_ras_factor_entries(pc,&out->entries);
        if(err) goto done;
        err=lis_precon_ilut_ras_interface_size(pc,&out->ngamma);
        if(err) goto done;
        err=lis_precon_ilut_ras_probe(pc,b,&out->worst_domain,
                                    &out->worst_local_ratio,
                                    &out->min_local_pivot);
        if(err) goto done;
    }
    /*
     * Numerical ablation: inspect the uniterated preconditioner action
     * independently from FGMRES. This distinguishes an extreme local
     * factor response from a Krylov/solver-dispatch issue. Record the
     * norms, never infer instability from the final residual alone.
     */
    err=lis_vector_nrm2(b,&bn);
    if(err || !(bn>0.0)) { err=LIS_ERR_ILL_ARG; goto done; }
    stage="initial_psolve";
    err=lis_psolve(solver,b,probe);
    if(err) goto done;
    err=lis_vector_nrm2(probe,&pn);
    if(err || !isfinite((double)pn)) { err=LIS_ERR_ILL_ARG; goto done; }
    out->probe_apply_ratio=pn/bn;
    err=lis_matvec(A,probe,r);
    if(err) goto done;
    err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
    if(err) goto done;
    err=lis_vector_nrm2(r,&rn);
    if(err || !isfinite((double)rn)) { err=LIS_ERR_ILL_ARG; goto done; }
    out->probe_true_relres=rn/bn;
    t1=lis_wtime();
    stage="solve_kernel";
    err=lis_solve_kernel(A,b,x,solver,pc);
    t2=lis_wtime();
    out->solve_wall=t2-t1;
    if(err) goto done;
    err=lis_solver_get_iter(solver,&out->iters);
    if(err) goto done;
    err=lis_matvec(A,x,r);
    if(err) goto done;
    err=lis_vector_axpy((LIS_SCALAR)-1.0,b,r);
    if(err) goto done;
    err=lis_vector_nrm2(r,&rn);
    if(err) goto done;
    out->true_relres=rn/bn;
    out->converged=isfinite((double)out->true_relres) &&
                    out->true_relres<=1.0e-8;
done:
    out->err=err;
    if(err)
        fprintf(stderr,
            "FEA_RAS_STAGE_ERROR mode=%d stage=%s rc=%ld "
            "baseline_threads=%s baseline_rcm=%ld omp_max=%d\n",
            mode,stage,(long)err,
            baseline_threads ? baseline_threads : "(same as OMP)",
            (long)baseline_rcm,
#ifdef _OPENMP
            omp_get_max_threads()
#else
            1
#endif
            );
    if(pc)
    {
        if(solver) solver->precon=NULL;
        lis_precon_destroy(pc);
    }
    if(solver) lis_solver_destroy(solver);
    if(probe) lis_vector_destroy(probe);
    if(r) lis_vector_destroy(r);
    if(x) lis_vector_destroy(x);
#ifdef _OPENMP
    if(override_threads) omp_set_num_threads(saved_threads);
#endif
}

static int
run_fea_ras_compare(const char *file,double base_rate,double ras_rate,
                    double schur_rate,double drop,LIS_INT domains,
                    LIS_INT overlap,LIS_INT rcm,LIS_INT parallel,LIS_INT row_scale,
                    LIS_INT local_rcm)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR truth=NULL,b=NULL;
    LIS_RAS_BENCH_RESULT m[3];
    double rates[3],gap;
    LIS_INT i,err=0;
    int failed=1,mode;
    int mode_count=getenv("LIS_FEA_ONLY_BASELINE") ? 1 :
                   (getenv("LIS_FEA_SKIP_SCHUR") ? 2 : 3);
    const char *names[3]={"baseline","ras","ras_schur"};
    if(!file || !isfinite(base_rate) || !(base_rate>0.0) ||
       !isfinite(ras_rate) || !(ras_rate>0.0) ||
       !isfinite(schur_rate) || !(schur_rate>0.0) ||
       !isfinite(drop) || drop<0.0 || drop>=1.0 ||
       domains<2 || domains>64 || overlap<0 || overlap>2 ||
       (rcm!=0 && rcm!=1) || (parallel!=0 && parallel!=1) ||
       (row_scale!=0 && row_scale!=1) ||
       (local_rcm!=0 && local_rcm!=1))
        return 1;
    rates[0]=base_rate; rates[1]=ras_rate; rates[2]=schur_rate;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_type(A,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_input_matrix(A,(char *)file);
    if(err) goto done;
    if(A->nprocs!=1 || A->matrix_type!=LIS_MATRIX_CSR ||
       A->n!=A->gn || A->n<domains)
    {
        err=LIS_ERR_NOT_IMPLEMENTED;
        goto done;
    }
    err=lis_vector_duplicate(A,&truth);
    if(err) goto done;
    err=lis_vector_duplicate(A,&b);
    if(err) goto done;
    for(i=0;i<A->n;i++)
        truth->value[i]=(LIS_SCALAR)(1.0+0.017*(double)((i*19)%29));
    err=lis_matvec(A,truth,b);
    if(err) goto done;
    for(mode=0;mode<mode_count;mode++)
        run_fea_ras_case(A,b,mode,rates[mode],drop,domains,
                         overlap,rcm,parallel,row_scale,local_rcm,&m[mode]);
    printf("FEA_RAS_MATRIX file=%s n=%ld nnz_A=%ld "
           "domains=%ld overlap=%ld rcm=%ld parallel=%ld row_scale=%ld "
           "local_rcm=%ld OMP_NUM_THREADS=%s baseline_threads=%s "
           "baseline_rcm=%s\n",
           file,(long)A->n,(long)A->nnz,
           (long)domains,(long)overlap,(long)rcm,(long)parallel,
           (long)row_scale,(long)local_rcm,
           getenv("OMP_NUM_THREADS") ? getenv("OMP_NUM_THREADS") : "(unset)",
           getenv("LIS_FEA_BASELINE_THREADS") ?
             getenv("LIS_FEA_BASELINE_THREADS") : "(same as OMP)",
           getenv("LIS_FEA_BASELINE_RCM") ?
             getenv("LIS_FEA_BASELINE_RCM") : "0");
    printf("variant,rate,stored_entry_proxy,ngamma,iterations,"
           "true_relres,setup_wall_s,solve_wall_s,total_wall_s,"
           "fill_gap_vs_baseline,fair_fill,status,lis_error,"
           "probe_norm_ratio,probe_one_apply_true_relres,"
           "worst_domain,worst_local_inverse_ratio,min_local_pivot\n");
    for(mode=0;mode<mode_count;mode++)
    {
        gap=m[0].entries>0 && m[mode].entries>0 ?
            fabs((double)m[mode].entries-(double)m[0].entries)/
            (double)m[0].entries : -1.0;
        printf("%s,%.9g,%ld,%ld,%ld,%.9g,%.6f,%.6f,%.6f,"
               "%.5f,%s,%s,%ld,%.9g,%.9g,%ld,%.9g,%.9g\n",
               names[mode],rates[mode],(long)m[mode].entries,
               (long)m[mode].ngamma,(long)m[mode].iters,
               (double)m[mode].true_relres,
               m[mode].setup_wall,m[mode].solve_wall,
               m[mode].setup_wall+m[mode].solve_wall,
               gap,gap>=0.0 && gap<=0.05 ? "yes":"no",
               m[mode].err==0 && m[mode].converged ? "converged" :
               (m[mode].err ? "solver_error" : "not_converged"),
               (long)m[mode].err,
               (double)m[mode].probe_apply_ratio,
               (double)m[mode].probe_true_relres,
               (long)m[mode].worst_domain,
               (double)m[mode].worst_local_ratio,
               (double)m[mode].min_local_pivot);
    }
    failed=0;
    for(mode=0;mode<mode_count;mode++)
        if(m[mode].err || !m[mode].converged) failed=1;
done:
    if(err) fprintf(stderr,"FEA_RAS_INPUT_ERROR file=%s rc=%ld\n",
                    file,(long)err);
    if(b) lis_vector_destroy(b);
    if(truth) lis_vector_destroy(truth);
    if(A) lis_matrix_destroy(A);
    return failed;
}


int
main(
	int argc,
	char **argv)
{
	LIS_MATRIX A;

	LIS_SCALAR *ref;
	LIS_SCALAR *refh;

	LIS_INT err;

	int strict;
	int failed;

	A = NULL;
	ref = NULL;
	refh = NULL;

	failed = 0;

	strict =
		argc>1 &&
		strcmp(
			argv[1],
			"fallback")==0;

	err =
		lis_initialize(
			&argc,
			&argv);

	if( err )
	{
		return 1;
	}

	if(argc==2 && (strcmp(argv[1],"--sparse-schur-smoke")==0 ||
	               strcmp(argv[1],"--ras-rcm-smoke")==0 ||
	               strcmp(argv[1],"--ras-half-smoke")==0 ||
	               strcmp(argv[1],"--ras-row-scale-smoke")==0 ||
	               strcmp(argv[1],"--ras-local-rcm-smoke")==0))
	{
		failed=check_ilut_sparse_schur_smoke(
		    strcmp(argv[1],"--ras-rcm-smoke")==0,
		    strcmp(argv[1],"--ras-half-smoke")==0,
		    strcmp(argv[1],"--ras-half-smoke")==0 ||
		    strcmp(argv[1],"--ras-row-scale-smoke")==0,
		    strcmp(argv[1],"--ras-local-rcm-smoke")==0);
		lis_finalize();
		return failed ? 1 : 0;
	}

	if((argc==11 || argc==12 || argc==13) && strcmp(argv[1],"--fea-ras")==0)
	{
		char *end[10],*value[10];
		double br,rr,sr,drop;
		long domains,overlap,rcm,parallel,row_scale=0,local_rcm=0;
		int j,valid=1;
		for(j=0;j<argc-3;j++) { value[j]=argv[j+3]; end[j]=NULL; }
		br=strtod(value[0],&end[0]);
		rr=strtod(value[1],&end[1]);
		sr=strtod(value[2],&end[2]);
		drop=strtod(value[3],&end[3]);
		domains=strtol(value[4],&end[4],10);
		overlap=strtol(value[5],&end[5],10);
		rcm=strtol(value[6],&end[6],10);
		parallel=strtol(value[7],&end[7],10);
		if(argc>=12) row_scale=strtol(value[8],&end[8],10);
		if(argc==13) local_rcm=strtol(value[9],&end[9],10);
		for(j=0;j<argc-3;j++)
		    if(end[j]==value[j] || *end[j]!='\0') valid=0;
		failed=valid ? run_fea_ras_compare(argv[2],br,rr,sr,drop,
                    (LIS_INT)domains,(LIS_INT)overlap,
                    (LIS_INT)rcm,(LIS_INT)parallel,(LIS_INT)row_scale,
                    (LIS_INT)local_rcm) : 1;
		lis_finalize();
		return failed ? 1 : 0;
	}

	if(argc>=5 && argc<=6 &&
	   strcmp(argv[1],"--fea-matrix")==0)
	{
		char *end1=NULL,*end2=NULL,*end3=NULL;
		double baseline_rate=strtod(argv[3],&end1);
		double rcm_rate=strtod(argv[4],&end2);
		double drop=argc==6 ? strtod(argv[5],&end3) : 0.005;
		if(end1==argv[3] || *end1!='\0' ||
		   end2==argv[4] || *end2!='\0' ||
		   (argc==6 && (end3==argv[5] || *end3!='\0')))
			failed=1;
		else
			failed=run_ilut_fea_matrix(argv[2],baseline_rate,rcm_rate,drop);
		lis_finalize();
		return failed ? 1 : 0;
	}

	err = create_matrix(&A);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	ref =
		(LIS_SCALAR *)malloc(
			A->n*sizeof(LIS_SCALAR));

	refh =
		(LIS_SCALAR *)malloc(
			A->n*sizeof(LIS_SCALAR));

	if( ref==NULL ||
		refh==NULL )
	{
		failed = 1;
		goto cleanup;
	}

	failed |= check_ilut_ordering(A);
	failed |= check_ilut_contact_kkt();
	failed |= check_ilut_kkt_unit_contrast();
	failed |= check_ilut_ras_schur(A);
	failed |= check_ilut_ras_unsafe_interior();
	failed |= check_ilut_ras_promoted_interior();
	failed |= check_ilut_ras_contact_kkt();

	failed |=
		run_mode(
			A,
			LIS_PRECON_HALF_OFF,
			0,
			ref,
			refh);

	if( strict )
	{
		failed |=
			run_mode(
				A,
				LIS_PRECON_HALF_ON,
				1,
				ref,
				refh);
	}
	else
	{
		failed |=
			run_mode(
				A,
				LIS_PRECON_HALF_ON,
				0,
				ref,
				refh);

		failed |=
			run_mode(
				A,
				LIS_PRECON_HALF_NOPROBE,
				0,
				ref,
				refh);
	}

	if( !failed )
	{
		printf(
			"HALF_ILUT_STAGE3 COMPLETE PASS\n");
	}


cleanup:

	free(ref);
	free(refh);

	if( A!=NULL )
		lis_matrix_destroy(A);

	lis_finalize();

	return failed ? 1 : 0;
}
