/* Generated continuation; native LIS source copyright and license apply.
 * Sources: src/esolver/lis_esolver.c, src/esolver/lis_esolver_block.c, src/esolver/lis_esolver_feast.c, src/esolver/lis_krylov.c, src/solver/lis_solver.c, src/solver/lis_solver_bicgstab.c, src/solver/lis_solver_cg.c, src/solver/lis_solver_fgcrodr.c, src/solver/lis_solver_gcrodr.c, src/solver/lis_solver_gmres.c
 */
/* Generated from native LIS solver functions. Copyright notices remain in their source files; LIS BSD terms apply. */
struct linear_frame_initial { int pc,parent;
    LIS_SOLVER solver;
    LIS_PRECON M;
    LIS_VECTOR t;
    LIS_VECTOR r;
    LIS_REAL * bnrm2;
    LIS_Comm comm;
    LIS_INT output;
    LIS_INT conv;
    LIS_INT err;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR xx;
    LIS_VECTOR b;
    LIS_VECTOR p;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL tol_w;
    LIS_REAL tol_switch;
};
struct linear_frame_install { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * v;
    LIS_VECTOR * z;
    LIS_INT jdim;
    LIS_FGCRODR_STATE state;
    LIS_INT err;
    LIS_INT target;
    LIS_INT candidate;
    LIS_INT keep;
    LIS_INT j;
    LIS_INT pass;
    LIS_SCALAR alpha;
    LIS_REAL nrm;
};
struct linear_frame_harmonic_build { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * v;
    LIS_VECTOR * z;
    LIS_INT jdim;
    const LIS_SCALAR * hraw;
    LIS_INT hld;
    LIS_FGCRODR_STATE state;
    LIS_SCALAR * G;
    LIS_SCALAR * F;
    LIS_SCALAR * AP;
    LIS_SCALAR * BP;
    LIS_SCALAR * T;
    LIS_SCALAR * TQR;
    LIS_SCALAR * QW;
    LIS_SCALAR * RW;
    LIS_SCALAR * AINV;
    LIS_FGCRODR_EIG * eig;
    LIS_REAL * vr;
    LIS_REAL * vi;
    LIS_VECTOR * neww;
    LIS_VECTOR * newz;
    LIS_SCALAR ag;
    LIS_SCALAR bg;
    LIS_SCALAR s;
    LIS_INT target;
    LIS_INT psmall;
    LIS_INT qsmall;
    LIS_INT i;
    LIS_INT r;
    LIS_INT c;
    LIS_INT l;
    LIS_INT neig;
    LIS_INT nsel;
    LIS_INT eidx;
    LIS_INT qriter;
    LIS_REAL qrerr;
    LIS_INT err;
};
struct linear_frame_replace { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * neww;
    LIS_VECTOR * newz;
    LIS_INT nsel;
    LIS_INT count_update;
    LIS_INT action;
    LIS_REAL quality;
    LIS_REAL quality_min;
    LIS_FGCRODR_STATE oldstate;
    LIS_FGCRODR_STATE trial;
    LIS_INT saved_cross_action;
    LIS_REAL saved_quality;
    LIS_REAL saved_quality_min;
    LIS_INT err;
};
struct linear_frame_commit { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * neww;
    LIS_VECTOR * newz;
    LIS_INT nsel;
    LIS_INT err;
};
struct linear_frame_harmonic_update { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * v;
    LIS_VECTOR * z;
    LIS_INT jdim;
    const LIS_SCALAR * hraw;
    LIS_INT hld;
    const LIS_SCALAR * bcoef;
    LIS_INT bld;
    LIS_FGCRODR_STATE state;
    LIS_SCALAR * G;
    LIS_SCALAR * F;
    LIS_SCALAR * AP;
    LIS_SCALAR * BP;
    LIS_SCALAR * T;
    LIS_SCALAR * TQR;
    LIS_SCALAR * QW;
    LIS_SCALAR * RW;
    LIS_SCALAR * AINV;
    LIS_FGCRODR_EIG * eig;
    LIS_REAL * vr;
    LIS_REAL * vi;
    LIS_VECTOR * neww;
    LIS_VECTOR * newz;
    LIS_VECTOR * overlap_q;
    LIS_VECTOR * overlap_z;
    LIS_SCALAR ag;
    LIS_SCALAR bg;
    LIS_SCALAR s;
    LIS_INT oldk;
    LIS_INT target;
    LIS_INT psmall;
    LIS_INT qsmall;
    LIS_INT i;
    LIS_INT r;
    LIS_INT c;
    LIS_INT l;
    LIS_INT neig;
    LIS_INT nsel;
    LIS_INT eidx;
    LIS_INT qriter;
    LIS_REAL qrerr;
    LIS_INT err;
};
struct linear_frame_build { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * v;
    LIS_VECTOR * z;
    LIS_INT jdim;
    const LIS_SCALAR * hraw;
    LIS_INT hld;
    const LIS_SCALAR * bcoef;
    LIS_INT bld;
    LIS_FGCRODR_STATE state;
    LIS_INT err;
};
struct linear_frame_prepare { int pc,parent;
    LIS_SOLVER solver;
    LIS_FGCRODR_STATE state;
    LIS_VECTOR * scratch_c;
    LIS_SCALAR alpha;
    LIS_REAL nrm;
    LIS_REAL proj2;
    LIS_REAL qi;
    LIS_REAL qsum;
    LIS_REAL qmin;
    LIS_REAL qavg;
    LIS_INT oldk;
    LIS_INT i;
    LIS_INT j;
    LIS_INT good;
    LIS_INT need;
    LIS_INT action;
    LIS_INT err;
};
struct linear_frame_core { int pc,parent;
    LIS_SOLVER solver;
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR s;
    LIS_VECTOR * z;
    LIS_VECTOR * v;
    LIS_SCALAR * h;
    LIS_SCALAR * hraw;
    LIS_SCALAR * bcoef;
    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR rr;
    LIS_REAL aa_abs;
    LIS_REAL bb_abs;
    LIS_REAL givens_scale;
    LIS_REAL givens_a;
    LIS_REAL givens_b;
    LIS_REAL givens_norm;
    LIS_SCALAR givens_alpha;
    LIS_SCALAR givens_c;
    LIS_SCALAR givens_s;
    LIS_SCALAR t;
    LIS_SCALAR alpha;
    LIS_SCALAR gamma;
    LIS_REAL bnrm2;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_REAL rnorm;
    LIS_INT err;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT output;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    LIS_INT noimp_stop;
    LIS_INT candidate;
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT m;
    LIS_INT ii;
    LIS_INT i1;
    LIS_INT iiv;
    LIS_INT i1v;
    LIS_INT iih;
    LIS_INT jj;
    LIS_INT h_dim;
    LIS_INT cs;
    LIS_INT sn;
    LIS_FGCRODR_STATE state;
    LIS_INT use_recycle;
    LIS_INT kcur;
    LIS_INT ri;
    LIS_INT intra_recycle;
    LIS_INT intra_trace;
    LIS_INT recycle_counted;
    LIS_INT normalize_projected;
    LIS_REAL residual_before_projection;
    double time;
    double ptime;
    LIS_INT rows;
    size_t capacity;
};
struct linear_frame_fgcro { int pc,parent;
    LIS_SOLVER solver;
    LIS_INT err;
    LIS_INT solve_cross_action;
    LIS_FGCRODR_STATE state;
};
struct lis_rci_native_fgcro_restore { int active,method,status,strict; LIS_INT retval,native_code; uint64_t token,matvecs,precons; LIS_SOLVER solver; lis_rci_native_linear_request request;
    struct linear_frame_initial initial;
    struct linear_frame_install install;
    struct linear_frame_harmonic_build harmonic_build;
    struct linear_frame_replace replace;
    struct linear_frame_commit commit;
    struct linear_frame_harmonic_update harmonic_update;
    struct linear_frame_build build;
    struct linear_frame_prepare prepare;
    struct linear_frame_core core;
    struct linear_frame_fgcro fgcro;
};
#define L_STEP_AGAIN 1
#define L_STEP_YIELD 2
#define L_STEP_STOP 3
static int linear_suspend(lis_rci_native_fgcro_restore *,uint32_t,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);
static int linear_return(lis_rci_native_fgcro_restore *,int,LIS_INT);
#undef __FUNC__
#define __FUNC__ "lis_solver_get_initial_residual"
static int linear_step_initial(lis_rci_native_fgcro_restore *K) {
struct linear_frame_initial *F=&K->initial;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}








	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	F->A  = F->solver->A;
	F->b  = F->solver->b;
	F->x  = F->solver->x;
	F->xx = F->solver->x;
	F->output     = F->solver->options[LIS_OPTIONS_OUTPUT];
	F->conv       = F->solver->options[LIS_OPTIONS_CONV_COND];
	F->tol        = F->solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	F->tol_w      = F->solver->params[LIS_PARAMS_RESID_WEIGHT-LIS_OPTIONS_LEN];
	F->tol_switch = F->solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];



	if( F->M==NULL )
	{
		F->p = F->r;
	}
	else
	{
		F->p = F->t;
	}

	if( !F->solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
			F->pc=1;
			return linear_suspend(K,LIS_RCI_OP_A,F->A,F->x,F->p);
			resume_1: F->pc=0; F->err=K->retval;
			if( F->err )
			{
				F->solver->retcode = F->err;
				F->solver->iter = 0;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_FAILS);
			}
			lis_vector_xpay(F->b,-1,F->p);
	}
	else
	{
			lis_vector_copy(F->b,F->p);
	}

	switch(F->conv)
	{
	case LIS_CONV_COND_NRM2_R:
		lis_vector_nrm2(F->p,&F->nrm2);
		*F->bnrm2 = F->nrm2;
		F->solver->tol = F->tol;
		F->solver->tol_switch = F->tol_switch;
		break;
	case LIS_CONV_COND_NRM2_B:
		lis_vector_nrm2(F->p,&F->nrm2);
		lis_vector_nrm2(F->b,F->bnrm2);
		F->solver->tol = F->tol;
		F->solver->tol_switch = F->tol_switch;
		break;
	case LIS_CONV_COND_NRM1_B:
		lis_vector_nrm1(F->p,&F->nrm2);
		lis_vector_nrm1(F->b,F->bnrm2);
		F->solver->tol = *F->bnrm2*F->tol_w + F->tol;
		F->solver->tol_switch = *F->bnrm2*F->tol_w + F->tol_switch;
		break;
	}
	if( *F->bnrm2 == 0.0 )
	{
		*F->bnrm2 = 1.0;
	}
	else
	{
		*F->bnrm2 = 1.0 / *F->bnrm2;
	}
	F->solver->bnrm = *F->bnrm2;
	F->nrm2 = F->nrm2 * *F->bnrm2;

	if( F->output && (F->r->precision==LIS_PRECISION_QUAD && F->solver->precision!=LIS_PRECISION_SWITCH) )
	{
	  if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[0] = F->nrm2;
	  if( F->output & LIS_PRINT_OUT ) lis_printf(F->comm,"iteration: %5d  relative residual = %e\n", 0, (double)F->nrm2);
	}
	lis_solver_monitor_residual(
		F->solver,
		0,
		F->nrm2,
		fabs(F->solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN]),
		LIS_SOLVER_RESIDUAL_TRUE);

	if(
		F->nrm2 <= fabs(F->solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN])
		&&
		(
			F->nrm2==0.0
			||
			lis_solver_test_convergence(
				F->solver,
				0,
				F->nrm2,
				fabs(F->solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN]),
				LIS_SOLVER_RESIDUAL_TRUE,
				LIS_TRUE)
		)
	)
	{
		F->solver->retcode = LIS_SUCCESS;
		F->solver->iter    = 1;
		F->solver->resid   = F->nrm2;
		LIS_DEBUG_FUNC_OUT;
		return linear_return(K,F->parent,LIS_FAILS);
	}

	if( F->M!=NULL )
	{

		F->pc=2;
		return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->p,F->r);
		resume_2: F->pc=0; F->err=K->retval;
		if( F->err )
		{
			F->solver->retcode = F->err;
			F->solver->iter = 0;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_FAILS);
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_install_recycle"
static int linear_step_install(lis_rci_native_fgcro_restore *K) {
struct linear_frame_install *F=&K->install;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}
















    if(
        F->solver==NULL ||
        F->v==NULL ||
        F->z==NULL ||
        F->jdim<0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid FGCRO-DR recycle install arguments\n");

        return linear_return(K,F->parent,LIS_ERR_ILL_ARG);
    }


    F->state =
        F->solver->fgcrodr_state;


    if( F->state==NULL )
        return linear_return(K,F->parent,LIS_SUCCESS);








    if( !F->state->update )
        return linear_return(K,F->parent,LIS_SUCCESS);







    if(
        F->state->valid &&
        F->state->dim>0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    if(
        F->state->target_dim<=0 ||
        F->jdim<=0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }






    F->target =
        F->solver->options[LIS_OPTIONS_RESTART]-1;

    if( F->target>F->state->target_dim )
        F->target = F->state->target_dim;

    if( F->target>F->jdim )
        F->target = F->jdim;


    if( F->target<=0 )
        return linear_return(K,F->parent,LIS_SUCCESS);


    F->err =
        lis_fgcrodr_ensure_recycle_vectors(
            F->solver,
            F->target);

    if( F->err )
        return linear_return(K,F->parent,F->err);


    F->state =
        F->solver->fgcrodr_state;


    F->keep = 0;






    for(
        F->candidate=0;
        F->candidate<F->jdim && F->keep<F->target;
        F->candidate++ )
    {







        F->err =
            lis_vector_copy(
                F->v[F->candidate],
                F->state->w[F->keep]);

        if( F->err )
            goto fail;


        F->err =
            lis_vector_copy(
                F->z[F->candidate],
                F->state->z[F->keep]);

        if( F->err )
            goto fail;







        F->pc=1;
        return linear_suspend(K,LIS_RCI_OP_A,F->solver->A,F->state->z[F->keep],F->state->c[F->keep]);
        resume_1: F->pc=0; F->err=K->retval;
        if(F->err) goto fail;















        for(F->pass=0;F->pass<2;F->pass++)
        {
            for(F->j=0;F->j<F->keep;F->j++)
            {
                lis_vector_dot(
                    F->state->c[F->keep],
                    F->state->c[F->j],
                    &F->alpha);


                lis_vector_axpy(
                    -F->alpha,
                    F->state->c[F->j],
                    F->state->c[F->keep]);


                lis_vector_axpy(
                    -F->alpha,
                    F->state->z[F->j],
                    F->state->z[F->keep]);


                lis_vector_axpy(
                    -F->alpha,
                    F->state->w[F->j],
                    F->state->w[F->keep]);
            }
        }


        lis_vector_nrm2(
            F->state->c[F->keep],
            &F->nrm);






        if(
            F->nrm<=FGCRODR_BUILD_TINY ||
            F->nrm!=F->nrm )
        {

            continue;
        }


        lis_vector_scale(
            1.0/F->nrm,
            F->state->c[F->keep]);


        lis_vector_scale(
            1.0/F->nrm,
            F->state->z[F->keep]);


        lis_vector_scale(
            1.0/F->nrm,
            F->state->w[F->keep]);


        F->keep++;

    }




    F->state->dim =
        F->keep;

    F->state->valid =
        (F->keep>0)
        ? LIS_TRUE
        : LIS_FALSE;


    if( F->keep>0 )
    {
        F->state->matrix =
            F->solver->A;

        F->state->structure_generation =
            F->solver->A->structure_generation;

        F->state->value_generation =
            F->solver->A->value_generation;
    }


    if( F->keep>0 )
    {
        F->state->last_action =
            LIS_FGCRODR_RECYCLE_BUILD;

        F->state->builds++;







        F->state->quality =
            1.0;

        F->state->quality_min =
            1.0;
    }
    else
    {
        F->state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->rejects++;

        F->state->quality =
            0.0;

        F->state->quality_min =
            0.0;
    }


    return linear_return(K,F->parent,LIS_SUCCESS);


fail:

    lis_fgcrodr_vectors_destroy(
        F->state);

    F->state->last_action =
        LIS_FGCRODR_RECYCLE_REJECT;

    F->state->rejects++;

    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_try_harmonic_build"
static int linear_step_harmonic_build(lis_rci_native_fgcro_restore *K) {
struct linear_frame_harmonic_build *F=&K->harmonic_build;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}



    F->G = NULL;
    F->F = NULL;

    F->AP = NULL;
    F->BP = NULL;

    F->T = NULL;
    F->TQR = NULL;

    F->QW = NULL;
    F->RW = NULL;
    F->AINV = NULL;

    F->eig = NULL;

    F->vr = NULL;
    F->vi = NULL;

    F->neww = NULL;
    F->newz = NULL;
























    if(
        F->solver==NULL ||
        F->v==NULL ||
        F->z==NULL ||
        F->hraw==NULL ||
        F->jdim<=0 ||
        F->hld<F->jdim+1 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state==NULL ||
        !F->state->update ||
        F->state->valid ||
        F->state->target_dim<=0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->target =
        F->solver->options[
            LIS_OPTIONS_RESTART]
        - 1;

    if( F->target>F->state->target_dim )
        F->target = F->state->target_dim;

    if( F->target>F->jdim )
        F->target = F->jdim;


    if( F->target<=0 )
        return linear_return(K,F->parent,LIS_SUCCESS);


    F->psmall =
        F->jdim;

    F->qsmall =
        F->jdim+1;


    F->G =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->qsmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::G");

    F->F =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->qsmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::F");

    F->AP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::AP");

    F->BP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::BP");

    F->T =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::T");

    F->TQR =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::TQR");

    F->QW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::QW");

    F->RW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::RW");

    F->AINV =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_build::AINV");

    F->eig =
        (LIS_FGCRODR_EIG *)lis_malloc(
            (size_t)F->psmall
            *sizeof(*F->eig),
            "lis_fgcrodr_try_harmonic_build::eig");

    F->vr =
        (LIS_REAL *)lis_malloc(
            (size_t)F->psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_build::vr");

    F->vi =
        (LIS_REAL *)lis_malloc(
            (size_t)F->psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_build::vi");

    F->neww =
        (LIS_VECTOR *)lis_malloc(
            (size_t)F->target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_build::neww");

    F->newz =
        (LIS_VECTOR *)lis_malloc(
            (size_t)F->target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_build::newz");


    if( F->neww )
    {
        for(F->i=0;F->i<F->target;F->i++)
            F->neww[F->i] = NULL;
    }

    if( F->newz )
    {
        for(F->i=0;F->i<F->target;F->i++)
            F->newz[F->i] = NULL;
    }


    if(
        F->G==NULL ||
        F->F==NULL ||
        F->AP==NULL ||
        F->BP==NULL ||
        F->T==NULL ||
        F->TQR==NULL ||
        F->QW==NULL ||
        F->RW==NULL ||
        F->AINV==NULL ||
        F->eig==NULL ||
        F->vr==NULL ||
        F->vi==NULL ||
        F->neww==NULL ||
        F->newz==NULL )
    {
        F->err =
            LIS_ERR_OUT_OF_MEMORY;

        goto cleanup;
    }


    memset(
        F->G,
        0,
        (size_t)F->qsmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));

    memset(
        F->F,
        0,
        (size_t)F->qsmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));







    for(F->c=0;F->c<F->psmall;F->c++)
    {
        for(F->r=0;F->r<F->qsmall;F->r++)
        {
            F->G[
                F->r+F->c*F->qsmall] =
                F->hraw[
                    F->r+F->c*F->hld];
        }
    }






    F->err=lis_fgcrodr_flexible_overlap(
        F->v,F->qsmall,F->z,F->psmall,F->F,F->qsmall,0,0);
    if(F->err)
        goto cleanup;






    for(F->c=0;F->c<F->psmall;F->c++)
    {
        for(F->r=0;F->r<F->psmall;F->r++)
        {
            F->ag =
                0.0;

            F->bg =
                0.0;


            for(F->l=0;F->l<F->qsmall;F->l++)
            {
                F->ag +=
                    lis_fgcrodr_conj_scalar(
                        F->G[
                            F->l+F->r*F->qsmall])
                    *F->G[
                        F->l+F->c*F->qsmall];

                F->bg +=
                    lis_fgcrodr_conj_scalar(
                        F->G[
                            F->l+F->r*F->qsmall])
                    *F->F[
                        F->l+F->c*F->qsmall];
            }


            F->AP[
                F->r+F->c*F->psmall] =
                F->ag;

            F->BP[
                F->r+F->c*F->psmall] =
                F->bg;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->AP)
        ||
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->BP) )
    {
        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        F->AINV,
        F->AP,
        (size_t)F->psmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));


    F->err =
        lis_array_ge(
            F->psmall,
            F->AINV);

    if(
        F->err
        ||
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->AINV) )
    {




        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }





    for(F->c=0;F->c<F->psmall;F->c++)
    {
        for(F->r=0;F->r<F->psmall;F->r++)
        {
            F->s =
                0.0;


            for(F->l=0;F->l<F->psmall;F->l++)
            {
                F->s +=
                    F->AINV[
                        F->r+F->l*F->psmall]
                    *F->BP[
                        F->l+F->c*F->psmall];
            }


            F->T[
                F->r+F->c*F->psmall] =
                F->s;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->T) )
    {
        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        F->TQR,
        F->T,
        (size_t)F->psmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));









    F->qriter = 0;
    F->qrerr = 0.0;

    F->err =
        lis_eigen_general_complex_solve(
            F->psmall,
            F->TQR,
            F->RW,
            F->QW,
            &F->qriter,
            &F->qrerr);

    if( F->err )
    {
        if( F->err==LIS_ERR_OUT_OF_MEMORY )
            goto cleanup;

        F->err = LIS_SUCCESS;
        goto cleanup;
    }


    if(
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->TQR)
        ||
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->QW) )
    {
        F->err = LIS_SUCCESS;
        goto cleanup;
    }


    lis_fgcrodr_sort_complex_eigenpairs(
        F->psmall,
        F->RW,
        F->QW);

    F->neig = F->psmall;



    for(F->i=0;F->i<F->target;F->i++)
    {
        F->err =
            lis_vector_duplicate(
                F->solver->A,
                &F->neww[F->i]);

        if( F->err )
            goto cleanup;


        F->err =
            lis_vector_duplicate(
                F->solver->A,
                &F->newz[F->i]);

        if( F->err )
            goto cleanup;
    }


    F->nsel = 0;


    for(
        F->eidx=0;
        F->eidx<F->neig && F->nsel<F->target;
        F->eidx++ )
    {
        F->err =
            lis_fgcrodr_build_harmonic_vector(
                F->v,
                F->z,
                F->jdim,
                &F->QW[F->eidx*F->psmall],
                F->neww[F->nsel],
                F->newz[F->nsel]);

        if( F->err )
        {
            if( F->err==LIS_ERR_OUT_OF_MEMORY )
                goto cleanup;

            F->err = LIS_SUCCESS;
            continue;
        }


        F->nsel++;
    }










    if( F->nsel!=F->target )
    {
        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    K->install.solver = F->solver;
    K->install.v = F->neww;
    K->install.z = F->newz;
    K->install.jdim = F->nsel;
    K->install.pc=0; K->install.parent=3;
    F->pc=1; K->active=2; return L_STEP_AGAIN;
    resume_1: F->pc=0; F->err=K->retval;

    if( F->err )
        goto cleanup;


    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state!=NULL &&
        F->state->valid &&
        F->state->dim>0 )
    {
    }


cleanup:

    if( F->neww )
    {
        for(F->i=0;F->i<F->target;F->i++)
        {
            if( F->neww[F->i] )
                lis_vector_destroy(
                    F->neww[F->i]);
        }
    }


    if( F->newz )
    {
        for(F->i=0;F->i<F->target;F->i++)
        {
            if( F->newz[F->i] )
                lis_vector_destroy(
                    F->newz[F->i]);
        }
    }


    lis_free(F->neww);
    lis_free(F->newz);

    lis_free(F->G);
    lis_free(F->F);

    lis_free(F->AP);
    lis_free(F->BP);

    lis_free(F->T);
    lis_free(F->TQR);

    lis_free(F->QW);
    lis_free(F->RW);
    lis_free(F->AINV);

    lis_free(F->eig);

    lis_free(F->vr);
    lis_free(F->vi);


    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_replace_recycle_basis"
static int linear_step_replace(lis_rci_native_fgcro_restore *K) {
struct linear_frame_replace *F=&K->replace;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}












    F->oldstate =
        F->solver->fgcrodr_state;


    F->saved_cross_action =
        F->oldstate
        ? F->oldstate->cross_action
        : LIS_FGCRODR_RECYCLE_NONE;

    F->saved_quality =
        F->oldstate
        ? F->oldstate->quality
        : 0.0;

    F->saved_quality_min =
        F->oldstate
        ? F->oldstate->quality_min
        : 0.0;


    if(
        F->oldstate==NULL ||
        F->nsel<=0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->trial =
        (LIS_FGCRODR_STATE)lis_malloc(
            sizeof(
                struct LIS_FGCRODR_STATE_STRUCT),
            "lis_fgcrodr_replace_recycle_basis::trial");

    if( F->trial==NULL )
    {
        LIS_SETERR_MEM(
            sizeof(
                struct LIS_FGCRODR_STATE_STRUCT));

        return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
    }


    lis_fgcrodr_state_defaults(
        F->trial);





    F->trial->target_dim =
        F->oldstate->target_dim;

    F->trial->keep_threshold =
        F->oldstate->keep_threshold;

    F->trial->reject_threshold =
        F->oldstate->reject_threshold;

    F->trial->update =
        F->oldstate->update;





    F->solver->fgcrodr_state =
        F->trial;


    K->install.solver = F->solver;
    K->install.v = F->neww;
    K->install.z = F->newz;
    K->install.jdim = F->nsel;
    K->install.pc=0; K->install.parent=4;
    F->pc=1; K->active=2; return L_STEP_AGAIN;
    resume_1: F->pc=0; F->err=K->retval;


    F->solver->fgcrodr_state =
        F->oldstate;


    if( F->err )
    {
        lis_fgcrodr_vectors_destroy(
            F->trial);

        lis_free(
            F->trial);

        return linear_return(K,F->parent,F->err);
    }






    if(
        !F->trial->valid ||
        F->trial->dim<=0 )
    {
        lis_fgcrodr_vectors_destroy(
            F->trial);

        lis_free(
            F->trial);

        return linear_return(K,F->parent,LIS_SUCCESS);
    }








    lis_fgcrodr_vectors_destroy(
        F->oldstate);


    F->oldstate->w =
        F->trial->w;

    F->oldstate->z =
        F->trial->z;

    F->oldstate->c =
        F->trial->c;

    F->oldstate->dim =
        F->trial->dim;

    F->oldstate->alloc =
        F->trial->alloc;

    F->oldstate->valid =
        F->trial->valid;


    F->oldstate->matrix =
        F->trial->matrix;

    F->oldstate->structure_generation =
        F->trial->structure_generation;

    F->oldstate->value_generation =
        F->trial->value_generation;







    if(
        F->saved_cross_action
        !=LIS_FGCRODR_RECYCLE_NONE )
    {
        F->oldstate->quality =
            F->saved_quality;

        F->oldstate->quality_min =
            F->saved_quality_min;

        F->oldstate->last_action =
            F->saved_cross_action;

        F->oldstate->cross_action =
            F->saved_cross_action;
    }
    else
    {
        F->oldstate->quality =
            F->quality;

        F->oldstate->quality_min =
            F->quality_min;

        F->oldstate->last_action =
            F->action;

        F->oldstate->cross_action =
            LIS_FGCRODR_RECYCLE_NONE;
    }


    if( F->count_update )
        F->oldstate->updates++;





    F->trial->w = NULL;
    F->trial->z = NULL;
    F->trial->c = NULL;

    F->trial->dim = 0;
    F->trial->alloc = 0;
    F->trial->valid = LIS_FALSE;


    lis_free(
        F->trial);


    return linear_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_commit_recycle_update"
static int linear_step_commit(lis_rci_native_fgcro_restore *K) {
struct linear_frame_commit *F=&K->commit;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}



    K->replace.solver = F->solver;
    K->replace.neww = F->neww;
    K->replace.newz = F->newz;
    K->replace.nsel = F->nsel;
    K->replace.count_update = LIS_TRUE;
    K->replace.action = LIS_FGCRODR_RECYCLE_KEEP;
    K->replace.quality = 1.0;
    K->replace.quality_min = 1.0;
    K->replace.pc=0; K->replace.parent=5;
    F->pc=1; K->active=4; return L_STEP_AGAIN;
    resume_1: F->pc=0; F->err=K->retval;
    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_try_harmonic_update"
static int linear_step_harmonic_update(lis_rci_native_fgcro_restore *K) {
struct linear_frame_harmonic_update *F=&K->harmonic_update;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}



    F->G = NULL;
    F->F = NULL;

    F->AP = NULL;
    F->BP = NULL;

    F->T = NULL;
    F->TQR = NULL;

    F->QW = NULL;
    F->RW = NULL;
    F->AINV = NULL;

    F->eig = NULL;

    F->vr = NULL;
    F->vi = NULL;

    F->neww = NULL;
    F->newz = NULL;
    F->overlap_q = NULL;
    F->overlap_z = NULL;


























    if(
        F->solver==NULL ||
        F->v==NULL ||
        F->z==NULL ||
        F->hraw==NULL ||
        F->bcoef==NULL ||
        F->jdim<=0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state==NULL ||
        !F->state->update ||
        !F->state->valid ||
        F->state->dim<=0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->oldk =
        F->state->dim;


    if(
        F->bld<F->oldk ||
        F->hld<F->jdim+1 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->target =
        F->solver->options[
            LIS_OPTIONS_RESTART]
        - 1;

    if( F->target>F->state->target_dim )
        F->target = F->state->target_dim;





    if( F->target<=0 )
        return linear_return(K,F->parent,LIS_SUCCESS);


    F->psmall =
        F->oldk+F->jdim;

    F->qsmall =
        F->oldk+F->jdim+1;


    if( F->target>F->psmall )
        F->target = F->psmall;


    F->G =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->qsmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::G");

    F->F =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->qsmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::F");

    F->AP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::AP");

    F->BP =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::BP");

    F->T =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::T");

    F->TQR =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::TQR");

    F->QW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::QW");

    F->RW =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::RW");

    F->AINV =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->psmall
            *(size_t)F->psmall
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_try_harmonic_update::AINV");

    F->eig =
        (LIS_FGCRODR_EIG *)lis_malloc(
            (size_t)F->psmall
            *sizeof(*F->eig),
            "lis_fgcrodr_try_harmonic_update::eig");

    F->vr =
        (LIS_REAL *)lis_malloc(
            (size_t)F->psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_update::vr");

    F->vi =
        (LIS_REAL *)lis_malloc(
            (size_t)F->psmall
            *sizeof(LIS_REAL),
            "lis_fgcrodr_try_harmonic_update::vi");

    F->neww =
        (LIS_VECTOR *)lis_malloc(
            (size_t)F->target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_update::neww");

    F->newz =
        (LIS_VECTOR *)lis_malloc(
            (size_t)F->target
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_try_harmonic_update::newz");


    if( F->neww )
    {
        for(F->i=0;F->i<F->target;F->i++)
            F->neww[F->i] = NULL;
    }

    if( F->newz )
    {
        for(F->i=0;F->i<F->target;F->i++)
            F->newz[F->i] = NULL;
    }


    if(
        F->G==NULL ||
        F->F==NULL ||
        F->AP==NULL ||
        F->BP==NULL ||
        F->T==NULL ||
        F->TQR==NULL ||
        F->QW==NULL ||
        F->RW==NULL ||
        F->AINV==NULL ||
        F->eig==NULL ||
        F->vr==NULL ||
        F->vi==NULL ||
        F->neww==NULL ||
        F->newz==NULL )
    {
        F->err =
            LIS_ERR_OUT_OF_MEMORY;

        goto cleanup;
    }


    memset(
        F->G,
        0,
        (size_t)F->qsmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));

    memset(
        F->F,
        0,
        (size_t)F->qsmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));







    for(F->i=0;F->i<F->oldk;F->i++)
    {
        F->G[
            F->i+F->i*F->qsmall] =
            1.0;
    }









    for(F->c=0;F->c<F->jdim;F->c++)
    {
        for(F->r=0;F->r<F->oldk;F->r++)
        {
            F->G[
                F->r
                +(F->oldk+F->c)*F->qsmall] =
                F->bcoef[
                    F->r+F->c*F->bld];
        }


        for(F->r=0;F->r<=F->jdim;F->r++)
        {
            F->G[
                F->oldk+F->r
                +(F->oldk+F->c)*F->qsmall] =
                F->hraw[
                    F->r+F->c*F->hld];
        }
    }






    F->overlap_q=(LIS_VECTOR *)lis_malloc(
        (size_t)F->qsmall*sizeof(LIS_VECTOR),
        "lis_fgcrodr_try_harmonic_update::overlap_q");
    F->overlap_z=(LIS_VECTOR *)lis_malloc(
        (size_t)F->psmall*sizeof(LIS_VECTOR),
        "lis_fgcrodr_try_harmonic_update::overlap_z");
    if(F->overlap_q==NULL || F->overlap_z==NULL)
    {
        F->err=LIS_ERR_OUT_OF_MEMORY;
        goto cleanup;
    }
    for(F->i=0;F->i<F->oldk;F->i++)
    {
        F->overlap_q[F->i]=F->state->c[F->i];
        F->overlap_z[F->i]=F->state->z[F->i];
    }
    for(F->i=0;F->i<F->jdim;F->i++)
    {
        F->overlap_q[F->oldk+F->i]=F->v[F->i];
        F->overlap_z[F->oldk+F->i]=F->z[F->i];
    }
    F->overlap_q[F->oldk+F->jdim]=F->v[F->jdim];

    F->err=lis_fgcrodr_flexible_overlap(
        F->overlap_q,F->qsmall,F->overlap_z,F->psmall,F->F,F->qsmall,0,0);
    if(F->err)
        goto cleanup;






    for(F->c=0;F->c<F->psmall;F->c++)
    {
        for(F->r=0;F->r<F->psmall;F->r++)
        {
            F->ag = 0.0;
            F->bg = 0.0;


            for(F->l=0;F->l<F->qsmall;F->l++)
            {
                F->ag +=
                    lis_fgcrodr_conj_scalar(
                        F->G[
                            F->l+F->r*F->qsmall])
                    *F->G[
                        F->l+F->c*F->qsmall];

                F->bg +=
                    lis_fgcrodr_conj_scalar(
                        F->G[
                            F->l+F->r*F->qsmall])
                    *F->F[
                        F->l+F->c*F->qsmall];
            }


            F->AP[
                F->r+F->c*F->psmall] =
                F->ag;

            F->BP[
                F->r+F->c*F->psmall] =
                F->bg;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->AP)
        ||
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->BP) )
    {
        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        F->AINV,
        F->AP,
        (size_t)F->psmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));


    F->err =
        lis_array_ge(
            F->psmall,
            F->AINV);

    if(
        F->err
        ||
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->AINV) )
    {



        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }





    for(F->c=0;F->c<F->psmall;F->c++)
    {
        for(F->r=0;F->r<F->psmall;F->r++)
        {
            F->s = 0.0;


            for(F->l=0;F->l<F->psmall;F->l++)
            {
                F->s +=
                    F->AINV[
                        F->r+F->l*F->psmall]
                    *F->BP[
                        F->l+F->c*F->psmall];
            }


            F->T[
                F->r+F->c*F->psmall] =
                F->s;
        }
    }


    if(
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->T) )
    {
        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    memcpy(
        F->TQR,
        F->T,
        (size_t)F->psmall
        *(size_t)F->psmall
        *sizeof(LIS_SCALAR));









    F->qriter = 0;
    F->qrerr = 0.0;

    F->err =
        lis_eigen_general_complex_solve(
            F->psmall,
            F->TQR,
            F->RW,
            F->QW,
            &F->qriter,
            &F->qrerr);

    if( F->err )
    {
        if( F->err==LIS_ERR_OUT_OF_MEMORY )
            goto cleanup;

        F->err = LIS_SUCCESS;
        goto cleanup;
    }


    if(
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->TQR)
        ||
        !lis_fgcrodr_dense_finite(
            F->psmall,
            F->QW) )
    {
        F->err = LIS_SUCCESS;
        goto cleanup;
    }


    lis_fgcrodr_sort_complex_eigenpairs(
        F->psmall,
        F->RW,
        F->QW);

    F->neig = F->psmall;



    for(F->i=0;F->i<F->target;F->i++)
    {
        F->err =
            lis_vector_duplicate(
                F->solver->A,
                &F->neww[F->i]);

        if( F->err )
            goto cleanup;


        F->err =
            lis_vector_duplicate(
                F->solver->A,
                &F->newz[F->i]);

        if( F->err )
            goto cleanup;
    }


    F->nsel = 0;


    for(
        F->eidx=0;
        F->eidx<F->neig && F->nsel<F->target;
        F->eidx++ )
    {
        F->err =
            lis_fgcrodr_build_augmented_vector(
                F->state,
                F->v,
                F->z,
                F->jdim,
                F->oldk,
                &F->QW[F->eidx*F->psmall],
                F->neww[F->nsel],
                F->newz[F->nsel]);

        if( F->err )
        {
            if( F->err==LIS_ERR_OUT_OF_MEMORY )
                goto cleanup;

            F->err = LIS_SUCCESS;
            continue;
        }


        F->nsel++;
    }







    if( F->nsel!=F->target )
    {
        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    K->commit.solver = F->solver;
    K->commit.neww = F->neww;
    K->commit.newz = F->newz;
    K->commit.nsel = F->nsel;
    K->commit.pc=0; K->commit.parent=6;
    F->pc=1; K->active=5; return L_STEP_AGAIN;
    resume_1: F->pc=0; F->err=K->retval;

    if( F->err )
        goto cleanup;


    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state!=NULL &&
        F->state->valid &&
        F->state->dim>0 )
    {
    }


cleanup:

    if( F->neww )
    {
        for(F->i=0;F->i<F->target;F->i++)
        {
            if( F->neww[F->i] )
                lis_vector_destroy(
                    F->neww[F->i]);
        }
    }


    if( F->newz )
    {
        for(F->i=0;F->i<F->target;F->i++)
        {
            if( F->newz[F->i] )
                lis_vector_destroy(
                    F->newz[F->i]);
        }
    }


    lis_free(F->neww);
    lis_free(F->newz);
    lis_free(F->overlap_q);
    lis_free(F->overlap_z);

    lis_free(F->G);
    lis_free(F->F);

    lis_free(F->AP);
    lis_free(F->BP);

    lis_free(F->T);
    lis_free(F->TQR);

    lis_free(F->QW);
    lis_free(F->RW);
    lis_free(F->AINV);

    lis_free(F->eig);

    lis_free(F->vr);
    lis_free(F->vi);


    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_build_recycle"
static int linear_step_build(lis_rci_native_fgcro_restore *K) {
struct linear_frame_build *F=&K->build;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}






    if( F->solver==NULL )
        return linear_return(K,F->parent,LIS_ERR_ILL_ARG);


    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state==NULL ||
        !F->state->update )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }







    if(
        F->state->valid &&
        F->state->dim>0 )
    {

        K->harmonic_update.solver = F->solver;
        K->harmonic_update.v = F->v;
        K->harmonic_update.z = F->z;
        K->harmonic_update.jdim = F->jdim;
        K->harmonic_update.hraw = F->hraw;
        K->harmonic_update.hld = F->hld;
        K->harmonic_update.bcoef = F->bcoef;
        K->harmonic_update.bld = F->bld;
        K->harmonic_update.pc=0; K->harmonic_update.parent=7;
        F->pc=1; K->active=6; return L_STEP_AGAIN;
        resume_1: F->pc=0; F->err=K->retval;
    return linear_return(K,F->parent,F->err);

    }








    K->harmonic_build.solver = F->solver;
    K->harmonic_build.v = F->v;
    K->harmonic_build.z = F->z;
    K->harmonic_build.jdim = F->jdim;
    K->harmonic_build.hraw = F->hraw;
    K->harmonic_build.hld = F->hld;
    K->harmonic_build.pc=0; K->harmonic_build.parent=7;
    F->pc=2; K->active=3; return L_STEP_AGAIN;
    resume_2: F->pc=0; F->err=K->retval;


    if( F->err )
        return linear_return(K,F->parent,F->err);


    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state!=NULL &&
        F->state->valid &&
        F->state->dim>0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }





    K->install.solver = F->solver;
    K->install.v = F->v;
    K->install.z = F->z;
    K->install.jdim = F->jdim;
    K->install.pc=0; K->install.parent=7;
    F->pc=3; K->active=2; return L_STEP_AGAIN;
    resume_3: F->pc=0; F->err=K->retval;
    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_prepare_recycle"
static int linear_step_prepare(lis_rci_native_fgcro_restore *K) {
struct linear_frame_prepare *F=&K->prepare;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}



    F->scratch_c = NULL;



















    F->state =
        F->solver->fgcrodr_state;


    if(
        F->state==NULL ||
        !F->state->valid ||
        F->state->dim<=0 )
    {
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->state->cross_action =
        LIS_FGCRODR_RECYCLE_NONE;






    if(
        F->state->matrix!=F->solver->A ||
        F->state->structure_generation!=
            F->solver->A->structure_generation )
    {
        F->err =
            lis_fgcrodr_clear_recycle(
                F->solver);

        if( F->err )
            return linear_return(K,F->parent,F->err);


        F->state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->cross_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->rejects++;
        F->state->policy_rejects++;

        F->state->quality =
            0.0;

        F->state->quality_min =
            0.0;


        return linear_return(K,F->parent,LIS_SUCCESS);
    }





    if(
        F->state->value_generation==
        F->solver->A->value_generation )
    {
        F->state->last_action =
            LIS_FGCRODR_RECYCLE_NONE;

        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->oldk =
        F->state->dim;


    F->scratch_c =
        (LIS_VECTOR *)lis_malloc(
            (size_t)F->oldk
            *sizeof(LIS_VECTOR),
            "lis_fgcrodr_prepare_recycle::scratch_c");

    if( F->scratch_c==NULL )
    {
        LIS_SETERR_MEM(
            (size_t)F->oldk
            *sizeof(LIS_VECTOR));

        return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
    }


    for(F->i=0;F->i<F->oldk;F->i++)
        F->scratch_c[F->i] = NULL;


    for(F->i=0;F->i<F->oldk;F->i++)
    {
        F->err =
            lis_vector_duplicate(
                F->solver->A,
                &F->scratch_c[F->i]);

        if( F->err )
            goto cleanup;
    }


    F->qsum =
        0.0;

    F->qmin =
        1.0;

    F->good =
        0;





    for(F->i=0;F->i<F->oldk;F->i++)
    {
        F->pc=1;
        return linear_suspend(K,LIS_RCI_OP_A,F->solver->A,F->state->z[F->i],F->scratch_c[F->i]);
        resume_1: F->pc=0; F->err=K->retval;
        if(F->err) goto cleanup;


        lis_vector_nrm2(
            F->scratch_c[F->i],
            &F->nrm);


        if(
            F->nrm<=FGCRODR_BUILD_TINY ||
            F->nrm!=F->nrm )
        {
            F->qi =
                0.0;
        }
        else
        {
            F->proj2 =
                0.0;


            for(F->j=0;F->j<F->oldk;F->j++)
            {
                lis_vector_dot(
                    F->state->c[F->j],
                    F->scratch_c[F->i],
                    &F->alpha);


                F->proj2 +=
                    (LIS_REAL)(
                        fabs(F->alpha)
                        *fabs(F->alpha));
            }


            if( F->proj2<0.0 )
                F->proj2 = 0.0;


            F->qi =
                sqrt(F->proj2)
                /F->nrm;


            if( F->qi<0.0 )
                F->qi = 0.0;

            if( F->qi>1.0 )
                F->qi = 1.0;
        }


        F->qsum +=
            F->qi;


        if( F->qi<F->qmin )
            F->qmin = F->qi;


        if(
            F->qi>=
            F->state->reject_threshold )
        {
            F->good++;
        }
    }


    F->qavg =
        F->qsum
        /(LIS_REAL)F->oldk;


    F->need =
        _max(
            (LIS_INT)1,
            (F->oldk+1)/2);


    if(
        F->qmin>=
        F->state->keep_threshold )
    {
        F->action =
            LIS_FGCRODR_RECYCLE_KEEP;
    }
    else if(
        F->qavg<
            F->state->reject_threshold ||
        F->good<F->need )
    {
        F->action =
            LIS_FGCRODR_RECYCLE_REJECT;
    }
    else
    {
        F->action =
            LIS_FGCRODR_RECYCLE_TRANSPORT;
    }


    if(
        F->action==
        LIS_FGCRODR_RECYCLE_REJECT )
    {
        F->err =
            lis_fgcrodr_clear_recycle(
                F->solver);

        if( F->err )
            goto cleanup;


        F->state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->cross_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->rejects++;
        F->state->policy_rejects++;

        F->state->quality =
            F->qavg;

        F->state->quality_min =
            F->qmin;


        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }





    F->state->quality =
        F->qavg;

    F->state->quality_min =
        F->qmin;

    F->state->cross_action =
        F->action;


    K->replace.solver = F->solver;
    K->replace.neww = F->state->w;
    K->replace.newz = F->state->z;
    K->replace.nsel = F->oldk;
    K->replace.count_update = LIS_FALSE;
    K->replace.action = F->action;
    K->replace.quality = F->qavg;
    K->replace.quality_min = F->qmin;
    K->replace.pc=0; K->replace.parent=8;
    F->pc=2; K->active=4; return L_STEP_AGAIN;
    resume_2: F->pc=0; F->err=K->retval;

    if( F->err )
        goto cleanup;


    F->state =
        F->solver->fgcrodr_state;


    if(
        !F->state->valid ||
        F->state->dim<=0 ||
        F->state->matrix!=F->solver->A ||
        F->state->structure_generation!=
            F->solver->A->structure_generation ||
        F->state->value_generation!=
            F->solver->A->value_generation )
    {
        F->err =
            lis_fgcrodr_clear_recycle(
                F->solver);

        if( F->err )
            goto cleanup;


        F->state->last_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->cross_action =
            LIS_FGCRODR_RECYCLE_REJECT;

        F->state->rejects++;
        F->state->policy_rejects++;

        F->state->quality =
            F->qavg;

        F->state->quality_min =
            F->qmin;

        F->state->cross_action =
            LIS_FGCRODR_RECYCLE_NONE;


        F->err =
            LIS_SUCCESS;

        goto cleanup;
    }


    F->state->cross_uses++;

    F->state->last_action =
        F->action;

    F->state->quality =
        F->qavg;

    F->state->quality_min =
        F->qmin;


    if(
        F->action==
        LIS_FGCRODR_RECYCLE_KEEP )
    {
        F->state->policy_keeps++;
    }
    else
    {
        F->state->policy_transports++;
    }


cleanup:

    if( F->scratch_c )
    {
        for(F->i=0;F->i<F->oldk;F->i++)
        {
            if( F->scratch_c[F->i] )
                lis_vector_destroy(
                    F->scratch_c[F->i]);
        }
    }


    lis_free(
        F->scratch_c);


    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr_flexible_core"
static int linear_step_core(lis_rci_native_fgcro_restore *K) {
struct linear_frame_core *F=&K->core;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;case 5:goto resume_5;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}

















































































    F->comm = LIS_COMM_WORLD;

    F->A = F->solver->A;
    F->b = F->solver->b;
    F->x = F->solver->x;

    F->state =
        F->solver->fgcrodr_state;

    F->use_recycle =
        (
            F->state!=NULL &&
            F->state->valid &&
            F->state->dim>0 &&
            F->state->matrix==F->A &&
            F->state->structure_generation==
                F->A->structure_generation &&
            F->state->value_generation==
                F->A->value_generation
        )
        ? LIS_TRUE
        : LIS_FALSE;

    F->kcur =
        F->use_recycle
        ? F->state->dim
        : 0;


    F->intra_recycle=F->state->intra_recycle;
    F->intra_trace=(F->intra_recycle && F->state->intra_trace);
    F->normalize_projected=F->state->projected_relative;
    F->recycle_counted=F->use_recycle;

    F->hraw = NULL;
    F->bcoef = NULL;

    F->maxiter =
        F->solver->options[LIS_OPTIONS_MAXITER];

    F->output =
        F->solver->options[LIS_OPTIONS_OUTPUT];

    F->m =
        F->solver->options[LIS_OPTIONS_RESTART];

    F->maxiter_noimp =
        F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];

    F->noimp_count = 0;
    F->min_nrm2 = LIS_SCALAR_MAX;
    F->noimp_stop = LIS_FALSE;

    F->h_dim = F->m+1;

    F->ptime = 0.0;









    F->s = F->solver->work[0];
    F->z = &F->solver->work[2];
    F->v = &F->solver->work[F->m+2];


    F->h = (LIS_SCALAR *)lis_malloc(
        sizeof(LIS_SCALAR)
        *(F->h_dim+1)
        *(F->h_dim+2),
        "lis_fgcrodr_flexible_core::h");

    if( F->h==NULL )
    {
        LIS_SETERR_MEM(
            sizeof(LIS_SCALAR)
            *(F->h_dim+1)
            *(F->h_dim+2));

        return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
    }


    F->cs = (F->m+1)*F->h_dim;
    F->sn = (F->m+2)*F->h_dim;


    F->hraw =
        (LIS_SCALAR *)lis_malloc(
            (size_t)F->h_dim
            *(size_t)F->m
            *sizeof(LIS_SCALAR),
            "lis_fgcrodr_flexible_core::hraw");

    if( F->hraw==NULL )
    {
        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

        LIS_SETERR_MEM(
            (size_t)F->h_dim
            *(size_t)F->m
            *sizeof(LIS_SCALAR));

        return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
    }


    memset(
        F->hraw,
        0,
        (size_t)F->h_dim
        *(size_t)F->m
        *sizeof(LIS_SCALAR));


    if( F->use_recycle || F->intra_recycle )
    {





        F->rows = F->intra_recycle ? F->state->target_dim : F->kcur;

        if(F->intra_recycle && F->rows>F->m-1)
            F->rows=F->m-1;
        if(F->rows<F->kcur)
            F->rows=F->kcur;
        F->capacity=(size_t)F->rows*(size_t)F->m;
        F->bcoef=(LIS_SCALAR *)lis_malloc(
            F->capacity*sizeof(LIS_SCALAR),
            "lis_fgcrodr_flexible_core::bcoef");
        if(F->bcoef==NULL)
        {
            lis_free(F->hraw);
            lis_free(F->h); F->h=NULL;
            LIS_SETERR_MEM(F->capacity*sizeof(LIS_SCALAR));
            return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
        }
        memset(F->bcoef,0,F->capacity*sizeof(LIS_SCALAR));
    }






    K->initial.solver = F->solver;
    K->initial.M = NULL;
    K->initial.t = NULL;
    K->initial.r = F->v[0];
    K->initial.bnrm2 = &F->bnrm2;
    K->initial.pc=0; K->initial.parent=9;
    F->pc=1; K->active=1; return L_STEP_AGAIN;
    resume_1: F->pc=0;
    if(K->retval)
    {
        if( F->bcoef )
            lis_free(F->bcoef);

        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

        if(F->solver->retcode!=LIS_SUCCESS)
        {
            F->solver->resid = NAN;
            lis_fgcrodr_clear_recycle(F->solver);
            return linear_return(K,F->parent,F->solver->retcode);
        }
        return linear_return(K,F->parent,LIS_SUCCESS);
    }


    F->tol = F->solver->tol;


    if( F->use_recycle )
    {
        for(F->ri=0;F->ri<F->kcur;F->ri++)
        {
            lis_vector_dot(
                F->state->c[F->ri],
                F->v[0],
                &F->gamma);

            lis_vector_axpy(
                F->gamma,
                F->state->z[F->ri],
                F->x);

            lis_vector_axpy(
                -F->gamma,
                F->state->c[F->ri],
                F->v[0]);
        }


        lis_vector_nrm2(
            F->v[0],
            &F->rnorm);


        if(
            F->rnorm==0.0 ||
            F->rnorm!=F->rnorm )
        {
            F->solver->retcode =
                (F->rnorm==0.0)
                ? LIS_SUCCESS
                : LIS_BREAKDOWN;

            F->solver->iter = 0;
            F->solver->resid = F->rnorm;
            F->solver->ptime = F->ptime;

            if( F->bcoef )
                lis_free(F->bcoef);

            if( F->hraw )
                lis_free(F->hraw);

            lis_free(F->h); F->h=NULL;

            return linear_return(K,F->parent,F->solver->retcode);
        }






        F->bnrm2 =
            1.0/F->rnorm;

        if(
            F->state->cross_action==
            LIS_FGCRODR_RECYCLE_NONE )
        {
            F->state->last_action =
                LIS_FGCRODR_RECYCLE_KEEP;
        }
        else
        {
            F->state->last_action =
                F->state->cross_action;
        }

        F->state->uses++;
    }
    else
    {
        F->rnorm =
            1.0/F->bnrm2;
    }


    F->iter = 0;
    F->nrm2 = F->rnorm * F->solver->bnrm;


    while( F->iter<F->maxiter )
    {
        F->candidate = LIS_FALSE;









        if(F->intra_recycle && F->iter>0)
        {
            F->use_recycle=(F->state->valid && F->state->dim>0 &&
                F->state->matrix==F->A &&
                F->state->structure_generation==F->A->structure_generation &&
                F->state->value_generation==F->A->value_generation);
            F->kcur=F->use_recycle ? F->state->dim : 0;
            if(F->use_recycle)
            {
                if(F->kcur>=F->m || F->kcur>F->state->target_dim)
                {
                    F->solver->retcode=LIS_ERR_ILL_ARG;
                    F->solver->iter=F->iter;
                    F->solver->ptime=F->ptime;
                    if(F->bcoef) lis_free(F->bcoef);
                    lis_free(F->hraw);
                    lis_free(F->h); F->h=NULL;
                    return linear_return(K,F->parent,LIS_ERR_ILL_ARG);
                }




                lis_vector_nrm2(F->v[0],&F->residual_before_projection);
                for(F->ri=0;F->ri<F->kcur;F->ri++)
                {
                    lis_vector_dot(F->state->c[F->ri],F->v[0],&F->gamma);
                    lis_vector_axpy(F->gamma,F->state->z[F->ri],F->x);
                    lis_vector_axpy(-F->gamma,F->state->c[F->ri],F->v[0]);
                }
                lis_vector_nrm2(F->v[0],&F->rnorm);


                if(!isfinite(F->residual_before_projection) ||
                   !isfinite(F->rnorm) ||
                   F->rnorm>F->residual_before_projection*(1.0+1.0e-8))
                {
                    F->solver->retcode=LIS_BREAKDOWN;
                    F->solver->iter=F->iter;
                    F->solver->resid=F->rnorm;
                    F->solver->ptime=F->ptime;
                    if(F->bcoef) lis_free(F->bcoef);
                    if(F->hraw) lis_free(F->hraw);
                    lis_free(F->h); F->h=NULL;
                    return linear_return(K,F->parent,LIS_BREAKDOWN);
                }
                if(F->rnorm==0.0 || F->rnorm!=F->rnorm)
                {
                    F->solver->retcode=(F->rnorm==0.0) ?
                        LIS_SUCCESS : LIS_BREAKDOWN;
                    F->solver->iter=F->iter;
                    F->solver->resid=F->rnorm;
                    F->solver->ptime=F->ptime;
                    if(F->bcoef) lis_free(F->bcoef);
                    lis_free(F->hraw);
                    lis_free(F->h); F->h=NULL;
                    return linear_return(K,F->parent,F->solver->retcode);
                }
                F->bnrm2=1.0/F->rnorm;
                memset(F->bcoef,0,
                    (size_t)F->kcur*(size_t)F->m*sizeof(LIS_SCALAR));
                if(!F->recycle_counted)
                {
                    F->state->uses++;
                    F->recycle_counted=LIS_TRUE;
                }
                if(F->intra_trace)
                    fprintf(stderr,
                        "LIS_FGCRODR_INTRA_RECYCLE iter=%d k=%d "
                        "relative=%.12e before=%.12e after=%.12e\n",
                        (int)F->iter,(int)F->kcur,
                        (double)(F->rnorm*F->solver->bnrm),
                        (double)F->residual_before_projection,
                        (double)F->rnorm);
            }
        }






        lis_vector_scale(
            F->bnrm2,
            F->v[0]);







        lis_vector_set_all(
            0,
            F->s);

        F->s->value[0] = F->rnorm;


        F->i = 0;


        do
        {
            F->iter++;
            F->i++;

            F->ii  = F->i-1;
            F->i1  = F->i;

            F->iiv = F->i-1;
            F->i1v = F->i;

            F->iih = (F->i-1)*F->h_dim;








            F->time = lis_wtime();

            F->pc=2;
            return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->v[F->iiv],F->z[F->iiv]);
            resume_2: F->pc=0; F->err=K->retval;

            F->ptime +=
                lis_wtime()-F->time;


            if( F->err )
            {
                F->solver->retcode = F->err;
                F->solver->iter = F->iter;
                F->solver->ptime = F->ptime;
                F->solver->resid = F->nrm2;

                if( F->bcoef )
            lis_free(F->bcoef);

        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

                lis_fgcrodr_clear_recycle(F->solver);
                return linear_return(K,F->parent,F->err);
            }













            F->pc=3;
            return linear_suspend(K,LIS_RCI_OP_A,F->A,F->z[F->iiv],F->v[F->i1v]);
            resume_3: F->pc=0; F->err=K->retval;
            if(F->err)
            {
                F->solver->retcode=F->err;
                F->solver->iter=F->iter;
                F->solver->ptime=F->ptime;
                F->solver->resid=F->nrm2;
                if(F->bcoef) lis_free(F->bcoef);
                if(F->hraw) lis_free(F->hraw);
                lis_free(F->h); F->h=NULL;
                lis_fgcrodr_clear_recycle(F->solver);
                return linear_return(K,F->parent,F->err);
            }








            if( F->use_recycle )
            {
                for(F->ri=0;F->ri<F->kcur;F->ri++)
                {
                    lis_vector_dot(
                        F->state->c[F->ri],
                        F->v[F->i1v],
                        &F->alpha);

                    F->bcoef[
                        F->ri
                        + F->iiv*F->kcur] =
                        F->alpha;

                    lis_vector_axpy(
                        -F->alpha,
                        F->state->c[F->ri],
                        F->v[F->i1v]);
                }
            }






            for(F->k=0;F->k<F->i;F->k++)
            {
                lis_vector_dot(
                    F->v[F->i1v],
                    F->v[F->k],
                    &F->t);

                F->h[F->k+F->iih] = F->t;

                lis_vector_axpy(
                    -F->t,
                    F->v[F->k],
                    F->v[F->i1v]);
            }






            lis_vector_nrm2(
                F->v[F->i1v],
                (LIS_REAL *)&F->t);

            F->h[F->i1+F->iih] = F->t;






            for(F->k=0;F->k<=F->i;F->k++)
            {
                F->hraw[
                    F->k
                    +F->iiv*F->h_dim] =
                    F->h[
                        F->k
                        +F->iih];
            }







            if( F->t!=0.0 && F->t==F->t )
            {
                lis_vector_scale(
                    1.0/F->t,
                    F->v[F->i1v]);
            }






            for(F->k=1;F->k<=F->ii;F->k++)
            {
                F->jj = F->k-1;

                F->t = F->h[F->jj+F->iih];

                F->aa =
                    F->h[F->jj+F->cs]*F->t;

                F->aa +=
                    F->h[F->jj+F->sn]
                    *F->h[F->k+F->iih];

                F->bb =
                    -conj(
                        F->h[F->jj+F->sn])
                    *F->t;

                F->bb +=
                    F->h[F->jj+F->cs]
                    *F->h[F->k+F->iih];

                F->h[F->jj+F->iih] = F->aa;
                F->h[F->k+F->iih] = F->bb;
            }






            F->aa = F->h[F->ii+F->iih];
            F->bb = F->h[F->i1+F->iih];




            F->aa_abs = cabs(F->aa);
            F->bb_abs = cabs(F->bb);

            F->givens_scale =
                (F->aa_abs>F->bb_abs)
                ? F->aa_abs
                : F->bb_abs;

            if( F->givens_scale==0.0 )
            {
                F->rr = 0.0;
                F->givens_c = 1.0;
                F->givens_s = 0.0;
            }
            else
            {
                F->givens_a =
                    F->aa_abs
                    /F->givens_scale;

                F->givens_b =
                    F->bb_abs
                    /F->givens_scale;

                F->givens_norm =
                    F->givens_scale
                    *sqrt(
                        F->givens_a*F->givens_a
                        +
                        F->givens_b*F->givens_b);

                F->rr =
                    (LIS_SCALAR)
                    F->givens_norm;

                if( F->aa_abs==0.0 )
                {
                    F->givens_c = 0.0;

                    F->givens_s =
                        conj(F->bb)
                        /F->bb_abs;
                }
                else
                {
                    F->givens_alpha =
                        F->aa
                        /F->aa_abs;

                    F->givens_c =
                        F->aa_abs
                        /F->givens_norm;

                    F->givens_s =
                        F->givens_alpha
                        *conj(F->bb)
                        /F->givens_norm;
                }
            }


            if( F->rr==0.0 || F->rr!=F->rr )
            {
                F->nrm2 =
                    fabs(F->s->value[F->ii]);

                F->solver->retcode =
                    LIS_BREAKDOWN;

                F->solver->iter =
                    F->iter;

                F->solver->resid =
                    F->nrm2;

                F->solver->ptime =
                    F->ptime;

                if( F->bcoef )
            lis_free(F->bcoef);

        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

                return linear_return(K,F->parent,LIS_BREAKDOWN);
            }


            F->h[F->ii+F->cs] =
                F->givens_c;

            F->h[F->ii+F->sn] =
                F->givens_s;


            F->s->value[F->i1] =
                -conj(
                    F->h[F->ii+F->sn])
                *F->s->value[F->ii];

            F->s->value[F->ii] =
                F->h[F->ii+F->cs]
                *F->s->value[F->ii];


            F->aa =
                F->h[F->ii+F->cs]
                *F->h[F->ii+F->iih];

            F->aa +=
                F->h[F->ii+F->sn]
                *F->h[F->i1+F->iih];

            F->h[F->ii+F->iih] = F->aa;

            F->h[F->i1+F->iih] =
                0.0;













            F->nrm2 =
                fabs(F->s->value[F->i1]);
            if(F->normalize_projected)
                F->nrm2 *= F->solver->bnrm;


            if( F->output )
            {
                if( F->output & LIS_PRINT_MEM )
                {
                    F->solver->rhistory[F->iter] =
                        F->nrm2;
                }

                if( F->output & LIS_PRINT_OUT )
                {
                    lis_print_rhistory(
                        F->comm,
                        F->iter,
                        F->nrm2);
                }
            }


            lis_solver_monitor_residual(
                F->solver,F->iter,F->nrm2,F->tol,
                LIS_SOLVER_RESIDUAL_PROJECTED);

            F->candidate = lis_solver_test_convergence(
                F->solver,F->iter,F->nrm2,F->tol,
                LIS_SOLVER_RESIDUAL_PROJECTED,
                F->tol>=F->nrm2);

            if( F->candidate )
                break;


            if( F->maxiter_noimp )
            {
                if( F->min_nrm2>F->nrm2 )
                {
                    F->min_nrm2 =
                        F->nrm2;

                    F->noimp_count =
                        0;
                }
                else if(
                    ++F->noimp_count >
                    F->maxiter_noimp )
                {





                    F->noimp_stop =
                        LIS_TRUE;

                    break;
                }
            }


        } while(
            F->i<F->m &&
            F->iter<F->maxiter );






        if(
            F->h[F->ii+F->iih]==0.0 ||
            F->h[F->ii+F->iih]!=F->h[F->ii+F->iih] )
        {
            F->solver->retcode =
                LIS_BREAKDOWN;

            F->solver->iter =
                F->iter;

            F->solver->resid =
                F->nrm2;

            F->solver->ptime =
                F->ptime;

            if( F->bcoef )
            lis_free(F->bcoef);

        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

            return linear_return(K,F->parent,LIS_BREAKDOWN);
        }


        F->s->value[F->ii] =
            F->s->value[F->ii]
            /F->h[F->ii+F->iih];







        if(
            F->s->value[F->ii]!=F->s->value[F->ii] ||
            fabs(F->s->value[F->ii])>LIS_SCALAR_MAX )
        {
            F->solver->retcode =
                LIS_BREAKDOWN;

            F->solver->iter =
                F->iter;

            F->solver->resid =
                F->nrm2;

            F->solver->ptime =
                F->ptime;

            if( F->bcoef )
                lis_free(
                    F->bcoef);

            if( F->hraw )
                lis_free(
                    F->hraw);

            lis_free(
                F->h);

            return linear_return(K,F->parent,LIS_BREAKDOWN);
        }


        for(F->k=1;F->k<=F->ii;F->k++)
        {
            F->jj = F->ii-F->k;

            F->t =
                F->s->value[F->jj];


            for(F->j=F->jj+1;F->j<=F->ii;F->j++)
            {
                F->t -=
                    F->h[F->jj+F->j*F->h_dim]
                    *F->s->value[F->j];
            }


            if(
                F->h[F->jj+F->jj*F->h_dim]==0.0 ||
                F->h[F->jj+F->jj*F->h_dim]
                    !=F->h[F->jj+F->jj*F->h_dim] )
            {
                F->solver->retcode =
                    LIS_BREAKDOWN;

                F->solver->iter =
                    F->iter;

                F->solver->resid =
                    F->nrm2;

                F->solver->ptime =
                    F->ptime;

                if( F->bcoef )
            lis_free(F->bcoef);

        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

                return linear_return(K,F->parent,LIS_BREAKDOWN);
            }


            F->s->value[F->jj] =
                F->t
                /F->h[F->jj+F->jj*F->h_dim];







            if(
                F->s->value[F->jj]!=F->s->value[F->jj] ||
                fabs(F->s->value[F->jj])>LIS_SCALAR_MAX )
            {
                F->solver->retcode =
                    LIS_BREAKDOWN;

                F->solver->iter =
                    F->iter;

                F->solver->resid =
                    F->nrm2;

                F->solver->ptime =
                    F->ptime;

                if( F->bcoef )
                    lis_free(
                        F->bcoef);

                if( F->hraw )
                    lis_free(
                        F->hraw);

                lis_free(
                    F->h);

                return linear_return(K,F->parent,LIS_BREAKDOWN);
            }
        }










        for(F->j=0;F->j<=F->ii;F->j++)
        {
            lis_vector_axpy(
                F->s->value[F->j],
                F->z[F->j],
                F->x);
        }


        if( F->use_recycle )
        {
            for(F->ri=0;F->ri<F->kcur;F->ri++)
            {
                F->gamma = 0.0;

                for(F->j=0;F->j<=F->ii;F->j++)
                {
                    F->gamma +=
                        F->bcoef[
                            F->ri
                            + F->j*F->kcur]
                        *F->s->value[F->j];
                }

                lis_vector_axpy(
                    -F->gamma,
                    F->state->z[F->ri],
                    F->x);
            }
        }










        K->build.solver = F->solver;
        K->build.v = F->v;
        K->build.z = F->z;
        K->build.jdim = F->ii+1;
        K->build.hraw = F->hraw;
        K->build.hld = F->h_dim;
        K->build.bcoef = F->bcoef;
        K->build.bld = F->kcur;
        K->build.pc=0; K->build.parent=9;
        F->pc=4; K->active=7; return L_STEP_AGAIN;
        resume_4: F->pc=0; F->err=K->retval;

        if( F->err )
        {
            F->solver->retcode =
                F->err;

            F->solver->iter =
                F->iter;

            F->solver->resid =
                F->nrm2;

            F->solver->ptime =
                F->ptime;

            if( F->bcoef )
            lis_free(F->bcoef);

        if( F->hraw )
            lis_free(F->hraw);

        lis_free(F->h); F->h=NULL;

            lis_fgcrodr_clear_recycle(F->solver);
            return linear_return(K,F->parent,F->err);
        }







        if( F->noimp_stop )
        {
            F->solver->retcode =
                LIS_MAXITER;

            F->solver->iter =
                F->iter;

            F->solver->resid =
                F->nrm2;

            F->solver->ptime =
                F->ptime;

            if( F->bcoef )
                lis_free(
                    F->bcoef);

            if( F->hraw )
                lis_free(
                    F->hraw);

            lis_free(
                F->h);

            return linear_return(K,F->parent,LIS_MAXITER);
        }










        F->pc=5;
        return linear_suspend(K,LIS_RCI_OP_A,F->A,F->x,F->v[0]);
        resume_5: F->pc=0; F->err=K->retval;
        if(!F->err) F->err=lis_vector_xpay(
            F->b,
            -1.0,
            F->v[0]);
        if(F->err)
        {
            F->solver->retcode=F->err;
            F->solver->iter=F->iter;
            F->solver->ptime=F->ptime;
            F->solver->resid=NAN;
            if(F->bcoef) lis_free(F->bcoef);
            if(F->hraw) lis_free(F->hraw);
            lis_free(F->h); F->h=NULL;
            lis_fgcrodr_clear_recycle(F->solver);
            return linear_return(K,F->parent,F->err);
        }


        if( F->candidate )
        {
            lis_solver_get_residual[0](
                F->v[0],
                F->solver,
                &F->nrm2);


            if(
                F->nrm2!=F->nrm2 ||
                fabs(F->nrm2)>LIS_SCALAR_MAX )
            {
                F->solver->retcode =
                    LIS_BREAKDOWN;

                F->solver->iter =
                    F->iter;

                F->solver->resid =
                    F->nrm2;

                F->solver->ptime =
                    F->ptime;

                if( F->bcoef )
                    lis_free(
                        F->bcoef);

                if( F->hraw )
                    lis_free(
                        F->hraw);

                lis_free(
                    F->h);

                return linear_return(K,F->parent,LIS_BREAKDOWN);
            }


            lis_solver_monitor_residual(
                F->solver,
                F->iter,
                F->nrm2,
                F->tol,
                LIS_SOLVER_RESIDUAL_TRUE);

            if( lis_solver_test_convergence(
                    F->solver,
                    F->iter,
                    F->nrm2,
                    F->tol,
                    LIS_SOLVER_RESIDUAL_TRUE,
                    F->tol>=F->nrm2) )
            {
                F->solver->retcode =
                    LIS_SUCCESS;

                F->solver->iter =
                    F->iter;

                F->solver->resid =
                    F->nrm2;

                F->solver->ptime =
                    F->ptime;

                if( F->bcoef )
                    lis_free(
                        F->bcoef);

                if( F->hraw )
                    lis_free(
                        F->hraw);

                lis_free(
                    F->h);

                return linear_return(K,F->parent,LIS_SUCCESS);
            }
        }










        if( F->use_recycle && !F->intra_recycle )
        {
            for(F->ri=0;F->ri<F->kcur;F->ri++)
            {
                lis_vector_dot(
                    F->state->c[F->ri],
                    F->v[0],
                    &F->gamma);

                lis_vector_axpy(
                    F->gamma,
                    F->state->z[F->ri],
                    F->x);

                lis_vector_axpy(
                    -F->gamma,
                    F->state->c[F->ri],
                    F->v[0]);
            }
        }


        lis_vector_nrm2(
            F->v[0],
            &F->rnorm);











        if(
            F->rnorm==0.0 ||
            F->rnorm!=F->rnorm )
        {
            F->solver->retcode =
                (F->rnorm==0.0)
                ? LIS_SUCCESS
                : LIS_BREAKDOWN;

            F->solver->iter =
                F->iter;

            F->solver->resid =
                F->rnorm;

            F->solver->ptime =
                F->ptime;

            if( F->bcoef )
                lis_free(
                    F->bcoef);

            if( F->hraw )
                lis_free(
                    F->hraw);

            lis_free(
                F->h);


            return linear_return(K,F->parent,F->solver->retcode);
        }


        F->bnrm2 =
            1.0/F->rnorm;
    }


    F->solver->retcode =
        LIS_MAXITER;

    F->solver->iter =
        F->iter+1;

    F->solver->resid =
        F->nrm2;



    if( F->bcoef )
        lis_free(F->bcoef);
    if( F->hraw )
        lis_free(F->hraw);

    lis_free(F->h); F->h=NULL;

    return linear_return(K,F->parent,LIS_MAXITER);

}
#undef __FUNC__
#define __FUNC__ "lis_fgcrodr"
static int linear_step_fgcro(lis_rci_native_fgcro_restore *K) {
struct linear_frame_fgcro *F=&K->fgcro;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}






    LIS_DEBUG_FUNC_IN;


    F->err =
        lis_fgcrodr_state_ensure(
            F->solver);

    if( F->err )
        return linear_return(K,F->parent,F->err);


    F->state =
        F->solver->fgcrodr_state;


    F->solve_cross_action =
        LIS_FGCRODR_RECYCLE_NONE;






    if( F->state->target_dim==0 )
    {
        LIS_DEBUG_FUNC_OUT;

        return linear_return(K,F->parent,LIS_ERR_ILL_ARG);
    }






    K->prepare.solver = F->solver;
    K->prepare.pc=0; K->prepare.parent=10;
    F->pc=1; K->active=8; return L_STEP_AGAIN;
    resume_1: F->pc=0; F->err=K->retval;

    if( F->err )
    {
        F->solver->retcode = F->err;
        F->solver->iter = 0;
        F->solver->resid = NAN;
        F->solver->ptime = 0.0;
        lis_fgcrodr_clear_recycle(F->solver);
        LIS_DEBUG_FUNC_OUT;
        return linear_return(K,F->parent,F->err);
    }













    if( F->state->valid &&
        F->state->dim >= F->solver->options[LIS_OPTIONS_RESTART] )
    {
        F->err = lis_fgcrodr_clear_recycle(F->solver);
        if( F->err )
            return linear_return(K,F->parent,F->err);

        F->state->rejects++;
        F->state->policy_rejects++;
        F->state->last_action = LIS_FGCRODR_RECYCLE_REJECT;
        F->state->cross_action = LIS_FGCRODR_RECYCLE_REJECT;
        F->solve_cross_action = LIS_FGCRODR_RECYCLE_REJECT;
    }
    else
    {
        F->solve_cross_action =
            F->state->cross_action;
    }


    K->core.solver = F->solver;
    K->core.pc=0; K->core.parent=10;
    F->pc=2; K->active=9; return L_STEP_AGAIN;
    resume_2: F->pc=0; F->err=K->retval;









    if(
        F->solve_cross_action
        !=LIS_FGCRODR_RECYCLE_NONE )
    {
        F->state->last_action =
            F->solve_cross_action;
    }


    F->state->cross_action =
        LIS_FGCRODR_RECYCLE_NONE;


    LIS_DEBUG_FUNC_OUT;

    return linear_return(K,F->parent,F->err);

}
static int linear_step(lis_rci_native_fgcro_restore *K){switch(K->active){
case 1:return linear_step_initial(K);
case 2:return linear_step_install(K);
case 3:return linear_step_harmonic_build(K);
case 4:return linear_step_replace(K);
case 5:return linear_step_commit(K);
case 6:return linear_step_harmonic_update(K);
case 7:return linear_step_build(K);
case 8:return linear_step_prepare(K);
case 9:return linear_step_core(K);
case 10:return linear_step_fgcro(K);
default:return L_STEP_STOP;}}
