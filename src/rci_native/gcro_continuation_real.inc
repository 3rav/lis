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
struct linear_frame_refresh { int pc,parent;
    LIS_SOLVER solver;
    LIS_INT i;
    LIS_INT j;
    LIS_INT pass;
    LIS_INT keep;
    LIS_INT err;
    LIS_SCALAR alpha;
    LIS_REAL nrm;
    LIS_VECTOR ui;
    LIS_VECTOR zui;
    LIS_VECTOR ci;
};
struct linear_frame_assess { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * scratch_z;
    LIS_VECTOR * scratch_c;
    LIS_INT i;
    LIS_INT j;
    LIS_INT good;
    LIS_INT need;
    LIS_INT err;
    LIS_INT action;
    LIS_INT oldk;
    LIS_SCALAR alpha;
    LIS_REAL nrm;
    LIS_REAL proj2;
    LIS_REAL qi;
    LIS_REAL qsum;
    LIS_REAL qmin;
};
struct linear_frame_update { int pc,parent;
    LIS_SOLVER solver;
    LIS_VECTOR * v;
    LIS_INT jdim;
    const LIS_SCALAR * hraw;
    LIS_INT hld;
    const LIS_SCALAR * bcoef;
    LIS_INT bld;
    LIS_INT oldk;
    LIS_INT kmax;
    LIS_INT knew;
    LIS_INT p;
    LIS_INT q;
    LIS_INT i;
    LIS_INT j;
    LIS_INT l;
    LIS_INT r;
    LIS_INT c;
    LIS_INT err;
    LIS_INT qriter;
    LIS_INT neig;
    LIS_INT nsel;
    LIS_INT eidx;
    LIS_REAL qrerr;
    LIS_SCALAR * G;
    LIS_SCALAR * F;
    LIS_SCALAR * AP;
    LIS_SCALAR * BP;
    LIS_SCALAR * T;
    LIS_SCALAR * TQR;
    LIS_SCALAR * QW;
    LIS_SCALAR * RW;
    LIS_SCALAR * AINV;
    LIS_GCRODR_EIG * eig;
    LIS_REAL * vr;
    LIS_REAL * vi;
    LIS_VECTOR * newu;
    LIS_SCALAR ag;
    LIS_SCALAR bg;
    LIS_SCALAR s;
};
struct linear_frame_gcro { int pc,parent;
    LIS_SOLVER solver;
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR s;
    LIS_VECTOR * z;
    LIS_VECTOR * v;
    LIS_VECTOR wtmp;
    LIS_SCALAR * h;
    LIS_SCALAR * hraw;
    LIS_SCALAR * bcoef;
    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR rr;
    LIS_SCALAR t;
    LIS_SCALAR alpha;
    LIS_SCALAR gamma;
    LIS_REAL bnrm;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_REAL rnorm;
    LIS_REAL wnorm;
    LIS_SCALAR a2;
    LIS_SCALAR b2;
    LIS_INT err;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT output;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    LIS_INT i;
    LIS_INT j;
    LIS_INT k;
    LIS_INT m;
    LIS_INT kmax;
    LIS_INT kcur;
    LIS_INT cycle_steps;
    LIS_INT jdim;
    LIS_INT ii;
    LIS_INT i1;
    LIS_INT jj;
    LIS_INT hld;
    LIS_INT cs;
    LIS_INT sn;
    LIS_INT pass;
    LIS_INT lucky;
    LIS_INT candidate;
    double time;
    double ptime;
};
struct lis_rci_native_gcro_restore { int active,method,status,strict; LIS_INT retval,native_code; uint64_t token,matvecs,precons; LIS_SOLVER solver; lis_rci_native_linear_request request;
    struct linear_frame_initial initial;
    struct linear_frame_refresh refresh;
    struct linear_frame_assess assess;
    struct linear_frame_update update;
    struct linear_frame_gcro gcro;
};
#define L_STEP_AGAIN 1
#define L_STEP_YIELD 2
#define L_STEP_STOP 3
static int linear_suspend(lis_rci_native_gcro_restore *,uint32_t,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);
static int linear_return(lis_rci_native_gcro_restore *,int,LIS_INT);
#undef __FUNC__
#define __FUNC__ "lis_solver_get_initial_residual"
static int linear_step_initial(lis_rci_native_gcro_restore *K) {
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
#define __FUNC__ "lis_gcrodr_refresh_recycle"
static int linear_step_refresh(lis_rci_native_gcro_restore *K) {
struct linear_frame_refresh *F=&K->refresh;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}






    if( !F->solver->recycle_valid || F->solver->recycle_dim<=0 ) return linear_return(K,F->parent,LIS_SUCCESS);
    if( !lis_gcrodr_recycle_compatible(F->solver) )
    {
        F->solver->recycle_rejects++;
        F->solver->recycle_dim = 0;
        F->solver->recycle_valid = LIS_FALSE;
        return linear_return(K,F->parent,LIS_SUCCESS);
    }

    F->keep = 0;
    for(F->i=0;F->i<F->solver->recycle_dim;F->i++)
    {
        F->ui = F->solver->recycle_u[F->i];
        F->zui = F->solver->recycle_zu[F->i];
        F->ci = F->solver->recycle_c[F->i];

        F->pc=1;
        return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->ui,F->zui);
        resume_1: F->pc=0; F->err=K->retval;
        if( F->err ) return linear_return(K,F->parent,F->err);
        F->pc=2;
        return linear_suspend(K,LIS_RCI_OP_A,F->solver->A,F->zui,F->ci);
        resume_2: F->pc=0; F->err=K->retval;
        if( F->err ) return linear_return(K,F->parent,F->err);

        for(F->pass=0;F->pass<2;F->pass++)
        {
            for(F->j=0;F->j<F->keep;F->j++)
            {
                lis_vector_dot(F->solver->recycle_c[F->j],F->ci,&F->alpha);
                lis_vector_axpy(-F->alpha,F->solver->recycle_c[F->j],F->ci);
                lis_vector_axpy(-F->alpha,F->solver->recycle_u[F->j],F->ui);
                lis_vector_axpy(-F->alpha,F->solver->recycle_zu[F->j],F->zui);
            }
        }

        lis_vector_nrm2(F->ci,&F->nrm);
        if( F->nrm<=GCRODR_REORTH_TINY || F->nrm!=F->nrm )
        {
            F->solver->recycle_rejects++;
            continue;
        }

        lis_vector_scale(1.0/F->nrm,F->ci);
        lis_vector_scale(1.0/F->nrm,F->ui);
        lis_vector_scale(1.0/F->nrm,F->zui);

        if( F->keep!=F->i )
        {
            lis_vector_copy(F->ui,F->solver->recycle_u[F->keep]);
            lis_vector_copy(F->zui,F->solver->recycle_zu[F->keep]);
            lis_vector_copy(F->ci,F->solver->recycle_c[F->keep]);
        }
        F->keep++;
    }

    F->solver->recycle_dim = F->keep;
    F->solver->recycle_valid = (F->keep>0) ? LIS_TRUE : LIS_FALSE;
    if( F->keep>0 ) F->solver->recycle_uses++;
    return linear_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_gcrodr_assess_refresh_recycle"
static int linear_step_assess(lis_rci_native_gcro_restore *K) {
struct linear_frame_assess *F=&K->assess;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}





    if( !F->solver->recycle_valid || F->solver->recycle_dim<=0 ) return linear_return(K,F->parent,LIS_SUCCESS);

    if( !lis_gcrodr_recycle_compatible(F->solver) )
    {
        F->solver->recycle_policy_rejects++;
        F->solver->recycle_rejects++;
        F->solver->recycle_last_action = LIS_GCRODR_RECYCLE_REJECT;
        F->solver->recycle_quality = 0.0;
        F->solver->recycle_quality_min = 0.0;
        F->solver->recycle_dim = 0;
        F->solver->recycle_valid = LIS_FALSE;
        return linear_return(K,F->parent,LIS_SUCCESS);
    }

    F->oldk = F->solver->recycle_dim;

    F->qsum = 0.0;
    F->qmin = 1.0;
    F->good = 0;


    for(F->i=0;F->i<F->oldk;F->i++)
    {
        F->pc=1;
        return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->solver->recycle_u[F->i],F->scratch_z[F->i]);
        resume_1: F->pc=0; F->err=K->retval;
        if( F->err )
        {
            return linear_return(K,F->parent,F->err);
        }
        F->pc=2;
        return linear_suspend(K,LIS_RCI_OP_A,F->solver->A,F->scratch_z[F->i],F->scratch_c[F->i]);
        resume_2: F->pc=0; F->err=K->retval;
        if( F->err ) return linear_return(K,F->parent,F->err);
        lis_vector_nrm2(F->scratch_c[F->i],&F->nrm);

        if( F->nrm<=GCRODR_REORTH_TINY || F->nrm!=F->nrm )
        {
            F->qi = 0.0;
        }
        else
        {
            F->proj2 = 0.0;
            for(F->j=0;F->j<F->oldk;F->j++)
            {
                lis_vector_dot(F->solver->recycle_c[F->j],F->scratch_c[F->i],&F->alpha);
                F->proj2 += lis_gcrodr_abs2_scalar(F->alpha);
            }
            if( F->proj2<0.0 ) F->proj2 = 0.0;
            F->qi = sqrt(F->proj2)/F->nrm;
            if( F->qi<0.0 ) F->qi = 0.0;
            if( F->qi>1.0 ) F->qi = 1.0;
        }

        F->qsum += F->qi;
        if( F->qi<F->qmin ) F->qmin = F->qi;
        if( F->qi>=F->solver->recycle_reject_threshold ) F->good++;
    }

    F->solver->recycle_quality = F->qsum/(LIS_REAL)F->oldk;
    F->solver->recycle_quality_min = F->qmin;

    F->need = _max((LIS_INT)1,(F->oldk+1)/2);

    if( F->qmin>=F->solver->recycle_keep_threshold )
    {
        F->action = LIS_GCRODR_RECYCLE_KEEP;
    }
    else if( F->solver->recycle_quality<F->solver->recycle_reject_threshold ||
             F->good<F->need )
    {
        F->action = LIS_GCRODR_RECYCLE_REJECT;
    }
    else
    {
        F->action = LIS_GCRODR_RECYCLE_REFRESH;
    }

    if( F->action==LIS_GCRODR_RECYCLE_REJECT )
    {
        F->solver->recycle_policy_rejects++;
        F->solver->recycle_rejects++;
        F->solver->recycle_last_action = F->action;
        F->solver->recycle_dim = 0;
        F->solver->recycle_valid = LIS_FALSE;
        return linear_return(K,F->parent,LIS_SUCCESS);
    }











    K->refresh.solver = F->solver;
    K->refresh.pc=0; K->refresh.parent=3;
    F->pc=3; K->active=2; return L_STEP_AGAIN;
    resume_3: F->pc=0; F->err=K->retval;
    if( F->err )
    {
        lis_gcrodr_clear_recycle(F->solver);
        return linear_return(K,F->parent,F->err);
    }

    if( !F->solver->recycle_valid || F->solver->recycle_dim<=0 )
    {
        F->solver->recycle_policy_rejects++;
        F->solver->recycle_last_action = LIS_GCRODR_RECYCLE_REJECT;
        F->solver->recycle_quality = 0.0;
        F->solver->recycle_quality_min = 0.0;
        return linear_return(K,F->parent,LIS_SUCCESS);
    }

    F->solver->recycle_cross_uses++;
    F->solver->recycle_last_action = F->action;

    if( F->action==LIS_GCRODR_RECYCLE_KEEP )
        F->solver->recycle_policy_keeps++;
    else
        F->solver->recycle_policy_refreshes++;

    return linear_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_gcrodr_update_recycle"
static int linear_step_update(lis_rci_native_gcro_restore *K) {
struct linear_frame_update *F=&K->update;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}









    F->kmax = _min(F->solver->recycle_target_dim, F->solver->options[LIS_OPTIONS_RESTART]-1);
    F->oldk = F->solver->recycle_valid ? F->solver->recycle_dim : 0;
    if( F->kmax<=0 || F->jdim<=0 ) return linear_return(K,F->parent,LIS_SUCCESS);

    F->p = F->oldk+F->jdim;
    F->q = F->oldk+F->jdim+1;
    if( F->p<=0 ) return linear_return(K,F->parent,LIS_SUCCESS);
    F->knew = _min(F->kmax,F->p);

    F->G    = (LIS_SCALAR *)lis_malloc((size_t)F->q*F->p*sizeof(LIS_SCALAR),"gcrodr::G");
    F->F    = (LIS_SCALAR *)lis_malloc((size_t)F->q*F->p*sizeof(LIS_SCALAR),"gcrodr::F");
    F->AP   = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::AP");
    F->BP   = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::BP");
    F->T    = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::T");
    F->TQR  = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::TQR");
    F->QW   = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::QW");
    F->RW   = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::RW");
    F->AINV = (LIS_SCALAR *)lis_malloc((size_t)F->p*F->p*sizeof(LIS_SCALAR),"gcrodr::AINV");
    F->eig  = (LIS_GCRODR_EIG *)lis_malloc((size_t)F->p*sizeof(*F->eig),"gcrodr::eig");
    F->vr   = (LIS_REAL *)lis_malloc((size_t)F->p*sizeof(LIS_REAL),"gcrodr::vr");
    F->vi   = (LIS_REAL *)lis_malloc((size_t)F->p*sizeof(LIS_REAL),"gcrodr::vi");
    F->newu = (LIS_VECTOR *)lis_malloc((size_t)F->knew*sizeof(LIS_VECTOR),"gcrodr::newu");
    if( F->newu!=NULL )
        for(F->i=0;F->i<F->knew;F->i++) F->newu[F->i]=NULL;

    if( F->G==NULL || F->F==NULL || F->AP==NULL || F->BP==NULL || F->T==NULL || F->TQR==NULL ||
        F->QW==NULL || F->RW==NULL || F->AINV==NULL ||
        F->eig==NULL || F->vr==NULL || F->vi==NULL ||
        F->newu==NULL )
    {
        F->err = LIS_ERR_OUT_OF_MEMORY;
        goto cleanup;
    }

    memset(F->G,0,(size_t)F->q*F->p*sizeof(LIS_SCALAR));
    memset(F->F,0,(size_t)F->q*F->p*sizeof(LIS_SCALAR));


    for(F->i=0;F->i<F->oldk;F->i++) F->G[F->i+F->i*F->q] = 1.0;
    for(F->j=0;F->j<F->jdim;F->j++)
    {
        for(F->i=0;F->i<F->oldk;F->i++) F->G[F->i+(F->oldk+F->j)*F->q] = F->bcoef[F->i+F->j*F->bld];
        for(F->i=0;F->i<=F->jdim;F->i++) F->G[(F->oldk+F->i)+(F->oldk+F->j)*F->q] = F->hraw[F->i+F->j*F->hld];
    }




    for(F->c=0;F->c<F->oldk;F->c++)
    {
        for(F->r=0;F->r<F->oldk;F->r++)
            lis_vector_dot(F->solver->recycle_c[F->r],F->solver->recycle_u[F->c],
                           &F->F[F->r+F->c*F->q]);
        for(F->r=0;F->r<=F->jdim;F->r++)
            lis_vector_dot(F->v[F->r],F->solver->recycle_u[F->c],
                           &F->F[(F->oldk+F->r)+F->c*F->q]);
    }
    for(F->c=0;F->c<F->jdim;F->c++) F->F[(F->oldk+F->c)+(F->oldk+F->c)*F->q] = 1.0;


    for(F->c=0;F->c<F->p;F->c++)
    {
        for(F->r=0;F->r<F->p;F->r++)
        {
            F->ag = 0.0;
            F->bg = 0.0;
            for(F->l=0;F->l<F->q;F->l++)
            {
                F->ag +=
                    lis_gcrodr_conj_scalar(
                        F->G[F->l+F->r*F->q])
                    *F->G[F->l+F->c*F->q];

                F->bg +=
                    lis_gcrodr_conj_scalar(
                        F->G[F->l+F->r*F->q])
                    *F->F[F->l+F->c*F->q];
            }
            F->AP[F->r+F->c*F->p] = F->ag;
            F->BP[F->r+F->c*F->p] = F->bg;
        }
    }

    if(
        !lis_gcrodr_dense_finite(
            F->p,
            F->AP)
        ||
        !lis_gcrodr_dense_finite(
            F->p,
            F->BP) )
    {
        goto reject;
    }

    memcpy(F->AINV,F->AP,(size_t)F->p*F->p*sizeof(LIS_SCALAR));

    F->err = lis_array_ge(F->p,F->AINV);

    if(
        F->err
        ||
        !lis_gcrodr_dense_finite(
            F->p,
            F->AINV) )
    {
        goto reject;
    }


    for(F->c=0;F->c<F->p;F->c++)
    {
        for(F->r=0;F->r<F->p;F->r++)
        {
            F->s = 0.0;
            for(F->l=0;F->l<F->p;F->l++) F->s += F->AINV[F->r+F->l*F->p]*F->BP[F->l+F->c*F->p];
            F->T[F->r+F->c*F->p] = F->s;
        }
    }

    if(
        !lis_gcrodr_dense_finite(
            F->p,
            F->T) )
    {
        goto reject;
    }

    memcpy(F->TQR,F->T,(size_t)F->p*F->p*sizeof(LIS_SCALAR));





    if( F->oldk>0 )
        F->qriter = 512;
    else
        F->qriter = 1024;

    F->err = lis_array_hqr_limit(
        F->p,
        F->TQR,
        F->QW,
        F->RW,
        F->qriter,
        &F->qriter,
        &F->qrerr);
    if( F->err ) goto reject;
    (void)F->qriter;
    (void)F->qrerr;

    F->err = lis_gcrodr_extract_eigs(F->p,F->TQR,F->eig,&F->neig);
    if( F->err ) goto reject;


    for(F->i=0;F->i<F->knew;F->i++)
    {
        F->err = lis_vector_duplicate(F->solver->A,&F->newu[F->i]);
        if( F->err ) goto cleanup;
    }

    F->nsel = 0;
    for(F->eidx=0;F->eidx<F->neig && F->nsel<F->knew;F->eidx++)
    {
        if( F->eig[F->eidx].pair && F->nsel+2>F->knew ) continue;
        F->err = lis_gcrodr_inverse_eigenvector(F->p,F->T,F->eig[F->eidx].lr,F->eig[F->eidx].li,F->vr,F->vi);
        if( F->err==LIS_ERR_OUT_OF_MEMORY ) goto cleanup;
        if( F->err ) continue;

        if( F->eig[F->eidx].pair )
        {
            F->err = lis_gcrodr_build_vector(F->solver,F->v,F->jdim,F->oldk,F->vr,F->newu[F->nsel]);
            if( F->err ) continue;
            F->err = lis_gcrodr_build_vector(F->solver,F->v,F->jdim,F->oldk,F->vi,F->newu[F->nsel+1]);
            if( F->err ) continue;
            F->nsel += 2;
        }
        else
        {
            F->err = lis_gcrodr_build_vector(F->solver,F->v,F->jdim,F->oldk,F->vr,F->newu[F->nsel]);
            if( F->err ) continue;
            F->nsel++;
        }
    }


    if( F->nsel<=0 ) goto reject;

    for(F->i=0;F->i<F->nsel;F->i++) lis_vector_copy(F->newu[F->i],F->solver->recycle_u[F->i]);
    F->solver->recycle_dim = F->nsel;
    F->solver->recycle_valid = LIS_TRUE;
    if( F->oldk>0 )
    {
        F->solver->recycle_updates++;
    }
    else
    {
        if( F->solver->recycle_builds==0 &&
            F->solver->recycle_last_action==LIS_GCRODR_RECYCLE_NONE )
            F->solver->recycle_last_action = LIS_GCRODR_RECYCLE_BUILD;
        F->solver->recycle_builds++;
    }

    K->refresh.solver = F->solver;
    K->refresh.pc=0; K->refresh.parent=4;
    F->pc=1; K->active=2; return L_STEP_AGAIN;
    resume_1: F->pc=0; F->err=K->retval;
    if( F->err )
    {
        lis_gcrodr_clear_recycle(F->solver);
        goto cleanup;
    }
    goto cleanup;

reject:
    F->solver->recycle_rejects++;
    F->err = LIS_SUCCESS;

cleanup:
    if( F->newu )
    {
        for(F->i=0;F->i<F->knew;F->i++) lis_vector_destroy(F->newu[F->i]);
    }
    lis_free(F->newu);
    lis_free(F->G); lis_free(F->F); lis_free(F->AP); lis_free(F->BP); lis_free(F->T);
    lis_free(F->TQR); lis_free(F->QW); lis_free(F->RW); lis_free(F->AINV);
    lis_free(F->eig); lis_free(F->vr); lis_free(F->vi);
    return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_gcrodr"
static int linear_step_gcro(lis_rci_native_gcro_restore *K) {
struct linear_frame_gcro *F=&K->gcro;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;case 5:goto resume_5;case 6:goto resume_6;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}













    LIS_DEBUG_FUNC_IN;
    F->comm = LIS_COMM_WORLD;
    F->A = F->solver->A;
    F->b = F->solver->b;
    F->x = F->solver->x;
    F->maxiter = F->solver->options[LIS_OPTIONS_MAXITER];
    F->output = F->solver->options[LIS_OPTIONS_OUTPUT];
    F->m = F->solver->options[LIS_OPTIONS_RESTART];
    if( F->solver->recycle_target_dim<0 )
        F->solver->recycle_target_dim = GCRODR_DEFAULT_RECYCLE_DIM;
    F->kmax = _min(F->solver->recycle_target_dim, F->m-1);
    F->maxiter_noimp = F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
    F->min_nrm2 = LIS_SCALAR_MAX;
    F->noimp_count = 0;
    F->ptime = 0.0;
    F->iter = 0;
    F->nrm2 = NAN;

    F->s = F->solver->work[0];
    F->wtmp = F->solver->work[1];
    F->z = &F->solver->work[2];
    F->v = &F->solver->work[F->m+2];
    F->hld = F->m+1;

    F->h = (LIS_SCALAR *)lis_malloc((size_t)(F->hld+1)*(F->hld+2)*sizeof(LIS_SCALAR),
                                 "lis_gcrodr::h");
    F->hraw = (LIS_SCALAR *)lis_malloc((size_t)(F->m+1)*F->m*sizeof(LIS_SCALAR),
                                    "lis_gcrodr::hraw");
    F->bcoef = (LIS_SCALAR *)lis_malloc((size_t)_max(1,F->kmax)*F->m*sizeof(LIS_SCALAR),
                                     "lis_gcrodr::bcoef");
    if( F->h==NULL || F->hraw==NULL || F->bcoef==NULL )
    {
        LIS_SETERR_MEM(sizeof(LIS_SCALAR));
        F->err = LIS_ERR_OUT_OF_MEMORY;
        goto fail;
    }
    F->cs = (F->m+1)*F->hld;
    F->sn = (F->m+2)*F->hld;

    F->err = lis_gcrodr_ensure_recycle(F->solver,F->kmax);
    if( F->err ) goto fail;





    if( F->solver->recycle_valid && F->solver->recycle_dim>0 )
    {
        K->assess.solver = F->solver;
        K->assess.scratch_z = F->z;
        K->assess.scratch_c = F->v;
        K->assess.pc=0; K->assess.parent=5;
        F->pc=1; K->active=3; return L_STEP_AGAIN;
        resume_1: F->pc=0; F->err=K->retval;
        if( F->err ) goto fail;
    }


    K->initial.solver = F->solver;
    K->initial.M = NULL;
    K->initial.t = NULL;
    K->initial.r = F->v[0];
    K->initial.bnrm2 = &F->bnrm;
    K->initial.pc=0; K->initial.parent=5;
    F->pc=2; K->active=1; return L_STEP_AGAIN;
    resume_2: F->pc=0;
    if(K->retval)
    {
        F->err = F->solver->retcode;
        if(F->err!=LIS_SUCCESS) goto fail;
        goto done;
    }
    F->tol = F->solver->tol;
    if(F->solver->options[LIS_OPTIONS_CONV_COND]==LIS_CONV_COND_NRM2_R &&
        isfinite(F->bnrm) && F->bnrm>0.0)
        F->nrm2 = 1.0;
    else
    {
        lis_vector_nrm2(F->v[0],&F->rnorm);
        F->nrm2 = F->rnorm*F->bnrm;
    }

    if( F->solver->recycle_valid && F->solver->recycle_dim>0 )
    {

        for(F->i=0;F->i<F->solver->recycle_dim;F->i++)
        {
            lis_vector_dot(F->solver->recycle_c[F->i],F->v[0],&F->alpha);
            lis_vector_axpy(F->alpha,F->solver->recycle_zu[F->i],F->x);
            lis_vector_axpy(-F->alpha,F->solver->recycle_c[F->i],F->v[0]);
        }
        lis_vector_nrm2(F->v[0],&F->rnorm);
        F->nrm2 = F->rnorm*F->bnrm;
        if( F->nrm2<=F->tol )
        {
            F->solver->retcode = LIS_SUCCESS;
            F->solver->iter = 0;
            F->solver->resid = F->nrm2;
            F->solver->ptime = F->ptime;
            F->err = LIS_SUCCESS;
            goto done;
        }
    }

    while( F->iter<F->maxiter )
    {
        F->kcur = F->solver->recycle_valid ? F->solver->recycle_dim : 0;
        memset(F->h,0,(size_t)(F->hld+1)*(F->hld+2)*sizeof(LIS_SCALAR));
        memset(F->hraw,0,(size_t)(F->m+1)*F->m*sizeof(LIS_SCALAR));
        memset(F->bcoef,0,(size_t)_max(1,F->kmax)*F->m*sizeof(LIS_SCALAR));

        lis_vector_nrm2(F->v[0],&F->rnorm);
        if( F->rnorm<=GCRODR_RESID_TINY || F->rnorm!=F->rnorm )
        {
            F->nrm2 = F->rnorm*F->bnrm;
            F->err = (F->nrm2<=F->tol) ? LIS_SUCCESS : LIS_BREAKDOWN;
            break;
        }
        lis_vector_scale(1.0/F->rnorm,F->v[0]);
        lis_vector_set_all(0.0,F->s);
        F->s->value[0] = F->rnorm;

        F->jdim = 0;
        F->lucky = LIS_FALSE;
        F->candidate = LIS_FALSE;
        F->cycle_steps = F->m-F->kcur;
        if( F->cycle_steps<1 ) F->cycle_steps = 1;
        for(F->j=0;F->j<F->cycle_steps && F->iter<F->maxiter;F->j++)
        {
            F->iter++;

            F->time = lis_wtime();
            F->pc=3;
            return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->v[F->j],F->z[F->j]);
            resume_3: F->pc=0; F->err=K->retval;
            F->ptime += lis_wtime()-F->time;
            if( F->err ) goto fail;

            F->pc=4;
            return linear_suspend(K,LIS_RCI_OP_A,F->A,F->z[F->j],F->v[F->j+1]);
            resume_4: F->pc=0; F->err=K->retval;
            if( F->err ) goto fail;


            for(F->pass=0;F->pass<2;F->pass++)
            {
                for(F->i=0;F->i<F->kcur;F->i++)
                {
                    lis_vector_dot(F->solver->recycle_c[F->i],F->v[F->j+1],&F->t);
                    F->bcoef[F->i+F->j*_max(1,F->kmax)] += F->t;
                    lis_vector_axpy(-F->t,F->solver->recycle_c[F->i],F->v[F->j+1]);
                }
            }


            for(F->pass=0;F->pass<2;F->pass++)
            {
                for(F->i=0;F->i<=F->j;F->i++)
                {
                    lis_vector_dot(F->v[F->i],F->v[F->j+1],&F->t);
                    F->hraw[F->i+F->j*F->hld] += F->t;
                    lis_vector_axpy(-F->t,F->v[F->i],F->v[F->j+1]);
                }
            }
            lis_vector_nrm2(F->v[F->j+1],&F->wnorm);
            F->hraw[(F->j+1)+F->j*F->hld] = F->wnorm;
            if( F->wnorm>GCRODR_TINY && F->wnorm==F->wnorm )
                lis_vector_scale(1.0/F->wnorm,F->v[F->j+1]);
            else
                F->lucky = LIS_TRUE;


            for(F->i=0;F->i<=F->j+1;F->i++) F->h[F->i+F->j*F->hld] = F->hraw[F->i+F->j*F->hld];


            for(F->k=1;F->k<=F->j;F->k++)
            {
                F->jj = F->k-1;
                F->t = F->h[F->jj+F->j*F->hld];
                F->aa =
                    F->h[F->jj+F->cs]*F->t
                    +F->h[F->jj+F->sn]*F->h[F->k+F->j*F->hld];

                F->bb =
                    -F->h[F->jj+F->sn]*F->t
                    +F->h[F->jj+F->cs]*F->h[F->k+F->j*F->hld];
                F->h[F->jj+F->j*F->hld] = F->aa;
                F->h[F->k+F->j*F->hld] = F->bb;
            }

            F->ii = F->j;
            F->i1 = F->j+1;
            F->aa = F->h[F->ii+F->j*F->hld];
            F->bb = F->h[F->i1+F->j*F->hld];
            F->a2 = F->aa*F->aa;
            F->b2 = F->bb*F->bb;
            F->rr = sqrt(F->a2+F->b2);
            if( F->rr==0.0 || F->rr!=F->rr )
            {
                if( F->lucky )
                {
                    F->jdim = F->j+1;
                    break;
                }
                F->err = LIS_BREAKDOWN;
                goto fail;
            }
            F->h[F->ii+F->cs] =
                F->aa/F->rr;

            F->h[F->ii+F->sn] =
                F->bb/F->rr;


            F->s->value[F->i1] =
                -F->h[F->ii+F->sn]
                *F->s->value[F->ii];


            F->s->value[F->ii] =
                F->h[F->ii+F->cs]
                *F->s->value[F->ii];


            F->h[F->ii+F->j*F->hld] =
                F->h[F->ii+F->cs]*F->h[F->ii+F->j*F->hld]
                +F->h[F->ii+F->sn]*F->h[F->i1+F->j*F->hld];



            F->nrm2 =
                fabs(
                    F->s->value[F->i1])
                *F->bnrm;
            if( F->output )
            {
                if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[F->iter] = F->nrm2;
                if( F->output & LIS_PRINT_OUT ) lis_print_rhistory(F->comm,F->iter,F->nrm2);
            }

            lis_solver_monitor_residual(
                F->solver,F->iter,F->nrm2,F->tol,
                LIS_SOLVER_RESIDUAL_PROJECTED);

            F->candidate = lis_solver_test_convergence(
                F->solver,F->iter,F->nrm2,F->tol,
                LIS_SOLVER_RESIDUAL_PROJECTED,
                F->nrm2<=F->tol);

            F->jdim = F->j+1;
            if( F->candidate || F->lucky ) break;

            if( F->maxiter_noimp )
            {
                if( F->min_nrm2>F->nrm2 )
                {
                    F->min_nrm2=F->nrm2;
                    F->noimp_count=0;
                }
                else if( ++F->noimp_count>F->maxiter_noimp )
                {
                    F->err = LIS_MAXITER;
                    goto fail;
                }
            }
        }

        if( F->jdim<=0 )
        {
            F->err = LIS_BREAKDOWN;
            goto fail;
        }


        F->ii = F->jdim-1;
        if( F->h[F->ii+F->ii*F->hld]==0.0 || F->h[F->ii+F->ii*F->hld]!=F->h[F->ii+F->ii*F->hld] )
        {
            F->err = LIS_BREAKDOWN;
            goto fail;
        }
        F->s->value[F->ii] /= F->h[F->ii+F->ii*F->hld];
        for(F->k=1;F->k<=F->ii;F->k++)
        {
            F->jj = F->ii-F->k;
            F->t = F->s->value[F->jj];
            for(F->j=F->jj+1;F->j<=F->ii;F->j++) F->t -= F->h[F->jj+F->j*F->hld]*F->s->value[F->j];
            if( F->h[F->jj+F->jj*F->hld]==0.0 || F->h[F->jj+F->jj*F->hld]!=F->h[F->jj+F->jj*F->hld] )
            {
                F->err = LIS_BREAKDOWN;
                goto fail;
            }
            F->s->value[F->jj] = F->t/F->h[F->jj+F->jj*F->hld];
        }




        F->nrm2 = NAN;
        for(F->j=0;F->j<F->jdim;F->j++) lis_vector_axpy(F->s->value[F->j],F->z[F->j],F->x);
        for(F->i=0;F->i<F->kcur;F->i++)
        {
            F->gamma = 0.0;
            for(F->j=0;F->j<F->jdim;F->j++)
                F->gamma += F->bcoef[F->i+F->j*_max(1,F->kmax)]*F->s->value[F->j];
            lis_vector_axpy(-F->gamma,F->solver->recycle_zu[F->i],F->x);
        }


        if( F->kmax>0 && !F->lucky )
        {
            K->update.solver = F->solver;
            K->update.v = F->v;
            K->update.jdim = F->jdim;
            K->update.hraw = F->hraw;
            K->update.hld = F->hld;
            K->update.bcoef = F->bcoef;
            K->update.bld = _max(1,F->kmax);
            K->update.pc=0; K->update.parent=5;
            F->pc=5; K->active=4; return L_STEP_AGAIN;
            resume_5: F->pc=0; F->err=K->retval;
            if( F->err ) goto fail;
        }


        F->pc=6;
        return linear_suspend(K,LIS_RCI_OP_A,F->A,F->x,F->wtmp);
        resume_6: F->pc=0; F->err=K->retval;
        if( F->err ) goto fail;
        lis_vector_xpay(F->b,-1.0,F->wtmp);
        lis_vector_copy(F->wtmp,F->v[0]);




        if( F->solver->recycle_valid )
        {
            for(F->i=0;F->i<F->solver->recycle_dim;F->i++)
            {
                lis_vector_dot(F->solver->recycle_c[F->i],F->v[0],&F->alpha);
                lis_vector_axpy(F->alpha,F->solver->recycle_zu[F->i],F->x);
                lis_vector_axpy(-F->alpha,F->solver->recycle_c[F->i],F->v[0]);
            }
        }
        lis_vector_nrm2(F->v[0],&F->rnorm);
        F->nrm2 = F->rnorm*F->bnrm;

        if( F->nrm2!=F->nrm2 || fabs(F->nrm2)>LIS_SCALAR_MAX )
        {
            F->err = LIS_BREAKDOWN;
            goto fail;
        }
        if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[F->iter] = F->nrm2;

        lis_solver_monitor_residual(
            F->solver,F->iter,F->nrm2,F->tol,
            LIS_SOLVER_RESIDUAL_TRUE);

        if( lis_solver_test_convergence(
                F->solver,F->iter,F->nrm2,F->tol,
                LIS_SOLVER_RESIDUAL_TRUE,
                F->nrm2<=F->tol) )
        {
            F->err = LIS_SUCCESS;
            break;
        }
    }

    if( F->iter>=F->maxiter ) F->err = LIS_MAXITER;
    F->solver->retcode = F->err;
    F->solver->iter = F->iter;
    F->solver->resid = F->nrm2;
    F->solver->ptime = F->ptime;
    goto done;

fail:

    lis_gcrodr_clear_recycle(F->solver);
    F->solver->retcode = F->err;
    F->solver->iter = F->iter;
    F->solver->resid = F->nrm2;
    F->solver->ptime = F->ptime;

done:
    lis_free(F->h); F->h=NULL;
    lis_free(F->hraw);
    lis_free(F->bcoef);
    LIS_DEBUG_FUNC_OUT;
    return linear_return(K,F->parent,F->err);

}
static int linear_step(lis_rci_native_gcro_restore *K){switch(K->active){
case 1:return linear_step_initial(K);
case 2:return linear_step_refresh(K);
case 3:return linear_step_assess(K);
case 4:return linear_step_update(K);
case 5:return linear_step_gcro(K);
default:return L_STEP_STOP;}}
