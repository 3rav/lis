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
struct linear_frame_cg { int pc,parent;
    LIS_SOLVER solver;
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR r;
    LIS_VECTOR z;
    LIS_VECTOR p;
    LIS_VECTOR q;
    LIS_SCALAR alpha;
    LIS_SCALAR beta;
    LIS_SCALAR rho;
    LIS_SCALAR rho_old;
    LIS_SCALAR dot_pq;
    LIS_REAL bnrm2;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT output;
    LIS_INT conv;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    LIS_INT err;
    double time;
    double ptime;
};
struct linear_frame_bi { int pc,parent;
    LIS_SOLVER solver;
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR r;
    LIS_VECTOR rtld;
    LIS_VECTOR t;
    LIS_VECTOR p;
    LIS_VECTOR v;
    LIS_VECTOR s;
    LIS_VECTOR phat;
    LIS_VECTOR shat;
    LIS_SCALAR alpha;
    LIS_SCALAR beta;
    LIS_SCALAR omega;
    LIS_SCALAR rho;
    LIS_SCALAR rho_old;
    LIS_SCALAR tmpdot1;
    LIS_SCALAR tmpdot2;
    LIS_REAL bnrm2;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_INT err;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT output;
    LIS_INT conv;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    double time;
    double ptime;
};
struct linear_frame_gm { int pc,parent;
    LIS_SOLVER solver;
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR r;
    LIS_VECTOR s;
    LIS_VECTOR z;
    LIS_VECTOR * v;
    LIS_SCALAR * h;
    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR rr;
    LIS_SCALAR a2;
    LIS_SCALAR b2;
    LIS_SCALAR t;
    LIS_REAL bnrm2;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_INT err;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT n;
    LIS_INT output;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    double time;
    double ptime;
    LIS_REAL rnorm;
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
};
struct linear_frame_fgm { int pc,parent;
    LIS_SOLVER solver;
    LIS_Comm comm;
    LIS_MATRIX A;
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR s;
    LIS_VECTOR * z;
    LIS_VECTOR * v;
    LIS_SCALAR * h;
    LIS_SCALAR aa;
    LIS_SCALAR bb;
    LIS_SCALAR rr;
    LIS_SCALAR a2;
    LIS_SCALAR b2;
    LIS_SCALAR t;
    LIS_REAL bnrm2;
    LIS_REAL nrm2;
    LIS_REAL tol;
    LIS_REAL min_nrm2;
    LIS_INT err;
    LIS_INT iter;
    LIS_INT maxiter;
    LIS_INT output;
    LIS_INT maxiter_noimp;
    LIS_INT noimp_count;
    LIS_INT candidate;
    LIS_INT noimp_stop;
    double time;
    double ptime;
    LIS_REAL rnorm;
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
};
struct lis_rci_native_linear_restore { int active,method,status,strict; LIS_INT retval,native_code; uint64_t token,matvecs,precons; LIS_SOLVER solver; lis_rci_native_linear_request request;
    struct linear_frame_initial initial;
    struct linear_frame_cg cg;
    struct linear_frame_bi bi;
    struct linear_frame_gm gm;
    struct linear_frame_fgm fgm;
};
#define L_STEP_AGAIN 1
#define L_STEP_YIELD 2
#define L_STEP_STOP 3
static int linear_suspend(lis_rci_native_linear_restore *,uint32_t,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);
static int linear_return(lis_rci_native_linear_restore *,int,LIS_INT);
#undef __FUNC__
#define __FUNC__ "lis_solver_get_initial_residual"
static int linear_step_initial(lis_rci_native_linear_restore *K) {
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
#define __FUNC__ "lis_cg"
static int linear_step_cg(lis_rci_native_linear_restore *K) {
struct linear_frame_cg *F=&K->cg;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}







	F->nrm2 = 1.0;



	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	F->A       = F->solver->A;
	F->x       = F->solver->x;
	F->maxiter = F->solver->options[LIS_OPTIONS_MAXITER];
	F->output  = F->solver->options[LIS_OPTIONS_OUTPUT];
	F->conv    = F->solver->options[LIS_OPTIONS_CONV_COND];
	F->maxiter_noimp = F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
	F->noimp_count   = 0;
	F->min_nrm2      = LIS_SCALAR_MAX;
	F->ptime   = 0.0;


	F->z       = F->solver->work[0];
	F->q       = F->solver->work[1];
	F->r       = F->solver->work[2];
	F->p       = F->solver->work[3];
	F->rho_old = (LIS_SCALAR)1.0;
	F->beta    = (LIS_SCALAR)0.0;

	K->initial.solver = F->solver;
	K->initial.M = NULL;
	K->initial.t = NULL;
	K->initial.r = F->r;
	K->initial.bnrm2 = &F->bnrm2;
	K->initial.pc=0; K->initial.parent=2;
	F->pc=1; K->active=1; return L_STEP_AGAIN;
	resume_1: F->pc=0;
	if(K->retval)
	{
		LIS_DEBUG_FUNC_OUT;
		if( F->solver->retcode != LIS_SUCCESS )
		{
			F->solver->resid = NAN;
			F->solver->ptime = F->ptime;
			return linear_return(K,F->parent,F->solver->retcode);
		}
		return linear_return(K,F->parent,LIS_SUCCESS);
	}
	F->tol     = F->solver->tol;

	lis_vector_set_all(0.0,F->p);

	for( F->iter=1; F->iter<=F->maxiter; F->iter++ )
	{

		F->time = lis_wtime();
		F->pc=2;
		return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->r,F->z);
		resume_2: F->pc=0; F->err=K->retval;
		F->ptime += lis_wtime() - F->time;
		if(F->err)goto callback_failure;


		lis_vector_dot(F->r,F->z,&F->rho);


		F->beta = F->rho / F->rho_old;


		lis_vector_xpay(F->z,F->beta,F->p);


		F->pc=3;
		return linear_suspend(K,LIS_RCI_OP_A,F->A,F->p,F->q);
		resume_3: F->pc=0; F->err=K->retval;
		if(F->err)goto callback_failure;


		lis_vector_dot(F->p,F->q,&F->dot_pq);


		if( F->dot_pq==0.0 )
		{
			F->solver->retcode   = LIS_BREAKDOWN;
			F->solver->iter      = F->iter;
			F->solver->resid     = F->nrm2;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}


		F->alpha = F->rho / F->dot_pq;


		lis_vector_axpy(F->alpha,F->p,F->x);


		lis_vector_axpy(-F->alpha,F->q,F->r);


		lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);
		if( F->output )
		{
			if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[F->iter] = F->nrm2;
			if( F->output & LIS_PRINT_OUT ) lis_print_rhistory(F->comm,F->iter,F->nrm2);
		}

		if( F->tol >= F->nrm2 )
		{
			F->solver->retcode    = LIS_SUCCESS;
			F->solver->iter       = F->iter;
			F->solver->resid      = F->nrm2;
			F->solver->ptime      = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_SUCCESS);
		}
		F->rho_old = F->rho;

		if ( F->maxiter_noimp )
		{
			if ( F->min_nrm2 > F->nrm2 )
			{
				F->min_nrm2 = F->nrm2;
				F->noimp_count = 0;
			}
			else if ( ++F->noimp_count > F->maxiter_noimp )
			{
				F->solver->retcode = LIS_MAXITER;
				F->solver->iter = F->iter;
				F->solver->resid = F->nrm2;
				F->solver->ptime = F->ptime;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_MAXITER);
			}
		}
	}

	F->solver->retcode   = LIS_MAXITER;
	F->solver->iter      = F->iter;
	F->solver->resid     = F->nrm2;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,LIS_MAXITER);

callback_failure:



	F->solver->retcode = F->err;
	F->solver->iter    = F->iter-1;
	if( F->iter==1 )
		lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);
	F->solver->resid   = F->nrm2;
	F->solver->ptime   = F->ptime;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_bicgstab"
static int linear_step_bi(lis_rci_native_linear_restore *K) {
struct linear_frame_bi *F=&K->bi;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;case 5:goto resume_5;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}










	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	F->A       = F->solver->A;
	F->x       = F->solver->x;
	F->maxiter = F->solver->options[LIS_OPTIONS_MAXITER];
	F->output  = F->solver->options[LIS_OPTIONS_OUTPUT];
	F->conv    = F->solver->options[LIS_OPTIONS_CONV_COND];
	F->maxiter_noimp = F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
	F->noimp_count   = 0;
	F->min_nrm2      = LIS_SCALAR_MAX;
	F->ptime   = 0.0;

	F->rtld    = F->solver->work[0];
	F->r       = F->solver->work[1];
	F->s       = F->solver->work[1];
	F->t       = F->solver->work[2];
	F->p       = F->solver->work[3];
	F->v       = F->solver->work[4];
	F->phat    = F->solver->work[5];
	F->shat    = F->solver->work[6];
	F->alpha   = (LIS_SCALAR)1.0;
	F->omega   = (LIS_SCALAR)1.0;
	F->rho_old = (LIS_SCALAR)1.0;

	lis_vector_set_all(0.0,F->p);
	lis_vector_set_all(0.0,F->phat);
	lis_vector_set_all(0.0,F->s);
	lis_vector_set_all(0.0,F->shat);


	K->initial.solver = F->solver;
	K->initial.M = NULL;
	K->initial.t = NULL;
	K->initial.r = F->r;
	K->initial.bnrm2 = &F->bnrm2;
	K->initial.pc=0; K->initial.parent=3;
	F->pc=1; K->active=1; return L_STEP_AGAIN;
	resume_1: F->pc=0;
	if(K->retval)
	{
		LIS_DEBUG_FUNC_OUT;
		return linear_return(K,F->parent,LIS_SUCCESS);
	}
	F->tol     = F->solver->tol;

	lis_solver_set_shadowresidual(F->solver,F->r,F->rtld);


	for( F->iter=1; F->iter<=F->maxiter; F->iter++ )
	{

		lis_vector_dot(F->rtld,F->r,&F->rho);


		if( F->rho==0.0 || F->rho!=F->rho )
		{
			lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);
			F->solver->retcode   = LIS_BREAKDOWN;
			F->solver->iter      = F->iter;
			F->solver->resid     = F->nrm2;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}

		if( F->iter==1 )
		{
			lis_vector_copy(F->r,F->p);
		}
		else
		{


			if( F->rho_old==0.0 || F->rho_old!=F->rho_old || F->omega==0.0 || F->omega!=F->omega )
			{
				lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);
				F->solver->retcode = LIS_BREAKDOWN;
				F->solver->iter    = F->iter;
				F->solver->resid   = F->nrm2;
				F->solver->ptime   = F->ptime;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_BREAKDOWN);
			}

			F->beta = (F->rho / F->rho_old) * (F->alpha / F->omega);


			lis_vector_axpy(-F->omega,F->v,F->p);
			lis_vector_xpay(F->r,F->beta,F->p);
		}


		F->time = lis_wtime();
		F->pc=2;
		return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->p,F->phat);
		resume_2: F->pc=0; F->err=K->retval;
		F->ptime += lis_wtime()-F->time;
		if(F->err)
		{
			F->solver->retcode=F->err;
			F->solver->iter=F->iter;
			F->solver->ptime=F->ptime;
			return linear_return(K,F->parent,F->err);
		}


		F->pc=3;
		return linear_suspend(K,LIS_RCI_OP_A,F->A,F->phat,F->v);
		resume_3: F->pc=0; F->err=K->retval;
		if(F->err)
		{
			F->solver->retcode=F->err;
			F->solver->iter=F->iter;
			F->solver->ptime=F->ptime;
			return linear_return(K,F->parent,F->err);
		}


		lis_vector_dot(F->rtld,F->v,&F->tmpdot1);





		if( F->tmpdot1==0.0 || F->tmpdot1!=F->tmpdot1 )
		{
			lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}

		F->alpha = F->rho / F->tmpdot1;


		if( F->alpha!=F->alpha )
		{
			lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}


		lis_vector_axpy(-F->alpha,F->v,F->r);


		lis_solver_get_residual[F->conv](F->s,F->solver,&F->nrm2);


		if( F->nrm2 <= F->tol )
		{
			if( F->output )
			{
				if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[F->iter] = F->nrm2;
				if( F->output & LIS_PRINT_OUT ) lis_print_rhistory(F->comm,F->iter,F->nrm2);
			}

			lis_vector_axpy(F->alpha,F->phat,F->x);
			F->solver->retcode    = LIS_SUCCESS;
			F->solver->iter       = F->iter;
			F->solver->resid      = F->nrm2;
			F->solver->ptime      = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_SUCCESS);
		}


		F->time = lis_wtime();
		F->pc=4;
		return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->s,F->shat);
		resume_4: F->pc=0; F->err=K->retval;
		F->ptime += lis_wtime()-F->time;
		if(F->err)
		{
			F->solver->retcode=F->err;
			F->solver->iter=F->iter;
			F->solver->ptime=F->ptime;
			return linear_return(K,F->parent,F->err);
		}


		F->pc=5;
		return linear_suspend(K,LIS_RCI_OP_A,F->A,F->shat,F->t);
		resume_5: F->pc=0; F->err=K->retval;
		if(F->err)
		{
			F->solver->retcode=F->err;
			F->solver->iter=F->iter;
			F->solver->ptime=F->ptime;
			return linear_return(K,F->parent,F->err);
		}




		lis_vector_dot(F->t,F->s,&F->tmpdot1);
		lis_vector_dot(F->t,F->t,&F->tmpdot2);

		if( F->tmpdot1==0.0 || F->tmpdot1!=F->tmpdot1 || F->tmpdot2==0.0 || F->tmpdot2!=F->tmpdot2 )
		{
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}

		F->omega   = F->tmpdot1 / F->tmpdot2;


		if( F->omega==0.0 || F->omega!=F->omega )
		{
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}


		lis_vector_axpy(F->alpha,F->phat,F->x);
		lis_vector_axpy(F->omega,F->shat,F->x);


		lis_vector_axpy(-F->omega,F->t,F->r);


		lis_solver_get_residual[F->conv](F->r,F->solver,&F->nrm2);



		if( F->output )
		{
			if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[F->iter] = F->nrm2;
			if( F->output & LIS_PRINT_OUT ) lis_print_rhistory(F->comm,F->iter,F->nrm2);
		}

		if( F->tol >= F->nrm2 )
		{
			F->solver->retcode    = LIS_SUCCESS;
			F->solver->iter       = F->iter;
			F->solver->resid      = F->nrm2;
			F->solver->ptime      = F->ptime;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_SUCCESS);
		}

		F->rho_old = F->rho;

		if ( F->maxiter_noimp )
		{
			if ( F->min_nrm2 > F->nrm2 )
			{
				F->min_nrm2 = F->nrm2;
				F->noimp_count = 0;
			}
			else if ( ++F->noimp_count > F->maxiter_noimp )
			{
				F->solver->retcode = LIS_MAXITER;
				F->solver->iter = F->iter;
				F->solver->resid = F->nrm2;
				F->solver->ptime = F->ptime;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_MAXITER);
			}
		}
	}

	F->solver->retcode   = LIS_MAXITER;
	F->solver->iter      = F->iter;
	F->solver->resid     = F->nrm2;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,LIS_MAXITER);

}
#undef __FUNC__
#define __FUNC__ "lis_gmres"
static int linear_step_gm(lis_rci_native_linear_restore *K) {
struct linear_frame_gm *F=&K->gm;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;case 5:goto resume_5;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}



















	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	F->A       = F->solver->A;
	F->b       = F->solver->b;
	F->x       = F->solver->x;
	F->n       = F->A->n;
	F->maxiter = F->solver->options[LIS_OPTIONS_MAXITER];
	F->output  = F->solver->options[LIS_OPTIONS_OUTPUT];
	F->m       = F->solver->options[LIS_OPTIONS_RESTART];
	F->maxiter_noimp = F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
	F->noimp_count   = 0;
	F->min_nrm2      = LIS_SCALAR_MAX;
	F->h_dim   = F->m+1;
	F->ptime   = 0.0;
	F->iter    = 0;
	F->nrm2    = NAN;

	F->s       = F->solver->work[0];
	F->r       = F->solver->work[1];
	F->z       = F->solver->work[2];
	F->v       = &F->solver->work[3];

	F->h       = (LIS_SCALAR *)lis_malloc( sizeof(LIS_SCALAR)*(F->h_dim+1)*(F->h_dim+2),"lis_gmres::h" );
	if(F->h==NULL)
	{
		LIS_SETERR_MEM(sizeof(LIS_SCALAR)*(F->h_dim+1)*(F->h_dim+2));
		F->solver->retcode = LIS_ERR_OUT_OF_MEMORY;
		F->solver->iter = 0;
		F->solver->resid = NAN;
		F->solver->ptime = 0.0;
		LIS_DEBUG_FUNC_OUT;
		return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
	}
	F->cs      = (F->m+1)*F->h_dim;
	F->sn      = (F->m+2)*F->h_dim;


	K->initial.solver = F->solver;
	K->initial.M = NULL;
	K->initial.t = NULL;
	K->initial.r = F->v[0];
	K->initial.bnrm2 = &F->bnrm2;
	K->initial.pc=0; K->initial.parent=4;
	F->pc=1; K->active=1; return L_STEP_AGAIN;
	resume_1: F->pc=0;
	if(K->retval)
	{
		lis_free(F->h); F->h=NULL;
		if(F->solver->retcode!=LIS_SUCCESS)
		{
			F->solver->resid = NAN;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,F->solver->retcode);
		}
		LIS_DEBUG_FUNC_OUT;
		return linear_return(K,F->parent,LIS_SUCCESS);
	}
	F->tol     = F->solver->tol;
	if(F->solver->options[LIS_OPTIONS_CONV_COND]==LIS_CONV_COND_NRM2_R &&
		isfinite(F->bnrm2) && F->bnrm2>0.0)
		F->nrm2 = 1.0;
	else
	{
		lis_vector_nrm2(F->v[0],&F->nrm2);
		F->nrm2 *= F->bnrm2;
	}


	F->iter=0;
	while( F->iter<F->maxiter )
	{


		lis_vector_nrm2(F->v[0],&F->rnorm);
		lis_vector_scale(1.0/F->rnorm,F->v[0]);


		lis_vector_set_all(0,F->s);
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
			return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->v[F->iiv],F->z);
			resume_2: F->pc=0; F->err=K->retval;
			F->ptime += lis_wtime()-F->time;
			if(F->err) goto callback_failure;


			F->pc=3;
			return linear_suspend(K,LIS_RCI_OP_A,F->A,F->z,F->v[F->i1v]);
			resume_3: F->pc=0; F->err=K->retval;
			if(F->err) goto callback_failure;

			for(F->k=0;F->k<F->i;F->k++)
			{


				lis_vector_dot(F->v[F->i1v],F->v[F->k],&F->t);
				F->h[F->k+F->iih] = F->t;
				lis_vector_axpy(-F->t,F->v[F->k],F->v[F->i1v]);
			}


			lis_vector_nrm2(F->v[F->i1v],(LIS_REAL *)&F->t);
			F->h[F->i1+F->iih] = F->t;

			if( F->t!=0.0 && F->t==F->t )
			{
				lis_vector_scale(1.0/F->t,F->v[F->i1v]);
			}

			for(F->k=1;F->k<=F->ii;F->k++)
			{
				F->jj  =  F->k-1;
				F->t   =  F->h[F->jj+F->iih];
				F->aa  =  F->h[F->jj+F->cs]*F->t;
				F->aa +=  F->h[F->jj+F->sn]*F->h[F->k+F->iih];
				F->bb  = -F->h[F->jj+F->sn]*F->t;
				F->bb +=  F->h[F->jj+F->cs]*F->h[F->k+F->iih];
				F->h[F->jj+F->iih] = F->aa;
				F->h[F->k+F->iih] = F->bb;
			}
			F->aa = F->h[F->ii+F->iih];
			F->bb = F->h[F->i1+F->iih];
			F->a2 = F->aa*F->aa;
			F->b2 = F->bb*F->bb;
			F->rr = sqrt(F->a2+F->b2);

			if( F->rr==0.0 || F->rr!=F->rr )
			{
				F->nrm2 = fabs(F->s->value[F->ii])*F->bnrm2;
				F->solver->retcode = LIS_BREAKDOWN;
				F->solver->iter    = F->iter;
				F->solver->resid   = F->nrm2;
				F->solver->ptime   = F->ptime;
				lis_free(F->h); F->h=NULL;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_BREAKDOWN);
			}
			F->h[F->ii+F->cs] = F->aa/F->rr;
			F->h[F->ii+F->sn] = F->bb/F->rr;
			F->s->value[F->i1] = -F->h[F->ii+F->sn]*F->s->value[F->ii];
			F->s->value[F->ii] =  F->h[F->ii+F->cs]*F->s->value[F->ii];

			F->aa  =  F->h[F->ii+F->cs]*F->h[F->ii+F->iih];
			F->aa +=  F->h[F->ii+F->sn]*F->h[F->i1+F->iih];
			F->h[F->ii+F->iih] = F->aa;


			F->nrm2 = fabs(F->s->value[F->i1])*F->bnrm2;

			if( F->output )
			{
				if( F->output & LIS_PRINT_MEM ) F->solver->rhistory[F->iter] = F->nrm2;
				if( F->output & LIS_PRINT_OUT ) lis_print_rhistory(F->comm,F->iter,F->nrm2);
			}

			lis_solver_monitor_residual(
				F->solver,F->iter,F->nrm2,F->tol,
				LIS_SOLVER_RESIDUAL_PROJECTED);
			if( lis_solver_test_convergence(
				F->solver,F->iter,F->nrm2,F->tol,
				LIS_SOLVER_RESIDUAL_PROJECTED,
				F->tol >= F->nrm2) ) break;

			if (F->maxiter_noimp)
			{
				if (F->min_nrm2 > F->nrm2)
				{
					F->min_nrm2 = F->nrm2;
					F->noimp_count = 0;
				}
				else if (++F->noimp_count > F->maxiter_noimp)
				{
					F->solver->retcode = LIS_MAXITER;
					F->solver->iter = F->iter;
					F->solver->resid = F->nrm2;
					F->solver->ptime = F->ptime;
					lis_free(F->h); F->h=NULL;
					LIS_DEBUG_FUNC_OUT;
					return linear_return(K,F->parent,LIS_MAXITER);
				}
			}

		} while( F->i<F->m && F->iter <F->maxiter );



		if( F->h[F->ii+F->iih]==0.0 || F->h[F->ii+F->iih]!=F->h[F->ii+F->iih] )
		{
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			lis_free(F->h); F->h=NULL;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}

		F->s->value[F->ii] = F->s->value[F->ii]/F->h[F->ii+F->iih];
		for(F->k=1;F->k<=F->ii;F->k++)
		{
			F->jj = F->ii-F->k;
			F->t  = F->s->value[F->jj];
			for(F->j=F->jj+1;F->j<=F->ii;F->j++)
			{
				F->t -= F->h[F->jj+F->j*F->h_dim]*F->s->value[F->j];
			}
			if( F->h[F->jj+F->jj*F->h_dim]==0.0 || F->h[F->jj+F->jj*F->h_dim]!=F->h[F->jj+F->jj*F->h_dim] )
			{
				F->solver->retcode = LIS_BREAKDOWN;
				F->solver->iter    = F->iter;
				F->solver->resid   = F->nrm2;
				F->solver->ptime   = F->ptime;
				lis_free(F->h); F->h=NULL;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_BREAKDOWN);
			}

			F->s->value[F->jj] = F->t/F->h[F->jj+F->jj*F->h_dim];
		}

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(LIS_INT local_k=0;local_k<F->n;local_k++)
		{
			F->z->value[local_k] = F->s->value[0]*F->v[0]->value[local_k];
		}
		for(F->j=1;F->j<=F->ii;F->j++)
		{
			lis_vector_axpy(F->s->value[F->j],F->v[F->j],F->z);
		}


		F->time = lis_wtime();
		F->pc=4;
		return linear_suspend(K,LIS_RCI_OP_PRECON,F->solver->A,F->z,F->r);
		resume_4: F->pc=0; F->err=K->retval;
		F->ptime += lis_wtime()-F->time;
		if(F->err) goto callback_failure;


		lis_vector_axpy(1,F->r,F->x);






		F->nrm2 = NAN;
		F->pc=5;
		return linear_suspend(K,LIS_RCI_OP_A,F->A,F->x,F->z);
		resume_5: F->pc=0; F->err=K->retval;
		if(F->err) goto callback_failure;
		lis_vector_xpay(F->b,-1.0,F->z);
		lis_vector_copy(F->z,F->v[0]);
		lis_vector_nrm2(F->v[0],&F->rnorm);
		F->nrm2 = F->rnorm * F->bnrm2;

		if( F->nrm2!=F->nrm2 || fabs(F->nrm2)>LIS_SCALAR_MAX )
		{
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			lis_free(F->h); F->h=NULL;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}

		if( F->output & LIS_PRINT_MEM )
		{
			F->solver->rhistory[F->iter] = F->nrm2;
		}

		lis_solver_monitor_residual(
			F->solver,F->iter,F->nrm2,F->tol,
			LIS_SOLVER_RESIDUAL_TRUE);
		if( lis_solver_test_convergence(
			F->solver,F->iter,F->nrm2,F->tol,
			LIS_SOLVER_RESIDUAL_TRUE,
			F->tol >= F->nrm2) )
		{
			F->solver->retcode = LIS_SUCCESS;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			lis_free(F->h); F->h=NULL;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_SUCCESS);
		}
	}

	F->solver->retcode   = LIS_MAXITER;
	F->solver->iter      = F->iter+1;
	F->solver->resid     = F->nrm2;
	lis_free(F->h); F->h=NULL;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,LIS_MAXITER);

callback_failure:
	F->solver->retcode = F->err;
	F->solver->iter = F->iter;
	F->solver->ptime = F->ptime;
	F->solver->resid = F->nrm2;
	lis_free(F->h); F->h=NULL;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_fgmres"
static int linear_step_fgm(lis_rci_native_linear_restore *K) {
struct linear_frame_fgm *F=&K->fgm;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;default:return linear_return(K,F->parent,LIS_BREAKDOWN);}




















	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	F->A       = F->solver->A;
	F->b       = F->solver->b;
	F->x       = F->solver->x;
	F->maxiter = F->solver->options[LIS_OPTIONS_MAXITER];
	F->output  = F->solver->options[LIS_OPTIONS_OUTPUT];
	F->m       = F->solver->options[LIS_OPTIONS_RESTART];
	F->maxiter_noimp = F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
	F->noimp_count   = 0;
	F->min_nrm2      = LIS_SCALAR_MAX;
     F->noimp_stop    = LIS_FALSE;
	F->h_dim   = F->m+1;
	F->ptime   = 0.0;
	F->iter    = 0;
	F->nrm2    = NAN;

	F->s       = F->solver->work[0];
	F->z       = &F->solver->work[2];
	F->v       = &F->solver->work[F->m+2];

	F->h       = (LIS_SCALAR *)lis_malloc( sizeof(LIS_SCALAR)*(F->h_dim+1)*(F->h_dim+2),"lis_gmres::h" );
	if(F->h==NULL)
	{
		LIS_SETERR_MEM(sizeof(LIS_SCALAR)*(F->h_dim+1)*(F->h_dim+2));
		F->solver->retcode = LIS_ERR_OUT_OF_MEMORY;
		F->solver->iter = 0;
		F->solver->resid = NAN;
		F->solver->ptime = 0.0;
		LIS_DEBUG_FUNC_OUT;
		return linear_return(K,F->parent,LIS_ERR_OUT_OF_MEMORY);
	}
	F->cs      = (F->m+1)*F->h_dim;
	F->sn      = (F->m+2)*F->h_dim;


	K->initial.solver = F->solver;
	K->initial.M = NULL;
	K->initial.t = NULL;
	K->initial.r = F->v[0];
	K->initial.bnrm2 = &F->bnrm2;
	K->initial.pc=0; K->initial.parent=5;
	F->pc=1; K->active=1; return L_STEP_AGAIN;
	resume_1: F->pc=0;
	if(K->retval)
	{
		lis_free(F->h); F->h=NULL;
		if(F->solver->retcode!=LIS_SUCCESS)
		{
			F->solver->resid = NAN;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,F->solver->retcode);
		}
		LIS_DEBUG_FUNC_OUT;
		return linear_return(K,F->parent,LIS_SUCCESS);
	}
	F->tol     = F->solver->tol;
	if(F->solver->options[LIS_OPTIONS_CONV_COND]==LIS_CONV_COND_NRM2_R &&
		isfinite(F->bnrm2) && F->bnrm2>0.0)
		F->nrm2 = 1.0;
	else
	{
		lis_vector_nrm2(F->v[0],&F->nrm2);
		F->nrm2 *= F->bnrm2;
	}
	F->rnorm   = 1.0/F->bnrm2;






	F->iter=0;
	while( F->iter<F->maxiter )
	{
		F->candidate = LIS_FALSE;


		lis_vector_scale(F->bnrm2,F->v[0]);


		lis_vector_set_all(0,F->s);
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
			F->ptime += lis_wtime()-F->time;

			if(F->err) goto callback_failure;


			F->pc=3;
			return linear_suspend(K,LIS_RCI_OP_A,F->A,F->z[F->iiv],F->v[F->i1v]);
			resume_3: F->pc=0; F->err=K->retval;
			if(F->err) goto callback_failure;

			for(F->k=0;F->k<F->i;F->k++)
			{


				lis_vector_dot(F->v[F->i1v],F->v[F->k],&F->t);
				F->h[F->k+F->iih] = F->t;
				lis_vector_axpy(-F->t,F->v[F->k],F->v[F->i1v]);
			}


			lis_vector_nrm2(F->v[F->i1v],(LIS_REAL *)&F->t);
			F->h[F->i1+F->iih] = F->t;

			if( F->t!=0.0 && F->t==F->t )
			{
				lis_vector_scale(1.0/F->t,F->v[F->i1v]);
			}

			for(F->k=1;F->k<=F->ii;F->k++)
			{
				F->jj  = F->k-1;
				F->t   =  F->h[F->jj+F->iih];
				F->aa  =  F->h[F->jj+F->cs]*F->t;
				F->aa +=  F->h[F->jj+F->sn]*F->h[F->k+F->iih];
				F->bb  = -F->h[F->jj+F->sn]*F->t;
				F->bb +=  F->h[F->jj+F->cs]*F->h[F->k+F->iih];
				F->h[F->jj+F->iih] = F->aa;
				F->h[F->k+F->iih] = F->bb;
			}
			F->aa = F->h[F->ii+F->iih];
			F->bb = F->h[F->i1+F->iih];
			F->a2 = F->aa*F->aa;
			F->b2 = F->bb*F->bb;
			F->rr = sqrt(F->a2+F->b2);

			if( F->rr==0.0 || F->rr!=F->rr )
			{
				F->nrm2 = fabs(F->s->value[F->ii]);
				F->solver->retcode = LIS_BREAKDOWN;
				F->solver->iter    = F->iter;
				F->solver->resid   = F->nrm2;
				F->solver->ptime   = F->ptime;
				lis_free(F->h); F->h=NULL;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_BREAKDOWN);
			}
			F->h[F->ii+F->cs] = F->aa/F->rr;
			F->h[F->ii+F->sn] = F->bb/F->rr;
			F->s->value[F->i1] = -F->h[F->ii+F->sn]*F->s->value[F->ii];
			F->s->value[F->ii] =  F->h[F->ii+F->cs]*F->s->value[F->ii];

			F->aa  =  F->h[F->ii+F->cs]*F->h[F->ii+F->iih];
			F->aa +=  F->h[F->ii+F->sn]*F->h[F->i1+F->iih];
			F->h[F->ii+F->iih] = F->aa;


			F->nrm2 = fabs(F->s->value[F->i1]) * F->solver->bnrm;

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
				F->tol >= F->nrm2);
			if( F->candidate ) break;

			if (F->maxiter_noimp)
			{
				if (F->min_nrm2 > F->nrm2)
				{
					F->min_nrm2 = F->nrm2;
					F->noimp_count = 0;
				}
				else if (++F->noimp_count > F->maxiter_noimp)
				{





					F->noimp_stop = LIS_TRUE;
					break;
				}
			}

		} while( F->i<F->m && F->iter <F->maxiter );



		if( F->h[F->ii+F->iih]==0.0 || F->h[F->ii+F->iih]!=F->h[F->ii+F->iih] )
		{
			F->solver->retcode = LIS_BREAKDOWN;
			F->solver->iter    = F->iter;
			F->solver->resid   = F->nrm2;
			F->solver->ptime   = F->ptime;
			lis_free(F->h); F->h=NULL;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_BREAKDOWN);
		}

		F->s->value[F->ii] = F->s->value[F->ii]/F->h[F->ii+F->iih];
		for(F->k=1;F->k<=F->ii;F->k++)
		{
			F->jj = F->ii-F->k;
			F->t  = F->s->value[F->jj];
			for(F->j=F->jj+1;F->j<=F->ii;F->j++)
			{
				F->t -= F->h[F->jj+F->j*F->h_dim]*F->s->value[F->j];
			}
			if( F->h[F->jj+F->jj*F->h_dim]==0.0 || F->h[F->jj+F->jj*F->h_dim]!=F->h[F->jj+F->jj*F->h_dim] )
			{
				F->solver->retcode = LIS_BREAKDOWN;
				F->solver->iter    = F->iter;
				F->solver->resid   = F->nrm2;
				F->solver->ptime   = F->ptime;
				lis_free(F->h); F->h=NULL;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_BREAKDOWN);
			}

			F->s->value[F->jj] = F->t/F->h[F->jj+F->jj*F->h_dim];
		}

		for(F->j=0;F->j<=F->ii;F->j++)
		{
			lis_vector_axpy(F->s->value[F->j],F->z[F->j],F->x);
		}





		if( F->noimp_stop )
		{
			F->solver->retcode = LIS_MAXITER;
			F->solver->iter = F->iter;
			F->solver->resid = F->nrm2;
			F->solver->ptime = F->ptime;
			lis_free(F->h); F->h=NULL;
			LIS_DEBUG_FUNC_OUT;
			return linear_return(K,F->parent,LIS_MAXITER);
		}






		F->nrm2 = NAN;
		F->pc=4;
		return linear_suspend(K,LIS_RCI_OP_A,F->A,F->x,F->v[0]);
		resume_4: F->pc=0; F->err=K->retval;
		if(F->err) goto callback_failure;
		lis_vector_xpay(F->b,-1.0,F->v[0]);
		lis_vector_nrm2(F->v[0],&F->rnorm);

		if( F->candidate )
		{
			lis_solver_get_residual[0](F->v[0],F->solver,&F->nrm2);

			if( F->nrm2!=F->nrm2 || fabs(F->nrm2)>LIS_SCALAR_MAX )
			{
				F->solver->retcode = LIS_BREAKDOWN;
				F->solver->iter    = F->iter;
				F->solver->resid   = F->nrm2;
				F->solver->ptime   = F->ptime;
				lis_free(F->h); F->h=NULL;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_BREAKDOWN);
			}

			lis_solver_monitor_residual(
				F->solver,F->iter,F->nrm2,F->tol,
				LIS_SOLVER_RESIDUAL_TRUE);

			if( lis_solver_test_convergence(
				F->solver,F->iter,F->nrm2,F->tol,
				LIS_SOLVER_RESIDUAL_TRUE,
				F->tol >= F->nrm2) )
			{
				F->solver->retcode    = LIS_SUCCESS;
				F->solver->iter       = F->iter;
				F->solver->resid      = F->nrm2;
				F->solver->ptime      = F->ptime;
				lis_free(F->h); F->h=NULL;
				LIS_DEBUG_FUNC_OUT;
				return linear_return(K,F->parent,LIS_SUCCESS);
			}
		}

		F->bnrm2 = 1.0/F->rnorm;
	}

	F->solver->retcode   = LIS_MAXITER;
	F->solver->iter      = F->iter+1;
	F->solver->resid     = F->nrm2;
	lis_free(F->h); F->h=NULL;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,LIS_MAXITER);

callback_failure:
	F->solver->retcode = F->err;
	F->solver->iter = F->iter;
	F->solver->ptime = F->ptime;
	F->solver->resid = F->nrm2;
	lis_free(F->h); F->h=NULL;
	LIS_DEBUG_FUNC_OUT;
	return linear_return(K,F->parent,F->err);

}
static int linear_step(lis_rci_native_linear_restore *K){switch(K->active){
case 1:return linear_step_initial(K);
case 2:return linear_step_cg(K);
case 3:return linear_step_bi(K);
case 4:return linear_step_gm(K);
case 5:return linear_step_fgm(K);
default:return L_STEP_STOP;}}
