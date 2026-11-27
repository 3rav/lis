/* Generated continuation; native LIS source copyright and license apply.
 * Sources: src/esolver/lis_esolver.c, src/esolver/lis_esolver_block.c, src/esolver/lis_esolver_feast.c, src/esolver/lis_krylov.c, src/solver/lis_solver.c, src/solver/lis_solver_bicgstab.c, src/solver/lis_solver_cg.c, src/solver/lis_solver_fgcrodr.c, src/solver/lis_solver_gcrodr.c, src/solver/lis_solver_gmres.c
 */
/* Generated from pinned original LIS functions, qualification only. */
struct full_frame_outer { int pc;
 LIS_MATRIX A;
 LIS_VECTOR b;
 LIS_VECTOR x;
 LIS_SOLVER solver;
 LIS_INT err;
 LIS_INT saved_num_threads;
 LIS_INT thread_scope_changed;
 LIS_PRECON precon;
 LIS_VECTOR effective_b;
 LIS_VECTOR projected_b;
};
struct full_frame_kernel { int pc;
 LIS_MATRIX A;
 LIS_VECTOR b;
 LIS_VECTOR x;
 LIS_SOLVER solver;
 LIS_PRECON precon;
 LIS_Comm comm;
 LIS_INT nsolver;
 LIS_INT precon_type;
 LIS_INT maxiter;
 LIS_INT maxiter_noimp;
 LIS_INT err;
 LIS_REAL * rhistory;
 LIS_VECTOR xx;
 LIS_INT output;
 LIS_INT scale;
 LIS_INT conv_cond;
 LIS_INT precision;
 LIS_INT is_use_at;
 LIS_INT storage;
 LIS_INT block;
 LIS_INT i;
 LIS_INT n;
 LIS_INT shell_scale;
 double p_c_time;
 double p_i_time;
 double itime;
 LIS_REAL nrm2;
 LIS_REAL tol;
 LIS_REAL tol_w;
 LIS_VECTOR t;
 LIS_INT residual_err;
 LIS_VECTOR bb;
 LIS_VECTOR shell_b;
 LIS_MATRIX AA;
 LIS_MATRIX B;
 LIS_MATRIX shell_A;
 LIS_MATRIX Ah;
 char buf[64];
 LIS_VECTOR r;
 LIS_VECTOR z;
};
struct full_frames_linear {
 struct full_frame_outer outer;
 struct full_frame_kernel kernel;
};
static int full_linear_parent(lis_rci_native_full *S,LIS_INT rc){S->completion=rc;S->stage=1;return FULL_AGAIN;}
#undef __FUNC__
#define __FUNC__ "lis_solve"
static int full_step_outer(lis_rci_native_full *S){
struct full_frames_linear *C=S->frame; struct full_frame_outer *F=&C->outer;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return full_native_finish(S,LIS_BREAKDOWN);}





F->effective_b = F->b;
F->projected_b = NULL;

	LIS_DEBUG_FUNC_IN;


	F->err = lis_solver_thread_scope_begin(
	    F->solver,&F->saved_num_threads,&F->thread_scope_changed);
	if( F->err ) return full_native_finish(S,F->err);
	F->solver->A = F->A;	                           
	if( F->solver->options[LIS_OPTIONS_PRECON] < 0 ||
	    F->solver->options[LIS_OPTIONS_PRECON] == LIS_PRECON_TYPE_ADDS ||
	    F->solver->options[LIS_OPTIONS_PRECON] >= precon_register_type )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Invalid LIS_OPTIONS_PRECON value %D\n",F->solver->options[LIS_OPTIONS_PRECON]);
		if( F->projected_b )
		{
		        lis_vector_destroy(F->projected_b);
		        F->projected_b = NULL;
		}
		lis_solver_thread_scope_end(F->saved_num_threads,F->thread_scope_changed);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}
	F->err = lis_solver_check_user_matrix(F->A,F->solver,NULL);
	if( F->err )
	{
		if( F->projected_b )
		{
		        lis_vector_destroy(F->projected_b);
		        F->projected_b = NULL;
		}
		lis_solver_thread_scope_end(
			F->saved_num_threads,F->thread_scope_changed);
		return full_native_finish(S,F->err);
	}

	F->err = lis_solver_prepare_rhs_for_solve(
                F->A,
                F->b,
                F->solver,
                &F->effective_b,
                &F->projected_b);

        if( F->err )
        {
                F->solver->retcode = F->err;
                if( F->projected_b )
                {
                        lis_vector_destroy(F->projected_b);
                        F->projected_b = NULL;
                }
                if( F->projected_b )
                {
                        lis_vector_destroy(F->projected_b);
                        F->projected_b = NULL;
                }
                lis_solver_thread_scope_end(
                        F->saved_num_threads,
                        F->thread_scope_changed);
                return full_native_finish(S,F->err);
        }

        F->err = lis_precon_create(F->solver, &F->precon);
	if( F->err )
	{
		lis_solver_work_destroy(F->solver);
		F->solver->retcode = F->err;
		if( F->projected_b )
		{
		        lis_vector_destroy(F->projected_b);
		        F->projected_b = NULL;
		}
		lis_solver_thread_scope_end(F->saved_num_threads,F->thread_scope_changed);
		return full_native_finish(S,F->err);
	}

	                                
	C->kernel.A=F->A; C->kernel.b=F->effective_b; C->kernel.x=F->x; C->kernel.solver=F->solver; C->kernel.precon=F->precon;
F->pc=1; S->stage=2; return FULL_AGAIN;
resume_1: F->pc=0; F->err=S->completion;
	if( F->err )
	{
		lis_solver_work_destroy(F->solver);	  
		if( F->solver->precon==F->precon ) F->solver->precon = NULL;
		lis_precon_destroy(F->precon);
		F->solver->retcode = F->err;
		if( F->projected_b )
		{
		        lis_vector_destroy(F->projected_b);
		        F->projected_b = NULL;
		}
		lis_solver_thread_scope_end(F->saved_num_threads,F->thread_scope_changed);
		return full_native_finish(S,F->err);
	}

	if( F->solver->precon==F->precon ) F->solver->precon = NULL;
	lis_precon_destroy(F->precon);

	LIS_DEBUG_FUNC_OUT;
	if( F->projected_b )
	{
	        lis_vector_destroy(F->projected_b);
	        F->projected_b = NULL;
	}
	lis_solver_thread_scope_end(F->saved_num_threads,F->thread_scope_changed);
	return full_native_finish(S,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_solve_kernel"
static int full_step_kernel(lis_rci_native_full *S){
struct full_frames_linear *C=S->frame; struct full_frame_kernel *F=&C->kernel;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return full_native_finish(S,LIS_BREAKDOWN);}

  





















	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	F->nsolver     = F->solver->options[LIS_OPTIONS_SOLVER];
	F->precon_type = F->solver->options[LIS_OPTIONS_PRECON];
	F->maxiter     = F->solver->options[LIS_OPTIONS_MAXITER];
	F->output      = F->solver->options[LIS_OPTIONS_OUTPUT];
	F->scale       = F->solver->options[LIS_OPTIONS_SCALE];
	F->precision   = F->solver->options[LIS_OPTIONS_PRECISION];
	F->is_use_at   = F->solver->options[LIS_OPTIONS_USE_AT];
	F->storage     = F->solver->options[LIS_OPTIONS_STORAGE];
	F->block       = F->solver->options[LIS_OPTIONS_STORAGE_BLOCK];
	F->conv_cond   = F->solver->options[LIS_OPTIONS_CONV_COND];
	F->maxiter_noimp = F->solver->options[LIS_OPTIONS_MAXITER_NO_IMP];
	F->tol         = F->solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	F->tol_w       = F->solver->params[LIS_PARAMS_RESID_WEIGHT-LIS_OPTIONS_LEN];
	F->solver->precision = F->precision;

	F->shell_scale =
	    (F->A->matrix_type==LIS_MATRIX_USER ||
	     F->A->matrix_type==LIS_MATRIX_OPERATOR) &&
	    F->scale==LIS_SCALE_SYMM_DIAG;

	F->shell_A = NULL;
	F->shell_b = NULL;

	if( F->nsolver < 1 || F->nsolver > LIS_SOLVERS_LEN )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SOLVER is %D (Set between 1 to %D)\n",F->nsolver, LIS_SOLVERS_LEN);
		return full_linear_parent(S,LIS_ERR_ILL_ARG);
	}
	if( F->precon_type < 0 ||
	    F->precon_type == LIS_PRECON_TYPE_ADDS ||
	    F->precon_type >= precon_register_type )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Invalid LIS_OPTIONS_PRECON value %D\n",F->precon_type);
		return full_linear_parent(S,LIS_ERR_ILL_ARG);
	}
	F->err = lis_solver_check_user_matrix(F->A,F->solver,F->precon);
	if( F->err ) return full_linear_parent(S,F->err);
	if( F->maxiter<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_MAXITER(=%D) is less than 0\n",F->maxiter);
		return full_linear_parent(S,LIS_ERR_ILL_ARG);
	}
	if (F->maxiter_noimp < 0)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG, "Parameter LIS_OPTIONS_MAXITER_NO_IMP(=%D) is less than 0\n", F->maxiter_noimp);
		return full_linear_parent(S,LIS_ERR_ILL_ARG);
	}
	if( F->conv_cond>0 && lis_solver_execute_conv_cond[F->nsolver]==NULL )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Option conv_cond is not implemented for solver %s\n",lis_solvername[F->nsolver]);
		return full_linear_parent(S,LIS_ERR_ILL_ARG);
	}
		if( F->precision==LIS_PRECISION_QUAD )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"Quad precision is not enabled\n");
			return full_linear_parent(S,LIS_ERR_ILL_ARG);
		}

	F->err = lis_solver_check_params[F->nsolver](F->solver);
	if( F->err )
	{
		F->solver->retcode = F->err;
		return full_linear_parent(S,F->err);
	}
	                         

     
                                                               
     
                                                              
                                                                  
                                       
      
   lis_solver_near_nullspace_coarse_destroy(F->solver);


      
                                       
      
                                                            
                
      
                        
      
                                                               
                                                                   
      
                                                           
                                                                 
                                           
      
                                                           
                                           
       
      
                                                                    
                                                                  
                                                                    
       
    {
        LIS_INT null_issues=0;
        F->err=lis_solver_check_nullspace_configuration(F->solver,&null_issues);
        if(F->err)
        {
            F->solver->retcode=F->err;
            return full_linear_parent(S,F->err);
        }
        if(null_issues & LIS_NULLSPACE_CONFIG_EXACT_SOLVER)
        {
            LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                "exact nullspace deflation currently requires GMRES, FGMRES, GCRO-DR or FGCRO-DR\n");
            F->solver->retcode=LIS_ERR_NOT_IMPLEMENTED;
            return full_linear_parent(S,LIS_ERR_NOT_IMPLEMENTED);
        }
        if(null_issues & LIS_NULLSPACE_CONFIG_EXACT_SCALE)
        {
            LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                "exact nullspace deflation currently requires -scale none\n");
            F->solver->retcode=LIS_ERR_NOT_IMPLEMENTED;
            return full_linear_parent(S,LIS_ERR_NOT_IMPLEMENTED);
        }
        if(null_issues & LIS_NULLSPACE_CONFIG_EXACT_PRECISION)
        {
            LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                "exact nullspace deflation currently requires default precision\n");
            F->solver->retcode=LIS_ERR_NOT_IMPLEMENTED;
            return full_linear_parent(S,LIS_ERR_NOT_IMPLEMENTED);
        }
        if(null_issues & LIS_NULLSPACE_CONFIG_NEAR_COMBINATION)
        {
            LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                "near-nullspace two-level correction currently requires FGMRES, FGCRO-DR or MFAMG\n");
            F->solver->retcode=LIS_ERR_NOT_IMPLEMENTED;
            return full_linear_parent(S,LIS_ERR_NOT_IMPLEMENTED);
        }
        if(null_issues & LIS_NULLSPACE_CONFIG_NEAR_PRECISION)
        {
            LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
                "near-nullspace two-level correction currently requires double precision\n");
            F->solver->retcode=LIS_ERR_NOT_IMPLEMENTED;
            return full_linear_parent(S,LIS_ERR_NOT_IMPLEMENTED);
        }
    }

    if(F->solver->nullspace_dim>0)
    {
        if(F->solver->nullspace==NULL ||
           F->solver->nullspace[0]==NULL ||
           !lis_solver_vector_compatible(F->solver->nullspace[0],F->x))
        {
            LIS_SETERR(LIS_ERR_ILL_ARG,
                "exact nullspace basis is incompatible with the solution vector\n");
            F->solver->retcode=LIS_ERR_ILL_ARG;
            return full_linear_parent(S,LIS_ERR_ILL_ARG);
        }
    }

	F->solver->A        = F->A;
	F->solver->b        = F->b;

	                           
		F->err = lis_vector_duplicate(F->A,&F->xx);
	if( F->err )
	{
		F->solver->retcode = F->err;
		return full_linear_parent(S,F->err);
	}
	if( F->solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
	  if( F->output ) lis_printf(F->comm,"initial vector x      : all components set to 0\n");
	                lis_vector_set_all(0.0,F->xx);
	}
	else
	{
	  if( F->output ) lis_printf(F->comm,"initial vector x      : user defined\n"); 
			lis_vector_copy(F->x,F->xx);
	}

	                                                      
	if( F->solver->rhistory )
	{
		lis_free(F->solver->rhistory);
		F->solver->rhistory = NULL;
	}


	F->n       = F->A->n;
	F->t       = NULL;
	F->Ah      = NULL;


	F->p_c_time = lis_wtime();

	if( F->shell_scale )
	{
	          
                                                          
                                                             
            
	        if( F->solver->d==NULL ||
	            !lis_solver_work_vector_compatible(
	                F->solver->d,F->A,LIS_PRECISION_DEFAULT) )
	        {
	                if( F->solver->d!=NULL )
	                {
	                        lis_vector_destroy(F->solver->d);
	                        F->solver->d = NULL;
	                }

	                F->err = lis_vector_duplicate(F->A,&F->solver->d);
	                if( F->err )
	                {
	                        lis_vector_destroy(F->xx);
	                        F->solver->retcode = F->err;
	                        return full_linear_parent(S,F->err);
	                }
	        }

	        F->err = lis_matrix_get_diagonal(F->A,F->solver->d);
	        if( F->err )
	        {
	                lis_vector_destroy(F->xx);
	                F->solver->retcode = F->err;
	                return full_linear_parent(S,F->err);
	        }

		  
                                                   
                                                     
     
		F->err = lis_matrix_prepare_symm_diag_scaling(F->solver->d);
		if( F->err )
		{
			lis_vector_destroy(F->xx);
			F->solver->retcode = F->err;
			return full_linear_parent(S,F->err);
		}

	        F->err = lis_vector_duplicate(F->A,&F->shell_b);
	        if( F->err )
	        {
	                lis_vector_destroy(F->xx);
	                F->solver->retcode = F->err;
	                return full_linear_parent(S,F->err);
	        }

	          
                         
            
	        #ifdef _OPENMP
	        #pragma omp parallel for
	        #endif
	        for(LIS_INT full_i=0;full_i<F->n;full_i++){
	                F->shell_b->value[full_i] =
	                    F->solver->d->value[full_i] *
	                    F->b->value[full_i];
	        }

	          
                                             
           
                               
            
	        F->err = lis_matrix_create_scaled_operator(
	                F->A,F->solver->d,&F->shell_A);
	        if( F->err )
	        {
	                lis_vector_destroy(F->shell_b);
	                F->shell_b = NULL;

	                lis_vector_destroy(F->xx);
	                F->solver->retcode = F->err;
	                return full_linear_parent(S,F->err);
	        }
	          
                                                       
                                                       
           
                             
            
	        if( !F->solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	        {
	                #ifdef _OPENMP
	                #pragma omp parallel for
	                #endif
	                for(LIS_INT full_i=0;full_i<F->n;full_i++){
	                        F->xx->value[full_i] =
	                            F->xx->value[full_i] /
	                            F->solver->d->value[full_i];
	                }
	        }
	}
	else if( F->precon_type==LIS_PRECON_TYPE_IS )
	{
		if( F->solver->d==NULL )
		{
			F->err = lis_vector_duplicate(F->A,&F->solver->d);
			if( F->err )
			{
				return full_linear_parent(S,F->err);
			}
		}
		if( !F->A->is_scaled )
		{
			lis_matrix_scale(F->A,F->b,F->solver->d,LIS_SCALE_JACOBI);
		}
		else if( !F->b->is_scaled )
		{
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(LIS_INT full_i=0;full_i<F->n;full_i++){
				F->b->value[full_i] = F->b->value[full_i]*F->solver->d->value[full_i];
			}
		}
		if( F->nsolver >= LIS_SOLVER_JACOBI && F->nsolver <= LIS_SOLVER_SOR )
		{
			F->solver->options[LIS_OPTIONS_ISLEVEL] = 0;
		}
	}
	else if( F->nsolver >= LIS_SOLVER_JACOBI && F->nsolver <= LIS_SOLVER_SOR && F->precon_type!=LIS_PRECON_TYPE_NONE )
	{
		if( F->solver->d==NULL )
		{
			F->err = lis_vector_duplicate(F->A,&F->solver->d);
			if( F->err )
			{
				return full_linear_parent(S,F->err);
			}
		}
		if( !F->A->is_scaled )
		{
			lis_matrix_scale(F->A,F->b,F->solver->d,LIS_SCALE_JACOBI);
		}
	}
	else if( F->scale )
	{
		if( F->storage==LIS_MATRIX_BSR && F->scale==LIS_SCALE_JACOBI )
		{
			if( F->A->matrix_type!=LIS_MATRIX_BSR )
			{
				F->err = lis_matrix_duplicate(F->A,&F->B);
				if( F->err ) return full_linear_parent(S,F->err);
				F->err = lis_matrix_set_blocksize(F->B,F->block,F->block,NULL,NULL);
				if( F->err )
				{
					lis_matrix_destroy(F->B);
					return full_linear_parent(S,F->err);
				}
				F->err = lis_matrix_set_type(F->B,F->storage);
				if( F->err )
				{
					lis_matrix_destroy(F->B);
					return full_linear_parent(S,F->err);
				}
				F->err = lis_matrix_convert(F->A,F->B);
				if( F->err )
				{
					lis_matrix_destroy(F->B);
					return full_linear_parent(S,F->err);
				}
				lis_matrix_storage_destroy(F->A);
				lis_matrix_DLU_destroy(F->A);
				lis_matrix_diag_destroy(F->A->WD);
				if( F->A->l2g_map ) lis_free( F->A->l2g_map );
				if( F->A->commtable ) lis_commtable_destroy( F->A->commtable );
				if( F->A->ranges ) lis_free( F->A->ranges );
				F->err = lis_matrix_copy_struct(F->B,F->A);
				if( F->err ) return full_linear_parent(S,F->err);
				lis_free(F->B);
			}
			F->err = lis_matrix_split(F->A);
			if( F->err ) return full_linear_parent(S,F->err);
			F->err = lis_matrix_diag_duplicate(F->A->D,&F->solver->WD);
			if( F->err ) return full_linear_parent(S,F->err);
			lis_matrix_diag_copy(F->A->D,F->solver->WD);
			lis_matrix_diag_inverse(F->solver->WD);
			lis_matrix_bscale_bsr(F->A,F->solver->WD);
			lis_vector_duplicate(F->A,&F->t);
			lis_matrix_diag_matvec(F->solver->WD,F->b,F->t);
			lis_vector_copy(F->t,F->b);
			lis_vector_destroy(F->t);
			F->t = NULL;
		}
		else
		{
			if( F->solver->d==NULL )
			{
				F->err = lis_vector_duplicate(F->A,&F->solver->d);
				if( F->err )
				{
					return full_linear_parent(S,F->err);
				}
			}
			if( F->scale==LIS_SCALE_JACOBI && F->nsolver==LIS_SOLVER_CG )
			{
				F->scale = LIS_SCALE_SYMM_DIAG;
			}
			if( !F->A->is_scaled )
			{
				F->err = lis_matrix_scale(F->A,F->b,F->solver->d,F->scale);
				if( F->err )
				{
					lis_vector_destroy(F->xx);
					F->solver->retcode = F->err;
					return full_linear_parent(S,F->err);
				}
			}
			else if( !F->b->is_scaled )
			{
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(LIS_INT full_i=0;full_i<F->n;full_i++){
					F->b->value[full_i] = F->b->value[full_i]*F->solver->d->value[full_i];
				}
			}
		}
	}

                                       
	if( F->shell_scale )
	{
		F->AA = F->shell_A;
		F->bb = F->shell_b;
	}
	else if( F->precon_type==LIS_PRECON_TYPE_IS )
	{
		if( F->nsolver < LIS_SOLVER_JACOBI || F->nsolver > LIS_SOLVER_SOR )
		{
			F->AA = F->solver->A;
			F->bb = F->solver->b;
		}
		else
		{
			F->AA = F->precon->A;
			F->bb = F->precon->Pb;
		}
	}
	else
	{
		F->AA = F->A;
		F->bb = F->b;
	}

	F->p_c_time = lis_wtime() - F->p_c_time;
	F->itime = lis_wtime();

	                    
	F->solver->A  = F->AA;
	F->solver->b  = F->bb;
	if( F->AA->matrix_type!=LIS_MATRIX_USER &&
	    F->AA->matrix_type!=LIS_MATRIX_OPERATOR )
	{
		F->err = lis_matrix_convert_self(F->solver);
		if( F->err )
		{
			lis_vector_destroy(F->xx);
			lis_solver_work_destroy(F->solver);

			if( F->shell_scale )
			{
			        F->solver->A = F->A;
			        F->solver->b = F->b;

			        lis_matrix_destroy(F->shell_A);
			        F->shell_A = NULL;

			        lis_vector_destroy(F->shell_b);
			        F->shell_b = NULL;
			}

			F->solver->retcode = F->err;
			return full_linear_parent(S,F->err);
		}
	}
	  
                                                               
                                                                   
              
   
                                                                      
                                                                     
                                                                    
                                  
    
	if( F->solver->near_nullspace_dim>0 &&
	    F->precon_type!=LIS_PRECON_TYPE_MFAMG &&
	    F->precon_type!=LIS_PRECON_TYPE_MFAMG_ILUTP &&
	    !(F->precon_type==LIS_PRECON_TYPE_HYBRID &&
	      F->solver->options[LIS_OPTIONS_PPRECON]==
	          LIS_PRECON_TYPE_MFAMG) )
	{
	        F->err =
	                lis_solver_near_nullspace_coarse_setup(
	                        F->solver,
	                        F->solver->A,
	                        F->scale);

	        if( F->err )
	        {
	                lis_vector_destroy(F->xx);
	                lis_solver_work_destroy(F->solver);

	                if( F->shell_scale )
	                {
	                        F->solver->A = F->A;
	                        F->solver->b = F->b;
	                        lis_matrix_destroy(F->shell_A);
	                        F->shell_A = NULL;
	                        lis_vector_destroy(F->shell_b);
	                        F->shell_b = NULL;
	                }

	                F->solver->retcode = F->err;
	                return full_linear_parent(S,F->err);
	        }
	}

	F->block = F->solver->A->bnr;

	if( F->A->my_rank==0 )
	{
#ifdef _LONG__DOUBLE
	  if( F->output ) printf("precision             : long double\n");
#else
 	  if( F->output ) printf("precision             : %s\n", lis_precisionname[F->precision]); 
#endif
	  if( F->output ) printf("linear solver         : %s\n", lis_solvername[F->nsolver]); 
		switch( F->precon_type )
		{
		case LIS_PRECON_TYPE_ILU:
			F->i = F->solver->options[LIS_OPTIONS_FILL];
			if( F->A->matrix_type==LIS_MATRIX_BSR || F->A->matrix_type==LIS_MATRIX_VBR )
			{
#ifdef _LONG__LONG
			  if( F->output ) sprintf(F->buf,"Block %s(%lld)",lis_preconname[F->precon_type],F->i); 
#else
			  if( F->output ) sprintf(F->buf,"Block %s(%d)",lis_preconname[F->precon_type],F->i); 
#endif
			}
			else
			{
#ifdef _LONG__LONG
			  if( F->output ) sprintf(F->buf,"%s(%lld)",lis_preconname[F->precon_type],F->i); 
#else
			  if( F->output ) sprintf(F->buf,"%s(%d)",lis_preconname[F->precon_type],F->i); 
#endif
			}
			break;		default:
		  if( F->output )
		  {
			if( F->precon_type>=LIS_PRECON_TYPE_USERDEF )
			  sprintf(F->buf,"%s",precon_register_top[F->precon_type-LIS_PRECON_TYPE_USERDEF].name);
			else
			  sprintf(F->buf,"%s",lis_preconname[F->precon_type]);
		  }
			break;
		}
		if( F->solver->options[LIS_OPTIONS_ADDS] && F->precon_type )
		{
 		  if( F->output ) printf("preconditioner        : %s + Additive Schwarz\n", F->buf);
		}
		else
		{
		  if( F->output ) printf("preconditioner        : %s\n", F->buf); 
		}
	}
	switch(F->conv_cond)
	{
	case LIS_CONV_COND_NRM2_R:
	  if( F->output ) lis_printf(F->comm,"convergence condition : ||b-Ax||_2 <= %6.1e * ||b-Ax_0||_2\n", (double)F->tol); 		  
		break;		
	case LIS_CONV_COND_NRM2_B:
		lis_vector_nrm2(F->b,&F->nrm2);
		F->nrm2 = F->nrm2*F->tol;
		if( F->output ) lis_printf(F->comm,"convergence condition : ||b-Ax||_2 <= %6.1e*||b||_2 = %6.1e\n", (double)F->tol,(double)F->nrm2);
		break;
	case LIS_CONV_COND_NRM1_B:
		lis_vector_nrm1(F->b,&F->nrm2);
		F->nrm2 = F->nrm2*F->tol_w + F->tol;
		if( F->output ) lis_printf(F->comm,"convergence condition : ||b-Ax||_1 <= %6.1e*||b||_1 + %6.1e = %6.1e\n", (double)F->tol_w,(double)F->tol,(double)F->nrm2);
		break;
	}
	if( F->AA->matrix_type==LIS_MATRIX_BSR || F->AA->matrix_type==LIS_MATRIX_BSC )
	  {
	    if( F->output ) lis_printf(F->comm,"matrix storage format : %s(%D x %D)\n", lis_storagename[F->AA->matrix_type-1],F->block,F->block); 
	  }	else if( F->AA->matrix_type==LIS_MATRIX_USER )
	  {
	    if( F->output ) lis_printf(F->comm,"matrix storage format : user/shell\n");
	  }
	  else if( F->AA->matrix_type==LIS_MATRIX_OPERATOR )
	  {
	    if( F->output ) lis_printf(F->comm,"matrix storage format : operator/shell\n");
	  }
	else
	  {
	    if( F->output ) lis_printf(F->comm,"matrix storage format : %s\n", lis_storagename[F->AA->matrix_type-1]);
	  }

	                                  
	if( F->solver->work!=NULL )
	{
	    if( lis_solver_work_compatible(F->solver) )
	    {
	        F->err = lis_solver_work_reset(F->solver);
	        if( F->err )
	        {
	            lis_vector_destroy(F->xx);

	            if( F->shell_scale )
	            {
	                    F->solver->A = F->A;
	                    F->solver->b = F->b;

	                    lis_matrix_destroy(F->shell_A);
	                    F->shell_A = NULL;

	                    lis_vector_destroy(F->shell_b);
	                    F->shell_b = NULL;
	            }
	            F->solver->retcode = F->err;
	            return full_linear_parent(S,F->err);
	        }
	    }
	    else
	    {
	        lis_solver_work_destroy(F->solver);
	    }
	}

	if( F->solver->work==NULL )
	{
	    F->err = lis_solver_malloc_work[F->nsolver](F->solver);
	    if( F->err )
	    {
	        lis_vector_destroy(F->xx);

	        if( F->shell_scale )
	        {
	                F->solver->A = F->A;
	                F->solver->b = F->b;

	                lis_matrix_destroy(F->shell_A);
	                F->shell_A = NULL;

	                lis_vector_destroy(F->shell_b);
	                F->shell_b = NULL;
	        }
	        F->solver->retcode = F->err;
	        return full_linear_parent(S,F->err);
	    }
	}
	if( F->nsolver==LIS_SOLVER_BICG && F->is_use_at )
	{
	  if( F->output ) lis_printf(F->comm,"Use Ah\n"); 
	  lis_matrix_duplicate(F->AA,&F->Ah);
	  lis_matrix_set_type(F->Ah,LIS_USE_AT_TYPE[F->AA->matrix_type]);
	  lis_matrix_convert(F->AA,F->Ah);
	  F->solver->Ah = F->Ah;
	}

	                                                         
	F->rhistory = (LIS_REAL *)lis_malloc((F->maxiter+2)*sizeof(LIS_REAL),"lis_solve::rhistory");
	if( F->rhistory==NULL )
	{
		LIS_SETERR_MEM((F->maxiter+2)*sizeof(LIS_REAL));
		lis_solver_work_destroy(F->solver);
		lis_vector_destroy(F->xx);

		if( F->shell_scale )
		{
		        F->solver->A = F->A;
		        F->solver->b = F->b;

		        lis_matrix_destroy(F->shell_A);
		        F->shell_A = NULL;

		        lis_vector_destroy(F->shell_b);
		        F->shell_b = NULL;
		}
		F->solver->retcode = LIS_ERR_OUT_OF_MEMORY;
		return full_linear_parent(S,LIS_ERR_OUT_OF_MEMORY);
	}
	F->rhistory[0] = 1.0;

	F->solver->x        = F->xx;
	F->solver->xx       = F->x;
	F->solver->precon   = F->precon;
	F->solver->rhistory = F->rhistory;

	                                                           
                                                                      
                                                     
                                 
   
	
	if (!F->solver->setup)
	  {
	                                                                    
	    F->solver->retcode = LIS_SUCCESS;
	    F->pc=1; return full_native_begin_linear(S,F->solver,F->nsolver);
resume_1: F->pc=0; F->err=S->completion;
	                                                                   
                                                                         
	    if( F->err==LIS_SUCCESS && F->solver->retcode!=LIS_SUCCESS )
	      F->err = F->solver->retcode;
	    F->solver->retcode = F->err;
	  }

	if( F->scale==LIS_SCALE_SYMM_DIAG && F->precon_type!=LIS_PRECON_TYPE_IS)
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(LIS_INT full_i=0;full_i<F->n;full_i++){
			F->x->value[full_i] = F->xx->value[full_i]*F->solver->d->value[full_i];
		}
	}
	else
	{
			lis_vector_copy(F->xx,F->x);
	}
	F->itime = lis_wtime() - F->itime - F->solver->ptime;
	F->p_i_time = F->solver->ptime;
	F->solver->ptime = F->p_c_time + F->p_i_time;
	F->solver->p_c_time = F->p_c_time;
	F->solver->p_i_time = F->p_i_time;
	F->solver->time  = F->solver->ptime + F->itime;
	F->solver->itime = F->itime;
	                                                                 
                                                                 
                                                                          
	F->t=NULL; F->xx->precision=LIS_PRECISION_DEFAULT;
if(S->forced){F->err=S->forced_code;F->solver->retcode=F->err;F->solver->resid=NAN;goto residual_done;}
	if( F->err!=LIS_SUCCESS &&
	    ((F->err!=LIS_MAXITER && F->err!=LIS_BREAKDOWN) ||
	     F->A->matrix_type==LIS_MATRIX_USER ||
	     F->A->matrix_type==LIS_MATRIX_OPERATOR ||
	     F->precon_type>=LIS_PRECON_TYPE_USERDEF ||
	     F->precon_type==LIS_PRECON_TYPE_HYBRID) )
	    goto residual_done;
	F->residual_err = lis_vector_duplicate(F->A,&F->t);
	if( F->residual_err==LIS_SUCCESS )
	{
		F->pc=2; return full_native_cert(S,F->A,F->shell_scale?F->x:F->xx,F->t);
resume_2: F->pc=0; F->residual_err=S->completion;
	}
	if( F->residual_err==LIS_SUCCESS )
		F->residual_err = lis_vector_xpay(F->b,-1.0,F->t);
	if( F->residual_err==LIS_SUCCESS && !F->shell_scale &&
	    F->scale==LIS_SCALE_SYMM_DIAG &&
	    F->precon_type!=LIS_PRECON_TYPE_IS )
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(LIS_INT full_i=0;full_i<F->n;full_i++){
			F->t->value[full_i] = F->t->value[full_i]/F->solver->d->value[full_i];
		}
	}

      
                                                             
                                                         
                                              
       
    if(
        F->err==LIS_SUCCESS
        &&
        F->residual_err==LIS_SUCCESS
        &&
        F->solver->nullspace_dim>0)
    {
        F->err =
            lis_solver_project_nullspace(
                F->solver,
                F->x);

        if(F->err)
        {
            F->solver->retcode = F->err;
        }
    }

	                                                                         
	if( F->residual_err==LIS_SUCCESS )
	{
		F->residual_err = lis_solver_get_residual[F->conv_cond](F->t,F->solver,&F->nrm2);
	}
	if( F->residual_err!=LIS_SUCCESS )
	{
		if( F->err==LIS_SUCCESS ) F->err=F->residual_err;
		F->solver->retcode = F->err;
		F->solver->resid = NAN;
		if(F->nsolver==LIS_SOLVER_GCRODR) lis_gcrodr_clear_recycle(F->solver);
		if(F->nsolver==LIS_SOLVER_FGCRODR) lis_fgcrodr_clear_recycle(F->solver);
	}

	if( !F->solver->setup && F->residual_err==LIS_SUCCESS )
	{
		F->solver->resid = F->nrm2;

		                                                              
                                                                  
		if( F->err==LIS_SUCCESS &&
		    (F->nrm2!=F->nrm2 || lis_solver_real_abs(F->nrm2)>LIS_SCALAR_MAX || F->nrm2>F->solver->tol) )
		{
			F->err = LIS_BREAKDOWN;
			F->solver->retcode = LIS_BREAKDOWN;
		}
	}

residual_done:
	                                                                 
                                                    
	if( F->err!=LIS_SUCCESS && F->err!=LIS_MAXITER && F->err!=LIS_BREAKDOWN )
	{
		if(F->nsolver==LIS_SOLVER_GCRODR) lis_gcrodr_clear_recycle(F->solver);
		if(F->nsolver==LIS_SOLVER_FGCRODR) lis_fgcrodr_clear_recycle(F->solver);
	}

	if( F->err )
	  {
	    if( F->output )
	    {
	                                                                   
                                                                     
                                               
	      const char *status_name =
	          (F->err>=LIS_SUCCESS &&
	           F->err<(LIS_INT)(sizeof(lis_returncode)/sizeof(lis_returncode[0])))
	          ? lis_returncode[F->err] : "LIS_CALLBACK_ERROR";
	      lis_printf(F->comm,"linear solver status  : %s(code=%D)\n\n",
	                 status_name,F->err);
	    }
	  }
	else
	  {
	    if( F->output ) lis_printf(F->comm,"linear solver status  : normal end\n\n"); 
	  }

	if( F->precision==LIS_PRECISION_DOUBLE )
	{
		F->solver->iter2 = F->solver->iter;
	}
	else if( F->precision==LIS_PRECISION_QUAD )
	{
		F->solver->iter2 = 0;
	}


	lis_vector_destroy(F->t);

	if( F->shell_scale )
	{
	        F->solver->A = F->A;
	        F->solver->b = F->b;

	        lis_matrix_destroy(F->shell_A);
	        F->shell_A = NULL;

	        lis_vector_destroy(F->shell_b);
	        F->shell_b = NULL;
	}

	                            
	lis_vector_destroy(F->xx);

	LIS_DEBUG_FUNC_OUT;
	return full_linear_parent(S,LIS_SUCCESS);

}
