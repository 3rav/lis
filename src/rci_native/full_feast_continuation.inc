/* Generated continuation; native LIS source copyright and license apply.
 * Sources: src/esolver/lis_esolver.c, src/esolver/lis_esolver_block.c, src/esolver/lis_esolver_feast.c, src/esolver/lis_krylov.c, src/solver/lis_solver.c, src/solver/lis_solver_bicgstab.c, src/solver/lis_solver_cg.c, src/solver/lis_solver_fgcrodr.c, src/solver/lis_solver_gcrodr.c, src/solver/lis_solver_gmres.c
 */
/* Generated from pinned original LIS functions, qualification only. */
struct full_frame_eigen { int pc;
 LIS_MATRIX A;
 LIS_MATRIX B;
 LIS_VECTOR x;
 LIS_SCALAR * evalue0;
 LIS_ESOLVER esolver;
 LIS_Comm comm;
 LIS_INT nesolver;
 LIS_INT niesolver;
 LIS_INT nigesolver;
 LIS_INT emaxiter;
 LIS_SCALAR * evalue;
 LIS_VECTOR * evector;
 LIS_REAL * resid;
 LIS_REAL * rhistory;
 LIS_INT * iter;
 LIS_INT * iter2;
 LIS_INT err;
 LIS_INT output;
 LIS_INT i;
 LIS_INT ss;
 LIS_INT mode;
 double time;
 LIS_SCALAR oshift;
 LIS_SCALAR ishift;
 LIS_INT estorage;
 LIS_INT eblock;
 LIS_MATRIX A0;
 LIS_INT eprecision;
 LIS_INT rval;
 LIS_VECTOR xx;
};
struct full_frames_feast {
 struct full_frame_eigen eigen;
};
#undef __FUNC__
#define __FUNC__ "lis_gesolve_once"
static int full_step_eigen(lis_rci_native_full *S){
struct full_frames_feast *C=S->frame; struct full_frame_eigen *F=&C->eigen;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return full_native_finish(S,LIS_BREAKDOWN);}

   
 
F->evalue = NULL;
F->evector = NULL;
F->resid = NULL;
F->rhistory = NULL;
F->iter = NULL;
F->iter2 = NULL;







F->A0 = NULL;


F->xx = NULL;

	LIS_DEBUG_FUNC_IN;

	F->comm = LIS_COMM_WORLD;

	                           
	F->err = lis_matrix_check(F->A,LIS_MATRIX_CHECK_ALL);
	if( F->B!=NULL ) F->err = lis_matrix_check(F->B,LIS_MATRIX_CHECK_ALL);	

	if( F->err ) return full_native_finish(S,F->err);
	if( F->x==NULL )
	{
		LIS_SETERR(LIS_ERR_ILL_ARG,"vector x is undefined\n");
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}
	if( F->A->n!=F->x->n )
	{
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}
	if( F->B!=NULL && F->B->n!=F->x->n )
	{
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}
	if( F->A->gn<=0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Size n(=%D) of matrix A is less than 0\n",F->A->gn);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}
	if( F->B!=NULL && F->B->gn<=0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Size n(=%D) of matrix B is less than 0\n",F->B->gn);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}

	F->nesolver = F->esolver->options[LIS_EOPTIONS_ESOLVER];
	F->niesolver = F->esolver->options[LIS_EOPTIONS_INNER_ESOLVER];
	F->nigesolver = F->esolver->options[LIS_EOPTIONS_INNER_GENERALIZED_ESOLVER];	
	F->ss = F->esolver->options[LIS_EOPTIONS_SUBSPACE];
	F->mode = F->esolver->options[LIS_EOPTIONS_MODE];
	F->emaxiter = F->esolver->options[LIS_EOPTIONS_MAXITER];
	F->oshift = F->esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
	F->ishift = F->esolver->ishift;
	F->output = F->esolver->options[LIS_EOPTIONS_OUTPUT];
	F->estorage = F->esolver->options[LIS_EOPTIONS_STORAGE];
	F->eblock = F->esolver->options[LIS_EOPTIONS_STORAGE_BLOCK];
	F->eprecision = F->esolver->options[LIS_EOPTIONS_PRECISION];
	F->esolver->eprecision = F->eprecision;
	F->rval = F->esolver->options[LIS_EOPTIONS_RVAL];

	if( F->esolver->transform < LIS_EIGEN_TRANSFORM_AUTO ||
	    F->esolver->transform > LIS_EIGEN_TRANSFORM_BUCKLING )
	{
		LIS_SETERR1(
			LIS_ERR_ILL_ARG,
			"Unknown spectral transformation %D\n",
			F->esolver->transform);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}

	if( F->esolver->transform != LIS_EIGEN_TRANSFORM_AUTO )
	{
		if( F->B==NULL )
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"spectral transformations require "
" a generalized eigenproblem\n");
			return full_native_finish(S,LIS_ERR_ILL_ARG);
		}

		  
                                                     
                                                 
                                        
     
		if( F->nesolver != LIS_ESOLVER_GII &&
		    F->nesolver != LIS_ESOLVER_GLI &&
            F->nesolver != LIS_ESOLVER_KRYLOV )
		{
			LIS_SETERR(
				LIS_ERR_NOT_IMPLEMENTED,
				"spectral transformation is not implemented "
"for this eigensolver\n");
			return full_native_finish(S,LIS_ERR_NOT_IMPLEMENTED);
		}


		if( F->esolver->transform ==
		        LIS_EIGEN_TRANSFORM_BUCKLING &&
		    F->oshift==0.0 && F->ishift==0.0 )
		{
			LIS_SETERR(
				LIS_ERR_ILL_ARG,
				"buckling transformation requires "
" a nonzero shift\n");
			return full_native_finish(S,LIS_ERR_ILL_ARG);
		}

	}

	  
                                                                 
                                                                        
    
	if( F->A->matrix_type==LIS_MATRIX_USER && F->estorage>0 )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "storage conversion is unavailable for LIS_MATRIX_USER\n");
		return full_native_finish(S,LIS_ERR_NOT_IMPLEMENTED);
	}
	if( (F->A->matrix_type==LIS_MATRIX_USER ||
	     (F->B!=NULL && F->B->matrix_type==LIS_MATRIX_USER)) &&
	    (F->oshift!=0.0 || F->ishift!=0.0) &&
	    !(F->B!=NULL &&
	  (F->nesolver==LIS_ESOLVER_GPI ||
	   F->nesolver==LIS_ESOLVER_GII ||
	   F->nesolver==LIS_ESOLVER_GCR ||
	   (F->nesolver==LIS_ESOLVER_GLI &&
	    F->esolver->transform!=LIS_EIGEN_TRANSFORM_AUTO) ||
           F->nesolver==LIS_ESOLVER_KRYLOV)) )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "matrix shifts are unavailable for LIS_MATRIX_USER\n");
		return full_native_finish(S,LIS_ERR_NOT_IMPLEMENTED);
	}

	if( F->nesolver < 1 || F->nesolver > LIS_ESOLVER_LEN )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %D (Set between 1 to %D)\n",F->nesolver, LIS_ESOLVER_LEN);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}
	  
                                                       
   
                                                            
                                                             
                                                        
                                                    
    
	if( F->esolver->block_size>1 &&
	    F->nesolver!=LIS_ESOLVER_LI &&
	    F->nesolver!=LIS_ESOLVER_GLI )
	{
		LIS_SETERR(
			LIS_ERR_NOT_IMPLEMENTED,
			"block eigensolver mode is currently available "
			"only for standard or generalized Lanczos\n");
		return full_native_finish(S,LIS_ERR_NOT_IMPLEMENTED);
	}

	if( F->esolver->block_size>1 &&
	    F->esolver->block_size>F->ss )
	{
		LIS_SETERR2(
			LIS_ERR_ILL_ARG,
			"Parameter -eblock is %D "
			"(Set less than or equal to subspace size %D)\n",
			F->esolver->block_size,
			F->ss);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}

	                                                                       
                                                 
	if( F->A->matrix_type==LIS_MATRIX_USER && F->nesolver==LIS_ESOLVER_RQI )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "Rayleigh quotient iteration requires matrix shifts for LIS_MATRIX_USER\n");
		return full_native_finish(S,LIS_ERR_NOT_IMPLEMENTED);
	}

	if( F->B!=NULL &&
	    (F->A->matrix_type==LIS_MATRIX_USER || F->B->matrix_type==LIS_MATRIX_USER) &&
	    F->nesolver==LIS_ESOLVER_GRQI )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "generalized Rayleigh quotient iteration requires matrix shifts for LIS_MATRIX_USER\n");
		return full_native_finish(S,LIS_ERR_NOT_IMPLEMENTED);
	}

	  
                                      
  
                                                                                                              
                         
  

                                                                                 
  
                                                                                                                             
                         
  

                                                                                               
  
                                                                                                                         
                         
  

                                                                                               
  
                                                                                                                         
                         
  
   

	if (( F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GSI ) && F->ss > F->A->gn )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_SUBSPACE is %D (Set less than or equal to matrix size %D for Subspace)\n", F->ss, F->A->gn);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}

	if (( F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_AI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GLI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GAI) && F->ss > F->A->gn )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_SUBSPACE is %D (Set less than or equal to matrix size %D for Lanczos and Arnoldi)\n", F->ss, F->A->gn);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}

	if (( F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GSI ) && F->mode >= F->ss )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_MODE is %D (Set less than subspace size %D for Subspace)\n", F->mode, F->ss);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}

	if (( F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_AI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GLI || F->esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_GAI) && F->mode >= F->ss )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_MODE is %D (Set less than subspace size %D for Lanczos and Arnoldi)\n", F->mode, F->ss);
		return full_native_finish(S,LIS_ERR_ILL_ARG);
	}


	F->err = lis_esolver_check_params[F->nesolver](F->esolver);
	if( F->err )
	{
		F->esolver->retcode = F->err;
		return full_native_finish(S,F->err);
	}
	                         

	                                                                    
                                                               
	F->evalue = (LIS_SCALAR *)lis_malloc((F->ss+2)*sizeof(LIS_SCALAR),"lis_esolve::evalue");
	if( F->evalue==NULL ) goto result_oom;
	for(F->i=0;F->i<F->ss+2;F->i++) F->evalue[F->i] = (LIS_SCALAR)NAN;
	F->evalue[0] = 1.0;
	F->evalue[F->ss-1] = 1.0;

	F->resid = (LIS_REAL *)lis_malloc((F->ss+2)*sizeof(LIS_REAL),"lis_esolve::resid");
	if( F->resid==NULL ) goto result_oom;
	F->iter = (LIS_INT *)lis_malloc((F->ss+2)*sizeof(LIS_INT),"lis_esolve::iter");
	if( F->iter==NULL ) goto result_oom;
	F->iter2 = (LIS_INT *)lis_malloc((F->ss+2)*sizeof(LIS_INT),"lis_esolve::iter2");
	if( F->iter2==NULL ) goto result_oom;
	for(F->i=0;F->i<F->ss+2;F->i++)
	{
		F->resid[F->i] = NAN;
		F->iter[F->i] = 0;
		F->iter2[F->i] = 0;
	}

	                           
		F->err = lis_vector_duplicate(F->A,&F->xx);
	if( F->err ) goto result_failure;
	if( F->esolver->options[LIS_EOPTIONS_INITGUESS_ONES] )
	{
	  if( F->output ) lis_printf(F->comm,"initial vector x      : all components set to 1\n");
			lis_vector_set_all(1.0,F->xx);
	}
	else
	{
	  if( F->output ) lis_printf(F->comm,"initial vector x      : user defined\n"); 
			lis_vector_copy(F->x,F->xx);
	}

	                                                      
	F->evector = (LIS_VECTOR *)lis_malloc((F->ss+2)*sizeof(LIS_VECTOR),"lis_esolve::evector");
	if( F->evector==NULL ) goto result_oom;
	memset(F->evector,0,(F->ss+2)*sizeof(LIS_VECTOR));
	F->rhistory = (LIS_REAL *)lis_malloc((F->emaxiter+2)*sizeof(LIS_REAL),"lis_esolve::rhistory");
	if( F->rhistory==NULL ) goto result_oom;

	                    
	if( F->estorage>0 && F->A->matrix_type!=F->estorage )
	{
		F->err = lis_matrix_duplicate(F->A,&F->A0);
		if( F->err ) goto result_failure;
		lis_matrix_set_blocksize(F->A0,F->eblock,F->eblock,NULL,NULL);
		lis_matrix_set_type(F->A0,F->estorage);
		F->err = lis_matrix_convert(F->A,F->A0);
		if( F->err ) goto result_failure;
		lis_matrix_storage_destroy(F->A);
		lis_matrix_DLU_destroy(F->A);
		lis_matrix_diag_destroy(F->A->WD);
		if( F->A->l2g_map ) lis_free( F->A->l2g_map );
		if( F->A->commtable ) lis_commtable_destroy( F->A->commtable );
		if( F->A->ranges ) lis_free( F->A->ranges );
		F->err = lis_matrix_copy_struct(F->A0,F->A);
		if( F->err ) goto result_failure;
		lis_free(F->A0);
		F->A0 = NULL;
	}

	                                                                       
	if( F->esolver->evalue ) lis_free(F->esolver->evalue);
	if( F->esolver->resid ) lis_free(F->esolver->resid);
	if( F->esolver->iter ) lis_free(F->esolver->iter);
	if( F->esolver->iter2 ) lis_free(F->esolver->iter2);
	if( F->esolver->rhistory ) lis_free(F->esolver->rhistory);
	if( F->esolver->evector )
	{
	          
                                                              
                                                 
            
	        for(F->i=0;F->i<F->esolver->evector_owned_count;F->i++)
	        {
	                if( F->esolver->evector[F->i] )
	                {
	                        lis_vector_destroy(
	                                F->esolver->evector[F->i]);
	                }
	        }

	        lis_free(F->esolver->evector);
	        F->esolver->evector = NULL;
	        F->esolver->evector_owned_count = 0;
	}
	F->esolver->A        = F->A;
	F->esolver->B        = F->B;
	F->esolver->evalue   = F->evalue;
	F->esolver->x        = F->x;
	F->esolver->evector  = F->evector;

	  
                                                           
                                                   
    
	if( F->nesolver==LIS_ESOLVER_LI ||
	    F->nesolver==LIS_ESOLVER_AI ||
	    F->nesolver==LIS_ESOLVER_SI ||
	    F->nesolver==LIS_ESOLVER_GLI ||
	    F->nesolver==LIS_ESOLVER_GSI ||
	    F->nesolver==LIS_ESOLVER_GAI ||
	    F->nesolver==LIS_ESOLVER_FEAST ||
            F->nesolver==LIS_ESOLVER_KRYLOV )
	{
	        F->esolver->evector_owned_count = F->ss+2;
	}
	else
	{
	        F->esolver->evector_owned_count = 0;
	}
	F->rhistory[0]       = 1.0;
	F->esolver->rhistory = F->rhistory;
	F->esolver->resid    = F->resid;
	F->esolver->iter     = F->iter;
	F->esolver->iter2    = F->iter2;

#ifdef _LONG__DOUBLE
	if ( F->output ) lis_printf(F->comm,"precision             : long double\n");
#else
	if ( F->output ) lis_printf(F->comm,"precision             : %s\n", lis_eprecisionname[F->eprecision]);
#endif
	if ( F->output ) lis_printf(F->comm,"eigensolver           : %s\n", lis_esolvername[F->nesolver]);

	if ( F->output ) lis_printf(F->comm,"convergence condition : ||lx-(B^-1)Ax||_2 <= %6.1e * ||lx||_2\n", (double)F->esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]); 

	if( F->A->matrix_type==LIS_MATRIX_BSR || F->A->matrix_type==LIS_MATRIX_BSC )
	  {
	    if ( F->output ) lis_printf(F->comm,"matrix storage format : %s(%D x %D)\n", lis_estoragename[F->A->matrix_type-1],F->eblock,F->eblock); 
	  }
	else
	  {
	    if( F->A->matrix_type==LIS_MATRIX_USER )
	      {
	        if ( F->output ) lis_printf(F->comm,"matrix storage format : user/shell\n");
	      }
	    else
	      {
	        if ( F->output ) lis_printf(F->comm,"matrix storage format : %s\n", lis_estoragename[F->A->matrix_type-1]);
	      }
	  }
	
	F->time = lis_wtime();

	F->esolver->ptime = 0;
	F->esolver->itime = 0;
	F->esolver->p_c_time = 0;
	F->esolver->p_i_time = 0;

	  
                                                                          
                                                                           
                                                                    
           
        if(F->esolver->work!=NULL)
        {
                F->err=lis_esolver_work_destroy(F->esolver);
                if(F->err)
                {
                        F->esolver->x=NULL;
                        F->esolver->xx=NULL;
                        lis_vector_destroy(F->xx);
                        F->esolver->retcode=F->err;
                        return full_native_finish(S,F->err);
                }
        }

	                        
	F->err = lis_esolver_malloc_work[F->nesolver](F->esolver);
	if( F->err )
	{
	  F->esolver->x        = NULL;
	  F->esolver->xx       = NULL;
	  lis_vector_destroy(F->xx);
	  F->esolver->retcode = F->err;
	  return full_native_finish(S,F->err);
	}

	F->esolver->x        = F->xx;
	F->esolver->xx       = F->x;

	                     
		F->pc=1; return full_native_begin_feast(S,F->esolver);
resume_1: F->pc=0; F->err=S->completion;
	if(lis_callback_scope_status()) F->err=lis_callback_scope_status();
	F->esolver->retcode = F->err;

	                                                                  
                                                                           
	if( lis_callback_scope_status() ||
	    (F->err!=LIS_SUCCESS && F->err!=LIS_MAXITER && F->err!=LIS_BREAKDOWN) )
	{
		for(F->i=0;F->i<F->ss+2;F->i++)
		{
			F->esolver->evalue[F->i] = (LIS_SCALAR)NAN;
			F->esolver->resid[F->i] = NAN;
			if(F->i!=0) F->esolver->iter[F->i] = 0;
			F->esolver->iter2[F->i] = 0;
		}
	}

	*F->evalue0 = F->esolver->evalue[0];
	lis_vector_copy(F->esolver->x, F->x);

	F->esolver->time = lis_wtime() - F->time; 

	if( F->err )
	  {
	    if( F->output )
	    {
	                                                                               
	      const char *name = (F->err>=LIS_SUCCESS &&
	          F->err<(LIS_INT)(sizeof(lis_ereturncode)/sizeof(lis_ereturncode[0])))
	          ? lis_ereturncode[F->err] : "LIS_CALLBACK_ERROR";
	      lis_printf(F->comm,"eigensolver status    : %s(code=%D)\n\n",name,F->err);
	    }
	  }
	else
	  {
	    if ( F->output ) lis_printf(F->comm,"eigensolver status    : normal end\n\n");
	  }

	if( F->eprecision==LIS_PRECISION_DOUBLE )
	  {
	    F->esolver->iter2[F->mode] = F->esolver->iter[F->mode];
	  }
	else if( F->eprecision==LIS_PRECISION_QUAD )
	  {
	    F->esolver->iter2[F->mode] = 0;
	  }

	                                            
	F->esolver->x        = NULL;
	F->esolver->xx       = NULL;
	lis_vector_destroy(F->xx);

	LIS_DEBUG_FUNC_OUT;
	return full_native_finish(S,LIS_SUCCESS);

result_oom:
	LIS_SETERR(LIS_OUT_OF_MEMORY,"eigensolver result allocation failed\n");
	F->err = LIS_OUT_OF_MEMORY;
result_failure:
	if(F->A0) lis_matrix_destroy(F->A0);
	if(F->xx) lis_vector_destroy(F->xx);
	if(F->evalue) lis_free(F->evalue);
	if(F->resid) lis_free(F->resid);
	if(F->iter) lis_free(F->iter);
	if(F->iter2) lis_free(F->iter2);
	if(F->evector) lis_free(F->evector);
	if(F->rhistory) lis_free(F->rhistory);
	F->esolver->retcode = F->err;
	if(F->evalue0) *F->evalue0 = (LIS_SCALAR)NAN;
	LIS_DEBUG_FUNC_OUT;
	return full_native_finish(S,F->err);

}
