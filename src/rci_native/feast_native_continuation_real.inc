/* Generated continuation; native LIS source copyright and license apply.
 * Sources: src/esolver/lis_esolver.c, src/esolver/lis_esolver_block.c, src/esolver/lis_esolver_feast.c, src/esolver/lis_krylov.c, src/solver/lis_solver.c, src/solver/lis_solver_bicgstab.c, src/solver/lis_solver_cg.c, src/solver/lis_solver_fgcrodr.c, src/solver/lis_solver_gcrodr.c, src/solver/lis_solver_gmres.c
 */
/* Generated, qualification only. Native source copyrights retained. */
extern LIS_RCI_PRIVATE LIS_INT lis_eigen_block_rank_transform(
    LIS_INT n,
    LIS_REAL rank_tol,
    LIS_SCALAR *gram,
    LIS_SCALAR *eval,
    LIS_SCALAR *evec,
    LIS_SCALAR *transform,
    LIS_INT *active_rank);
extern LIS_RCI_PRIVATE LIS_INT lis_eigen_block_zero_outputs(
    LIS_INT capacity,
    LIS_INT start,
    LIS_VECTOR source,
    LIS_VECTOR Y[]);
struct nf_projector { int pc,parent;
 LIS_MATRIX A;
 LIS_MATRIX B;
 LIS_REAL emin;
 LIS_REAL emax;
 LIS_REAL ratio;
 LIS_INT npoints;
 LIS_VECTOR x_real;
 LIS_VECTOR x_imag;
 LIS_VECTOR q_real;
 LIS_VECTOR q_imag;
 LIS_SOLVER solver;
 LIS_REAL z_real[1024];
 LIS_REAL z_imag[1024];
 LIS_REAL w_real[1024];
 LIS_REAL w_imag[1024];
 LIS_VECTOR rhs_real;
 LIS_VECTOR rhs_imag;
 LIS_VECTOR y_real;
 LIS_VECTOR y_imag;
 LIS_INT i;
 LIS_INT err;
};
struct nf_bgram { int pc,parent;
 LIS_MATRIX B;
 LIS_INT nx;
 LIS_VECTOR * X;
 LIS_INT ny;
 LIS_VECTOR * Y;
 LIS_VECTOR * BY;
 LIS_SCALAR * gram;
 LIS_INT i;
 LIS_INT j;
 LIS_INT err;
};
struct nf_brank { int pc,parent;
 LIS_MATRIX B;
 LIS_INT capacity;
 LIS_INT nominal_rank;
 LIS_VECTOR * W;
 LIS_VECTOR * BW;
 LIS_VECTOR * Q;
 LIS_VECTOR * BQ;
 LIS_REAL rank_tol;
 LIS_INT * active_rank;
 LIS_SCALAR * gram;
 LIS_SCALAR * eval;
 LIS_SCALAR * evec;
 LIS_SCALAR * transform;
 LIS_INT i;
 LIS_INT j;
 LIS_INT k;
 LIS_INT err;
 LIS_INT rank;
};
struct nf_filter { int pc,parent;
 LIS_MATRIX A;
 LIS_MATRIX B;
 LIS_REAL emin;
 LIS_REAL emax;
 LIS_REAL ratio;
 LIS_INT npoints;
 LIS_INT block_size;
 LIS_VECTOR * X;
 LIS_VECTOR * Q;
 LIS_VECTOR * BQ;
 LIS_REAL rank_tol;
 LIS_INT * active_rank;
 LIS_SOLVER solver;
 LIS_VECTOR * W;
 LIS_VECTOR * BW;
 LIS_VECTOR zero_imag;
 LIS_VECTOR q_real;
 LIS_VECTOR q_imag;
 LIS_SCALAR * gram;
 LIS_SCALAR * eval;
 LIS_SCALAR * evec;
 LIS_SCALAR * transform;
 size_t vector_bytes;
 size_t square_count;
 size_t square_bytes;
 size_t eval_bytes;
 LIS_INT i;
 LIS_INT err;
};
struct nf_rr { int pc,parent;
 LIS_MATRIX A;
 LIS_INT active_rank;
 LIS_VECTOR * Q;
 LIS_SCALAR * ritz_value;
 LIS_VECTOR * ritz_vector;
 LIS_INT * projected_niter;
 LIS_REAL * projected_error;
 LIS_VECTOR * AQ;
 LIS_SCALAR * projected;
 LIS_SCALAR * coeff;
 LIS_SCALAR value;
 size_t vector_bytes;
 size_t square_count;
 size_t square_bytes;
 LIS_INT i;
 LIS_INT j;
 LIS_INT k;
 LIS_INT err;
};
struct nf_standard { int pc,parent;
 LIS_MATRIX A;
 LIS_VECTOR x;
 LIS_SCALAR lambda;
 LIS_VECTOR Ax;
 LIS_VECTOR r;
 LIS_REAL * resid;
 LIS_INT err;
 LIS_REAL anrm;
 LIS_REAL xnrm;
 LIS_REAL rnrm;
 LIS_REAL denom;
 LIS_REAL lmag;
};
struct nf_generalized { int pc,parent;
 LIS_MATRIX A;
 LIS_MATRIX B;
 LIS_VECTOR x;
 LIS_SCALAR lambda;
 LIS_VECTOR Ax;
 LIS_VECTOR Bx;
 LIS_VECTOR r;
 LIS_REAL * resid;
 LIS_INT err;
 LIS_REAL anrm;
 LIS_REAL bnrm;
 LIS_REAL rnrm;
 LIS_REAL denom;
 LIS_REAL lmag;
};
struct nf_cert { int pc,parent;
 LIS_MATRIX A;
 LIS_MATRIX B;
 LIS_REAL emin;
 LIS_REAL emax;
 LIS_REAL tol;
 LIS_INT active_rank;
 LIS_SCALAR * ritz_value;
 LIS_VECTOR * ritz_vector;
 LIS_INT previous_count;
 LIS_REAL previous_trace;
 LIS_INT * selected_index;
 LIS_REAL * residual;
 LIS_INT * selected_count;
 LIS_REAL * trace;
 LIS_REAL * trace_error;
 LIS_REAL * max_residual;
 LIS_INT * converged;
 LIS_VECTOR Ax;
 LIS_VECTOR Bx;
 LIS_VECTOR r;
 LIS_REAL lambda_real;
 LIS_REAL lambda_imag;
 LIS_REAL imag_limit;
 LIS_REAL interval_scale;
 LIS_REAL interval_slack;
 LIS_REAL pair_residual;
 LIS_REAL local_trace;
 LIS_REAL local_trace_error;
 LIS_REAL local_max_residual;
 LIS_REAL trace_scale;
 LIS_INT count;
 LIS_INT i;
 LIS_INT j;
 LIS_INT err;
};
struct nf_feast { int pc,parent;
 LIS_ESOLVER esolver;
 LIS_MATRIX A;
 LIS_MATRIX B;
 LIS_VECTOR * X;
 LIS_VECTOR * Q;
 LIS_VECTOR * BQ;
 LIS_VECTOR * Ritz;
 LIS_SOLVER shifted_solver;
 LIS_SCALAR * ritz_value;
 LIS_INT * selected_index;
 LIS_REAL * selected_residual;
 LIS_REAL tol;
 LIS_REAL rank_tol;
 LIS_REAL previous_trace;
 LIS_REAL trace;
 LIS_REAL trace_error;
 LIS_REAL max_residual;
 LIS_REAL history_value;
 LIS_REAL projected_error;
 LIS_INT ss;
 LIS_INT emaxiter;
 LIS_INT previous_count;
 LIS_INT selected_count;
 LIS_INT active_rank;
 LIS_INT projected_niter;
 LIS_INT converged;
 LIS_INT first;
 LIS_INT iter;
 LIS_INT i;
 LIS_INT err;
 size_t scalar_bytes;
 size_t integer_bytes;
 size_t real_bytes;
};
struct lis_rci_native_feast { int active,status,strict,forced; LIS_INT completion,native; uint64_t token; LIS_ESOLVER esolver; LIS_MATRIX owner_a,owner_b; LIS_UNSIGNED_INT as, av, bs, bv; lis_rci_native_feast_request request; LIS_VECTOR one_in,one_out;
 struct nf_projector projector;
 struct nf_bgram bgram;
 struct nf_brank brank;
 struct nf_filter filter;
 struct nf_rr rr;
 struct nf_standard standard;
 struct nf_generalized generalized;
 struct nf_cert cert;
 struct nf_feast feast;
};
static int nf_return(lis_rci_native_feast *,int,LIS_INT);
static int nf_block(lis_rci_native_feast *,int,LIS_MATRIX,LIS_INT,LIS_VECTOR *,LIS_VECTOR *);
static int nf_one(lis_rci_native_feast *,int,LIS_MATRIX,LIS_VECTOR,LIS_VECTOR);
static int nf_shift(lis_rci_native_feast *,LIS_MATRIX,LIS_MATRIX,LIS_REAL,LIS_REAL,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER);
#undef __FUNC__
#define __FUNC__ "lis_efeast_projector_apply_split"
static int nf_step_projector(lis_rci_native_feast *K){
struct nf_projector *F=&K->projector;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;case 4:goto resume_4;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}















        LIS_DEBUG_FUNC_IN;

        F->rhs_real = NULL;
        F->rhs_imag = NULL;

        F->y_real = NULL;
        F->y_imag = NULL;

        if( F->A==NULL ||
            F->x_real==NULL ||
            F->x_imag==NULL ||
            F->q_real==NULL ||
            F->q_imag==NULL ||
            F->solver==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projector received an undefined argument\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( !lis_efeast_shifted_vector_compatible(
                    F->A,
                    F->x_real) ||
            !lis_efeast_shifted_vector_compatible(
                    F->A,
                    F->x_imag) ||
            !lis_efeast_shifted_vector_compatible(
                    F->A,
                    F->q_real) ||
            !lis_efeast_shifted_vector_compatible(
                    F->A,
                    F->q_imag) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projector requires compatible split vectors\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( F->B!=NULL &&
            (F->B->n!=F->A->n ||
             F->B->gn!=F->A->gn ||
             F->B->my_rank!=F->A->my_rank ||
             F->B->nprocs!=F->A->nprocs) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projector requires compatible A and B layouts\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        F->err = lis_efeast_contour_gauss(
                F->emin,
                F->emax,
                F->npoints,
                F->ratio,
                F->z_real,
                F->z_imag,
                F->w_real,
                F->w_imag);

        if( F->err )
        {
                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,F->err);
        }

        F->err = lis_vector_duplicate(
                F->A,
                &F->rhs_real);

        if( F->err ) goto cleanup;

        F->err = lis_vector_duplicate(
                F->A,
                &F->rhs_imag);

        if( F->err ) goto cleanup;

        F->err = lis_vector_duplicate(
                F->A,
                &F->y_real);

        if( F->err ) goto cleanup;

        F->err = lis_vector_duplicate(
                F->A,
                &F->y_imag);

        if( F->err ) goto cleanup;

          
                                                       
           
        if( F->B!=NULL )
        {
                F->pc=1;
return nf_one(K,NF_B,F->B,F->x_real,F->rhs_real);
resume_1: F->pc=0; F->err=K->completion;

                if( F->err ) goto cleanup;

                F->pc=2;
return nf_one(K,NF_B,F->B,F->x_imag,F->rhs_imag);
resume_2: F->pc=0; F->err=K->completion;

                if( F->err ) goto cleanup;
        }
        else
        {
                F->err = lis_vector_copy(
                        F->x_real,
                        F->rhs_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_copy(
                        F->x_imag,
                        F->rhs_imag);

                if( F->err ) goto cleanup;
        }

          
                                                                  
                                          
           
        F->err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                F->q_real);

        if( F->err ) goto cleanup;

        F->err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                F->q_imag);

        if( F->err ) goto cleanup;

        for(F->i=0;F->i<F->npoints;F->i++)
        {
                  
                              
                  
                               
                   
                F->err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        F->y_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        F->y_imag);

                if( F->err ) goto cleanup;

                F->pc=3;
return nf_shift(K,F->A,F->B,F->z_real[F->i],F->z_imag[F->i],F->rhs_real,F->rhs_imag,F->y_real,F->y_imag,F->solver);
resume_3: F->pc=0; F->err=K->completion;

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)F->w_real[F->i],
                        F->y_real,
                        F->q_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)(-F->w_imag[F->i]),
                        F->y_imag,
                        F->q_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)F->w_imag[F->i],
                        F->y_real,
                        F->q_imag);

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)F->w_real[F->i],
                        F->y_imag,
                        F->q_imag);

                if( F->err ) goto cleanup;

                  
                              
                  
                                           
                  
                                                                       
                                                                      
                                      
                   
                F->err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        F->y_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_set_all(
                        (LIS_SCALAR)0.0,
                        F->y_imag);

                if( F->err ) goto cleanup;

                F->pc=4;
return nf_shift(K,F->A,F->B,F->z_real[F->i],-F->z_imag[F->i],F->rhs_real,F->rhs_imag,F->y_real,F->y_imag,F->solver);
resume_4: F->pc=0; F->err=K->completion;

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)F->w_real[F->i],
                        F->y_real,
                        F->q_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)F->w_imag[F->i],
                        F->y_imag,
                        F->q_real);

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)(-F->w_imag[F->i]),
                        F->y_real,
                        F->q_imag);

                if( F->err ) goto cleanup;

                F->err = lis_vector_axpy(
                        (LIS_SCALAR)F->w_real[F->i],
                        F->y_imag,
                        F->q_imag);

                if( F->err ) goto cleanup;
        }

        F->err = LIS_SUCCESS;

cleanup:

        if( F->y_imag!=NULL )
        {
                lis_vector_destroy(
                        F->y_imag);
        }

        if( F->y_real!=NULL )
        {
                lis_vector_destroy(
                        F->y_real);
        }

        if( F->rhs_imag!=NULL )
        {
                lis_vector_destroy(
                        F->rhs_imag);
        }

        if( F->rhs_real!=NULL )
        {
                lis_vector_destroy(
                        F->rhs_real);
        }

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_eigen_block_bgram"
static int nf_step_bgram(lis_rci_native_feast *K){
struct nf_bgram *F=&K->bgram;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}



    LIS_DEBUG_FUNC_IN;

    if( F->nx<=0 || F->ny<=0 ||
        F->X==NULL || F->Y==NULL ||
        F->gram==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block B-Gram arguments\n");

        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
    }

    if( F->B==NULL )
    {
        return nf_return(K,F->parent,lis_eigen_block_gram(
            F->nx,
            F->X,
            F->ny,
            F->Y,
            F->gram));
    }

    if( F->BY==NULL )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "missing B-Gram work vectors\n");

        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
    }

    for(F->i=0;F->i<F->nx;F->i++)
    {
        if( F->X[F->i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null left block vector\n");

            return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }
    }

    for(F->j=0;F->j<F->ny;F->j++)
    {
        if( F->Y[F->j]==NULL || F->BY[F->j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null right or B-work block vector\n");

            return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }
    }

    F->pc=1;
return nf_block(K,NF_B,F->B,F->ny,F->Y,F->BY);
resume_1: F->pc=0; F->err=K->completion;

    if( F->err )
    {
        return nf_return(K,F->parent,F->err);
    }

    for(F->i=0;F->i<F->nx;F->i++)
    {
        for(F->j=0;F->j<F->ny;F->j++)
        {
            F->err = lis_vector_dot(
                F->X[F->i],
                F->BY[F->j],
                &F->gram[F->i*F->ny+F->j]);

            if( F->err )
            {
                return nf_return(K,F->parent,F->err);
            }
        }
    }

    LIS_DEBUG_FUNC_OUT;

    return nf_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_eigen_block_brank_deflate"
static int nf_step_brank(lis_rci_native_feast *K){
struct nf_brank *F=&K->brank;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}




    LIS_DEBUG_FUNC_IN;

    if( F->B==NULL ||
        F->capacity<=0 ||
        F->nominal_rank<=0 ||
        F->nominal_rank>F->capacity ||
        F->W==NULL ||
        F->BW==NULL ||
        F->Q==NULL ||
        F->BQ==NULL ||
        F->active_rank==NULL ||
        F->gram==NULL ||
        F->eval==NULL ||
        F->evec==NULL ||
        F->transform==NULL ||
        !isfinite((double)F->rank_tol) ||
        F->rank_tol<0.0 ||
        F->rank_tol>=1.0 )
    {
        LIS_SETERR(
            LIS_ERR_ILL_ARG,
            "invalid block B-rank-deflation arguments\n");

        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
    }

    *F->active_rank = 0;

    for(F->i=0;F->i<F->nominal_rank;F->i++)
    {
        if( F->W[F->i]==NULL ||
            F->BW[F->i]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null B-rank-deflation candidate or work vector\n");

            return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        for(F->k=0;F->k<F->nominal_rank;F->k++)
        {
            if( F->BW[F->i]==F->W[F->k] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "B-rank deflation requires separate BW vectors\n");

                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
            }
        }
    }

    for(F->j=0;F->j<F->capacity;F->j++)
    {
        if( F->Q[F->j]==NULL ||
            F->BQ[F->j]==NULL )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "null B-rank-deflation output vector\n");

            return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( F->Q[F->j]==F->BQ[F->j] )
        {
            LIS_SETERR(
                LIS_ERR_ILL_ARG,
                "B-rank deflation requires separate Q and BQ vectors\n");

            return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        for(F->i=0;F->i<F->nominal_rank;F->i++)
        {
            if( F->Q[F->j]==F->W[F->i] ||
                F->Q[F->j]==F->BW[F->i] ||
                F->BQ[F->j]==F->W[F->i] ||
                F->BQ[F->j]==F->BW[F->i] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "B-rank deflation requires "
                    "out-of-place output vectors\n");

                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
            }
        }

        for(F->k=0;F->k<F->j;F->k++)
        {
            if( F->Q[F->j]==F->Q[F->k] ||
                F->BQ[F->j]==F->BQ[F->k] ||
                F->Q[F->j]==F->BQ[F->k] ||
                F->BQ[F->j]==F->Q[F->k] )
            {
                LIS_SETERR(
                    LIS_ERR_ILL_ARG,
                    "B-rank deflation requires distinct outputs\n");

                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
            }
        }
    }

      
                                                            
       
    F->pc=1;
K->bgram.B = F->B;
K->bgram.nx = F->nominal_rank;
K->bgram.X = F->W;
K->bgram.ny = F->nominal_rank;
K->bgram.Y = F->W;
K->bgram.BY = F->BW;
K->bgram.gram = F->gram;
K->bgram.pc=0; K->bgram.parent=3;
K->active=2; return 3;
resume_1: F->pc=0; F->err=K->completion;

    if( F->err )
    {
        return nf_return(K,F->parent,F->err);
    }

    F->err = lis_eigen_block_rank_transform(
        F->nominal_rank,
        F->rank_tol,
        F->gram,
        F->eval,
        F->evec,
        F->transform,
        &F->rank);

    if( F->err )
    {
        return nf_return(K,F->parent,F->err);
    }

    if( F->rank>0 )
    {
        F->err =
            lis_eigen_block_linear_combination(
                F->nominal_rank,
                F->W,
                F->rank,
                F->Q,
                F->transform,
                F->nominal_rank);

        if( F->err )
        {
            return nf_return(K,F->parent,F->err);
        }

          
                                          
          
                                                            
           
        F->err =
            lis_eigen_block_linear_combination(
                F->nominal_rank,
                F->BW,
                F->rank,
                F->BQ,
                F->transform,
                F->nominal_rank);

        if( F->err )
        {
            return nf_return(K,F->parent,F->err);
        }
    }

    F->err = lis_eigen_block_zero_outputs(
        F->capacity,
        F->rank,
        F->W[0],
        F->Q);

    if( F->err )
    {
        return nf_return(K,F->parent,F->err);
    }

    F->err = lis_eigen_block_zero_outputs(
        F->capacity,
        F->rank,
        F->BW[0],
        F->BQ);

    if( F->err )
    {
        return nf_return(K,F->parent,F->err);
    }

    *F->active_rank = F->rank;

    LIS_DEBUG_FUNC_OUT;

    return nf_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_efeast_filtered_block_basis"
static int nf_step_filter(lis_rci_native_feast *K){
struct nf_filter *F=&K->filter;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}





















        LIS_DEBUG_FUNC_IN;

        F->W = NULL;
        F->BW = NULL;

        F->zero_imag = NULL;
        F->q_real = NULL;
        F->q_imag = NULL;

        F->gram = NULL;
        F->eval = NULL;
        F->evec = NULL;
        F->transform = NULL;

        if( F->active_rank!=NULL )
        {
                *F->active_rank = 0;
        }

        if( F->A==NULL ||
            F->block_size<=0 ||
            F->X==NULL ||
            F->Q==NULL ||
            F->active_rank==NULL ||
            F->solver==NULL ||
            !isfinite((double)F->rank_tol) ||
            F->rank_tol<0.0 ||
            F->rank_tol>=1.0 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST filtered block received invalid arguments\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( F->B!=NULL && F->BQ==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Generalized FEAST filtered block requires BQ outputs\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( F->B==NULL && F->BQ!=NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Standard FEAST filtered block requires BQ to be NULL\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( F->B!=NULL &&
            (F->B->n!=F->A->n ||
             F->B->gn!=F->A->gn ||
             F->B->my_rank!=F->A->my_rank ||
             F->B->nprocs!=F->A->nprocs) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST filtered block requires compatible A and B layouts\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        for(F->i=0;F->i<F->block_size;F->i++)
        {
                if( F->X[F->i]==NULL ||
                    F->Q[F->i]==NULL ||
                    !lis_efeast_shifted_vector_compatible(
                            F->A,
                            F->X[F->i]) ||
                    !lis_efeast_shifted_vector_compatible(
                            F->A,
                            F->Q[F->i]) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST filtered block contains an incompatible vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                }

                if( F->B!=NULL &&
                    (F->BQ[F->i]==NULL ||
                     !lis_efeast_shifted_vector_compatible(
                            F->A,
                            F->BQ[F->i])) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST filtered block contains an incompatible BQ vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                }
        }

        F->vector_bytes =
                (size_t)F->block_size
                *
                sizeof(LIS_VECTOR);

        F->square_count =
                (size_t)F->block_size
                *
                (size_t)F->block_size;

        F->square_bytes =
                F->square_count
                *
                sizeof(LIS_SCALAR);

        F->eval_bytes =
                (size_t)F->block_size
                *
                sizeof(LIS_SCALAR);

        F->W = (LIS_VECTOR *)lis_malloc(
                F->vector_bytes,
                "lis_efeast_filtered_block_basis::W");

        if( F->W==NULL )
        {
                LIS_SETERR_MEM(
                        F->vector_bytes);

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_OUT_OF_MEMORY);
        }

        memset(
                F->W,
                0,
                F->vector_bytes);

        if( F->B!=NULL )
        {
                F->BW = (LIS_VECTOR *)lis_malloc(
                        F->vector_bytes,
                        "lis_efeast_filtered_block_basis::BW");

                if( F->BW==NULL )
                {
                        lis_free(
                                F->W);

                        LIS_SETERR_MEM(
                                F->vector_bytes);

                        LIS_DEBUG_FUNC_OUT;
                        return nf_return(K,F->parent,LIS_OUT_OF_MEMORY);
                }

                memset(
                        F->BW,
                        0,
                        F->vector_bytes);
        }

        F->gram = (LIS_SCALAR *)lis_malloc(
                F->square_bytes,
                "lis_efeast_filtered_block_basis::gram");

        if( F->gram==NULL )
        {
                F->err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        F->square_bytes);
                goto cleanup;
        }

        F->eval = (LIS_SCALAR *)lis_malloc(
                F->eval_bytes,
                "lis_efeast_filtered_block_basis::eval");

        if( F->eval==NULL )
        {
                F->err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        F->eval_bytes);
                goto cleanup;
        }

        F->evec = (LIS_SCALAR *)lis_malloc(
                F->square_bytes,
                "lis_efeast_filtered_block_basis::evec");

        if( F->evec==NULL )
        {
                F->err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        F->square_bytes);
                goto cleanup;
        }

        F->transform = (LIS_SCALAR *)lis_malloc(
                F->square_bytes,
                "lis_efeast_filtered_block_basis::transform");

        if( F->transform==NULL )
        {
                F->err = LIS_OUT_OF_MEMORY;
                LIS_SETERR_MEM(
                        F->square_bytes);
                goto cleanup;
        }

        F->err = lis_vector_duplicate(
                F->A,
                &F->zero_imag);

        if( F->err ) goto cleanup;

        F->err = lis_vector_duplicate(
                F->A,
                &F->q_real);

        if( F->err ) goto cleanup;

        F->err = lis_vector_duplicate(
                F->A,
                &F->q_imag);

        if( F->err ) goto cleanup;

        F->err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                F->zero_imag);

        if( F->err ) goto cleanup;

        for(F->i=0;F->i<F->block_size;F->i++)
        {
                F->err = lis_vector_duplicate(
                        F->A,
                        &F->W[F->i]);

                if( F->err ) goto cleanup;

                if( F->B!=NULL )
                {
                        F->err = lis_vector_duplicate(
                                F->A,
                                &F->BW[F->i]);

                        if( F->err ) goto cleanup;
                }
        }

          
                                               
          
                                                                  
                                          
          
                           
           
        for(F->i=0;F->i<F->block_size;F->i++)
        {
                F->pc=1;
K->projector.A = F->A;
K->projector.B = F->B;
K->projector.emin = F->emin;
K->projector.emax = F->emax;
K->projector.ratio = F->ratio;
K->projector.npoints = F->npoints;
K->projector.x_real = F->X[F->i];
K->projector.x_imag = F->zero_imag;
K->projector.q_real = F->q_real;
K->projector.q_imag = F->q_imag;
K->projector.solver = F->solver;
K->projector.pc=0; K->projector.parent=4;
K->active=1; return 3;
resume_1: F->pc=0; F->err=K->completion;

                if( F->err ) goto cleanup;

                F->err = lis_vector_copy(
                        F->q_real,
                        F->W[F->i]);

                if( F->err ) goto cleanup;

        }

          
                                                         
          
                                                                    
                                                                     
           
        if( F->B!=NULL )
        {
                F->pc=2;
K->brank.B = F->B;
K->brank.capacity = F->block_size;
K->brank.nominal_rank = F->block_size;
K->brank.W = F->W;
K->brank.BW = F->BW;
K->brank.Q = F->Q;
K->brank.BQ = F->BQ;
K->brank.rank_tol = F->rank_tol;
K->brank.active_rank = F->active_rank;
K->brank.gram = F->gram;
K->brank.eval = F->eval;
K->brank.evec = F->evec;
K->brank.transform = F->transform;
K->brank.pc=0; K->brank.parent=4;
K->active=3; return 3;
resume_2: F->pc=0; F->err=K->completion;
        }
        else
        {
                F->err = lis_eigen_block_rank_deflate(
                        F->block_size,
                        F->block_size,
                        F->W,
                        F->Q,
                        F->rank_tol,
                        F->active_rank,
                        F->gram,
                        F->eval,
                        F->evec,
                        F->transform);
        }

cleanup:

        if( F->q_imag!=NULL )
        {
                lis_vector_destroy(
                        F->q_imag);
        }

        if( F->q_real!=NULL )
        {
                lis_vector_destroy(
                        F->q_real);
        }

        if( F->zero_imag!=NULL )
        {
                lis_vector_destroy(
                        F->zero_imag);
        }

        lis_efeast_vector_block_destroy(
                F->block_size,
                F->BW);

        lis_efeast_vector_block_destroy(
                F->block_size,
                F->W);

        if( F->transform!=NULL )
        {
                lis_free(
                        F->transform);
        }

        if( F->evec!=NULL )
        {
                lis_free(
                        F->evec);
        }

        if( F->eval!=NULL )
        {
                lis_free(
                        F->eval);
        }

        if( F->gram!=NULL )
        {
                lis_free(
                        F->gram);
        }

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_efeast_projected_rayleigh_ritz"
static int nf_step_rr(lis_rci_native_feast *K){
struct nf_rr *F=&K->rr;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}

















        LIS_DEBUG_FUNC_IN;

        F->AQ = NULL;

        F->projected = NULL;
        F->coeff = NULL;

        if( F->projected_niter!=NULL )
        {
                *F->projected_niter = 0;
        }

        if( F->projected_error!=NULL )
        {
                *F->projected_error = 0.0;
        }

        if( F->A==NULL ||
            F->active_rank<=0 ||
            F->Q==NULL ||
            F->ritz_value==NULL ||
            F->ritz_vector==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST projected Rayleigh-Ritz received "
                        "invalid arguments\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        for(F->i=0;F->i<F->active_rank;F->i++)
        {
                if( F->Q[F->i]==NULL ||
                    F->ritz_vector[F->i]==NULL ||
                    !lis_efeast_shifted_vector_compatible(
                            F->A,
                            F->Q[F->i]) ||
                    !lis_efeast_shifted_vector_compatible(
                            F->A,
                            F->ritz_vector[F->i]) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST projected Rayleigh-Ritz contains "
                                "an incompatible vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                }

                for(F->j=0;F->j<F->active_rank;F->j++)
                {
                        if( F->ritz_vector[F->i]==F->Q[F->j] )
                        {
                                LIS_SETERR(
                                        LIS_ERR_ILL_ARG,
                                        "FEAST Ritz reconstruction requires "
                                        "out-of-place vectors\n");

                                LIS_DEBUG_FUNC_OUT;
                                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                        }
                }

                for(F->j=0;F->j<F->i;F->j++)
                {
                        if( F->ritz_vector[F->i]==F->ritz_vector[F->j] )
                        {
                                LIS_SETERR(
                                        LIS_ERR_ILL_ARG,
                                        "FEAST Ritz vectors must be distinct\n");

                                LIS_DEBUG_FUNC_OUT;
                                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                        }
                }
        }

        F->vector_bytes =
                (size_t)F->active_rank
                *
                sizeof(LIS_VECTOR);

        F->square_count =
                (size_t)F->active_rank
                *
                (size_t)F->active_rank;

        F->square_bytes =
                F->square_count
                *
                sizeof(LIS_SCALAR);

        F->AQ = (LIS_VECTOR *)lis_malloc(
                F->vector_bytes,
                "lis_efeast_projected_rayleigh_ritz::AQ");

        if( F->AQ==NULL )
        {
                LIS_SETERR_MEM(
                        F->vector_bytes);

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_OUT_OF_MEMORY);
        }

        memset(
                F->AQ,
                0,
                F->vector_bytes);

        F->projected = (LIS_SCALAR *)lis_malloc(
                F->square_bytes,
                "lis_efeast_projected_rayleigh_ritz::projected");

        if( F->projected==NULL )
        {
                F->err = LIS_OUT_OF_MEMORY;

                LIS_SETERR_MEM(
                        F->square_bytes);

                goto cleanup;
        }

        F->coeff = (LIS_SCALAR *)lis_malloc(
                F->square_bytes,
                "lis_efeast_projected_rayleigh_ritz::coeff");

        if( F->coeff==NULL )
        {
                F->err = LIS_OUT_OF_MEMORY;

                LIS_SETERR_MEM(
                        F->square_bytes);

                goto cleanup;
        }

        for(F->i=0;F->i<F->active_rank;F->i++)
        {
                F->err = lis_vector_duplicate(
                        F->A,
                        &F->AQ[F->i]);

                if( F->err ) goto cleanup;
        }

          
                                                                  
                                                              
           
        F->pc=1;
return nf_block(K,NF_A,F->A,F->active_rank,F->Q,F->AQ);
resume_1: F->pc=0; F->err=K->completion;

        if( F->err ) goto cleanup;

          
                       
           
        F->err = lis_eigen_block_gram(
                F->active_rank,
                F->Q,
                F->active_rank,
                F->AQ,
                F->projected);

        if( F->err ) goto cleanup;

          
                                                                   
                                                             
          
                                                            
                                                               
                                                      
           
        for(F->i=0;F->i<F->active_rank;F->i++)
        {

                for(F->j=F->i+1;F->j<F->active_rank;F->j++)
                {

                        F->value =
                                (LIS_SCALAR)0.5
                                *
                                (
                                        F->projected[F->i*F->active_rank+F->j]
                                        +
                                        F->projected[F->j*F->active_rank+F->i]
                                );

                        F->projected[F->i*F->active_rank+F->j] =
                                F->value;

                        F->projected[F->j*F->active_rank+F->i] =
                                F->value;

                }
        }

          
                                                                  
                                                                  
           
        F->err = lis_eigen_projected_solve(
                F->active_rank,
                F->projected,
                F->ritz_value,
                F->coeff,
                F->projected_niter,
                F->projected_error);

        if( F->err ) goto cleanup;

          
                           
          
                                                   
           
        F->err = lis_eigen_ritz_reconstruct(
                F->active_rank,
                F->Q,
                F->active_rank,
                F->coeff,
                F->active_rank,
                F->ritz_vector);

        if( F->err ) goto cleanup;

          
                                                                   
                                             
           
        for(F->k=0;F->k<F->active_rank;F->k++)
        {

                if( !isfinite((double)F->ritz_value[F->k]) )

                {
                        LIS_SETERR(
                                LIS_BREAKDOWN,
                                "FEAST projected Rayleigh-Ritz produced "
                                "a non-finite eigenvalue\n");

                        F->err = LIS_BREAKDOWN;
                        goto cleanup;
                }
        }

        F->err = LIS_SUCCESS;

cleanup:

        lis_efeast_vector_block_destroy(
                F->active_rank,
                F->AQ);

        if( F->coeff!=NULL )
        {
                lis_free(
                        F->coeff);
        }

        if( F->projected!=NULL )
        {
                lis_free(
                        F->projected);
        }

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_esolver_standard_residual"
static int nf_step_standard(lis_rci_native_feast *K){
struct nf_standard *F=&K->standard;
switch(F->pc){case 0:break;case 1:goto resume_1;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}





        LIS_DEBUG_FUNC_IN;

        F->pc=1;
return nf_one(K,NF_A,F->A,F->x,F->Ax);
resume_1: F->pc=0; F->err=K->completion;
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->err = lis_vector_nrm2(F->Ax,&F->anrm);
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->err = lis_vector_nrm2(F->x,&F->xnrm);
        if( F->err ) return nf_return(K,F->parent,F->err);

        lis_vector_axpyz(-F->lambda,F->x,F->Ax,F->r);

        F->err = lis_vector_nrm2(F->r,&F->rnrm);
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->lmag = fabs(F->lambda);

        if( !isfinite((double)F->anrm) ||
            !isfinite((double)F->xnrm) ||
            !isfinite((double)F->rnrm) ||
            !isfinite((double)F->lmag) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite standard eigenpair certification data\n");
                return nf_return(K,F->parent,LIS_BREAKDOWN);
        }

        F->denom = F->anrm+F->lmag*F->xnrm;
        if( !isfinite((double)F->denom) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite standard eigenpair certification denominator\n");
                return nf_return(K,F->parent,LIS_BREAKDOWN);
        }

        if( F->denom>0.0 )
        {
                *F->resid = F->rnrm/F->denom;
        }
        else
        {
                *F->resid = F->rnrm;
        }

        if( !isfinite((double)*F->resid) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite standard eigenpair certification residual\n");
                return nf_return(K,F->parent,LIS_BREAKDOWN);
        }

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_esolver_generalized_residual"
static int nf_step_generalized(lis_rci_native_feast *K){
struct nf_generalized *F=&K->generalized;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}





        LIS_DEBUG_FUNC_IN;

        F->pc=1;
return nf_one(K,NF_A,F->A,F->x,F->Ax);
resume_1: F->pc=0; F->err=K->completion;
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->pc=2;
return nf_one(K,NF_B,F->B,F->x,F->Bx);
resume_2: F->pc=0; F->err=K->completion;
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->err = lis_vector_nrm2(F->Ax,&F->anrm);
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->err = lis_vector_nrm2(F->Bx,&F->bnrm);
        if( F->err ) return nf_return(K,F->parent,F->err);

        lis_vector_axpyz(-F->lambda,F->Bx,F->Ax,F->r);

        F->err = lis_vector_nrm2(F->r,&F->rnrm);
        if( F->err ) return nf_return(K,F->parent,F->err);

        F->lmag = fabs(F->lambda);

        if( !isfinite((double)F->anrm) ||
            !isfinite((double)F->bnrm) ||
            !isfinite((double)F->rnrm) ||
            !isfinite((double)F->lmag) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite generalized eigenpair certification data\n");
                return nf_return(K,F->parent,LIS_BREAKDOWN);
        }

        F->denom = F->anrm+F->lmag*F->bnrm;
        if( !isfinite((double)F->denom) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite generalized eigenpair certification denominator\n");
                return nf_return(K,F->parent,LIS_BREAKDOWN);
        }

        if( F->denom>0.0 )
        {
                *F->resid = F->rnrm/F->denom;
        }
        else
        {
                *F->resid = F->rnrm;
        }

        if( !isfinite((double)*F->resid) )
        {
                LIS_SETERR(LIS_BREAKDOWN,
                        "non-finite generalized eigenpair certification residual\n");
                return nf_return(K,F->parent,LIS_BREAKDOWN);
        }

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,LIS_SUCCESS);

}
#undef __FUNC__
#define __FUNC__ "lis_efeast_certify_ritz_pairs"
static int nf_step_cert(lis_rci_native_feast *K){
struct nf_cert *F=&K->cert;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}























        LIS_DEBUG_FUNC_IN;

        F->Ax = NULL;
        F->Bx = NULL;
        F->r = NULL;

        if( F->selected_count!=NULL )
        {
                *F->selected_count = 0;
        }

        if( F->trace!=NULL )
        {
                *F->trace = 0.0;
        }

        if( F->trace_error!=NULL )
        {
                *F->trace_error = 0.0;
        }

        if( F->max_residual!=NULL )
        {
                *F->max_residual = 0.0;
        }

        if( F->converged!=NULL )
        {
                *F->converged = LIS_FALSE;
        }

        if( F->A==NULL ||
            !(F->emin<F->emax) ||
            !(F->tol>0.0) ||
            !isfinite((double)F->emin) ||
            !isfinite((double)F->emax) ||
            !isfinite((double)F->tol) ||
            F->active_rank<=0 ||
            F->ritz_value==NULL ||
            F->ritz_vector==NULL ||
            F->selected_index==NULL ||
            F->residual==NULL ||
            F->selected_count==NULL ||
            F->trace==NULL ||
            F->trace_error==NULL ||
            F->max_residual==NULL ||
            F->converged==NULL ||
            F->previous_count < -1 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST Ritz certification received invalid arguments\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        if( F->previous_count>=0 &&
            !isfinite((double)F->previous_trace) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "FEAST previous trace must be finite\n");

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
        }

        for(F->i=0;F->i<F->active_rank;F->i++)
        {
                F->selected_index[F->i] = -1;
                F->residual[F->i] = 0.0;

                if( F->ritz_vector[F->i]==NULL ||
                    !lis_efeast_shifted_vector_compatible(
                            F->A,
                            F->ritz_vector[F->i]) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST Ritz certification contains "
                                "an incompatible Ritz vector\n");

                        LIS_DEBUG_FUNC_OUT;
                        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                }

                if( F->B!=NULL &&
                    !lis_efeast_shifted_vector_compatible(
                            F->B,
                            F->ritz_vector[F->i]) )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "FEAST generalized Ritz certification "
                                "contains an incompatible B layout\n");

                        LIS_DEBUG_FUNC_OUT;
                        return nf_return(K,F->parent,LIS_ERR_ILL_ARG);
                }
        }

        F->err = lis_vector_duplicate(
                F->A,
                &F->Ax);

        if( F->err ) goto cleanup;

        F->err = lis_vector_duplicate(
                F->A,
                &F->r);

        if( F->err ) goto cleanup;

        if( F->B!=NULL )
        {
                F->err = lis_vector_duplicate(
                        F->A,
                        &F->Bx);

                if( F->err ) goto cleanup;
        }

        F->count = 0;
        F->local_trace = 0.0;
        F->local_trace_error = 0.0;
        F->local_max_residual = 0.0;

        for(F->i=0;F->i<F->active_rank;F->i++)
        {

                F->lambda_real =
                        (LIS_REAL)F->ritz_value[F->i];

                if( !isfinite((double)F->lambda_real) )
                {
                        LIS_SETERR(
                                LIS_BREAKDOWN,
                                "FEAST Ritz certification received "
                                "a non-finite Ritz value\n");

                        F->err = LIS_BREAKDOWN;
                        goto cleanup;
                }


                  
                                                            
                  
                                                                  
                                                                      
                                                               
                                                                   
                                                                     
                                                                     
                                     
                   
                F->interval_scale = 1.0;

                if( lis_feast_real_abs((LIS_REAL)F->emin) >
                    (double)F->interval_scale )
                {
                        F->interval_scale =
                                (LIS_REAL)lis_feast_real_abs((LIS_REAL)F->emin);
                }

                if( lis_feast_real_abs((LIS_REAL)F->emax) >
                    (double)F->interval_scale )
                {
                        F->interval_scale =
                                (LIS_REAL)lis_feast_real_abs((LIS_REAL)F->emax);
                }

                if( lis_feast_real_abs((LIS_REAL)F->lambda_real) >
                    (double)F->interval_scale )
                {
                        F->interval_scale =
                                (LIS_REAL)lis_feast_real_abs((LIS_REAL)F->lambda_real);
                }

                F->interval_slack =
                        (LIS_REAL)(
                                64.0
                                *
                                DBL_EPSILON)
                        *
                        F->interval_scale;

                if( F->lambda_real <
                        F->emin-F->interval_slack ||
                    F->lambda_real >
                        F->emax+F->interval_slack )
                {
                        continue;
                }

                if( F->B==NULL )
                {
                        F->pc=1;
K->standard.A = F->A;
K->standard.x = F->ritz_vector[F->i];
K->standard.lambda = F->ritz_value[F->i];
K->standard.Ax = F->Ax;
K->standard.r = F->r;
K->standard.resid = &F->pair_residual;
K->standard.pc=0; K->standard.parent=8;
K->active=6; return 3;
resume_1: F->pc=0; F->err=K->completion;
                }
                else
                {
                        F->pc=2;
K->generalized.A = F->A;
K->generalized.B = F->B;
K->generalized.x = F->ritz_vector[F->i];
K->generalized.lambda = F->ritz_value[F->i];
K->generalized.Ax = F->Ax;
K->generalized.Bx = F->Bx;
K->generalized.r = F->r;
K->generalized.resid = &F->pair_residual;
K->generalized.pc=0; K->generalized.parent=8;
K->active=7; return 3;
resume_2: F->pc=0; F->err=K->completion;
                }

                if( F->err ) goto cleanup;

                if( !isfinite((double)F->pair_residual) )
                {
                        LIS_SETERR(
                                LIS_BREAKDOWN,
                                "FEAST Ritz certification produced "
                                "a non-finite physical residual\n");

                        F->err = LIS_BREAKDOWN;
                        goto cleanup;
                }

                F->selected_index[F->count] = F->i;
                F->residual[F->count] = F->pair_residual;

                F->local_trace +=
                        F->lambda_real;

                if( F->pair_residual >
                    F->local_max_residual )
                {
                        F->local_max_residual =
                                F->pair_residual;
                }

                F->count++;
        }

        if( !isfinite((double)F->local_trace) ||
            !isfinite((double)F->local_max_residual) )
        {
                LIS_SETERR(
                        LIS_BREAKDOWN,
                        "FEAST Ritz certification produced "
                        "non-finite convergence data\n");

                F->err = LIS_BREAKDOWN;
                goto cleanup;
        }

          
                                 
          
                            
                                      
           
        if( F->previous_count>=0 )
        {
                F->trace_scale = 1.0;

                if( fabs((double)F->local_trace) >
                    (double)F->trace_scale )
                {
                        F->trace_scale =
                                fabs(
                                        (double)F->local_trace);
                }

                if( fabs((double)F->previous_trace) >
                    (double)F->trace_scale )
                {
                        F->trace_scale =
                                fabs(
                                        (double)F->previous_trace);
                }

                F->local_trace_error =
                        fabs(
                                (double)(
                                        F->local_trace
                                        -
                                        F->previous_trace))
                        /
                        F->trace_scale;

                if( !isfinite(
                            (double)F->local_trace_error) )
                {
                        LIS_SETERR(
                                LIS_BREAKDOWN,
                                "FEAST trace convergence metric "
                                "is non-finite\n");

                        F->err = LIS_BREAKDOWN;
                        goto cleanup;
                }
        }

        *F->selected_count = F->count;
        *F->trace = F->local_trace;
        *F->trace_error = F->local_trace_error;
        *F->max_residual = F->local_max_residual;

          
                                                                   
                                                                  
                                                           
           
        if( F->previous_count>=0 &&
            F->previous_count==F->count &&
            F->local_trace_error<=F->tol &&
            F->local_max_residual<=F->tol )
        {
                *F->converged = LIS_TRUE;
        }
        else
        {
                *F->converged = LIS_FALSE;
        }

        F->err = LIS_SUCCESS;

cleanup:

        if( F->Bx!=NULL )
        {
                lis_vector_destroy(
                        F->Bx);
        }

        if( F->r!=NULL )
        {
                lis_vector_destroy(
                        F->r);
        }

        if( F->Ax!=NULL )
        {
                lis_vector_destroy(
                        F->Ax);
        }

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,F->err);

}
#undef __FUNC__
#define __FUNC__ "lis_efeast"
static int nf_step_feast(lis_rci_native_feast *K){
struct nf_feast *F=&K->feast;
switch(F->pc){case 0:break;case 1:goto resume_1;case 2:goto resume_2;case 3:goto resume_3;default:return nf_return(K,F->parent,LIS_BREAKDOWN);}















































        LIS_DEBUG_FUNC_IN;

        F->shifted_solver = NULL;

        F->ritz_value = NULL;
        F->selected_index = NULL;
        F->selected_residual = NULL;

        F->err = lis_efeast_check_params(
                F->esolver);

        if( F->err )
        {
                if( F->esolver!=NULL )
                {
                          
                                                                     
                                                       
                           
                        F->esolver->nev =
                                0;

                        F->esolver->retcode =
                                F->err;
                }

                LIS_DEBUG_FUNC_OUT;
                return nf_return(K,F->parent,F->err);
        }

        F->A = F->esolver->A;
        F->B = F->esolver->B;

        F->ss = F->esolver->options[
                LIS_EOPTIONS_SUBSPACE];

        F->emaxiter = F->esolver->options[
                LIS_EOPTIONS_MAXITER];

        F->tol = F->esolver->params[
                LIS_EPARAMS_RESID
                -
                LIS_EOPTIONS_LEN];

        if( F->A==NULL ||
            F->esolver->x==NULL ||
            F->esolver->evalue==NULL ||
            F->esolver->resid==NULL ||
            F->esolver->iter==NULL ||
            F->esolver->iter2==NULL ||
            F->esolver->rhistory==NULL ||
            F->esolver->evector==NULL ||
            F->esolver->work==NULL ||
            F->esolver->worklen
                !=LIS_EFEAST_WORK_BLOCKS*F->ss )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST generic eigensolver state "
                        "is incomplete\n");

                F->err = LIS_ERR_ILL_ARG;
                goto cleanup;
        }

        if( F->ss<=0 ||
            F->ss>F->A->gn ||
            F->emaxiter<=0 ||
            !isfinite((double)F->tol) ||
            F->tol<=0.0 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST outer-loop configuration "
                        "is invalid\n");

                F->err = LIS_ERR_ILL_ARG;
                goto cleanup;
        }

        if( F->B!=NULL &&
            (F->B->n!=F->A->n ||
             F->B->gn!=F->A->gn) )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "Native FEAST generalized matrices "
                        "have incompatible dimensions\n");

                F->err = LIS_ERR_ILL_ARG;
                goto cleanup;
        }

          
                                                           
           
        F->X =
                F->esolver->work;

        F->Q =
                F->esolver->work
                +
                F->ss;

        F->BQ =
                F->esolver->work
                +
                2*F->ss;

        F->Ritz =
                F->esolver->work
                +
                3*F->ss;

        F->scalar_bytes =
                (size_t)F->ss
                *
                sizeof(LIS_SCALAR);

        F->integer_bytes =
                (size_t)F->ss
                *
                sizeof(LIS_INT);

        F->real_bytes =
                (size_t)F->ss
                *
                sizeof(LIS_REAL);

        F->ritz_value =
                (LIS_SCALAR *)lis_malloc(
                        F->scalar_bytes,
                        "lis_efeast::ritz_value");

        if( F->ritz_value==NULL )
        {
                LIS_SETERR_MEM(
                        F->scalar_bytes);

                F->err = LIS_OUT_OF_MEMORY;
                goto cleanup;
        }

        F->selected_index =
                (LIS_INT *)lis_malloc(
                        F->integer_bytes,
                        "lis_efeast::selected_index");

        if( F->selected_index==NULL )
        {
                LIS_SETERR_MEM(
                        F->integer_bytes);

                F->err = LIS_OUT_OF_MEMORY;
                goto cleanup;
        }

        F->selected_residual =
                (LIS_REAL *)lis_malloc(
                        F->real_bytes,
                        "lis_efeast::selected_residual");

        if( F->selected_residual==NULL )
        {
                LIS_SETERR_MEM(
                        F->real_bytes);

                F->err = LIS_OUT_OF_MEMORY;
                goto cleanup;
        }

        memset(
                F->ritz_value,
                0,
                F->scalar_bytes);

        memset(
                F->selected_index,
                0,
                F->integer_bytes);

        memset(
                F->selected_residual,
                0,
                F->real_bytes);

          
                                                                  
                                                                     
                                                                   
                                                          
           
        for(F->i=0;F->i<F->ss+2;F->i++)
        {
                F->esolver->evalue[F->i] =
                        0.0;

                F->esolver->resid[F->i] =
                        0.0;

                F->esolver->iter[F->i] =
                        0;

                F->esolver->iter2[F->i] =
                        0;
        }

        F->esolver->nev = 0;
        F->esolver->tol = F->tol;

        F->esolver->rhistory[0] =
                1.0;

        F->err = lis_efeast_seed_subspace(
                F->esolver,
                F->ss,
                F->X);

        if( F->err )
                goto cleanup;

        F->err = lis_efeast_create_shifted_solver(
                F->tol,
                F->esolver->feast_solver,
                &F->shifted_solver);

        if( F->err )
                goto cleanup;

        F->previous_count = -1;
        F->previous_trace = 0.0;

          
                                                                  
                                                                    
                                                               
          
                                                   
           
        F->rank_tol =
                (LIS_REAL)1.0e-8;

        F->err = LIS_MAXITER;

        for(F->iter=1;
            F->iter<=F->emaxiter;
            F->iter++)
        {
                F->active_rank = 0;

                F->projected_niter = 0;
                F->projected_error = 0.0;

                F->selected_count = 0;
                F->trace = 0.0;
                F->trace_error = 0.0;
                F->max_residual = 0.0;
                F->converged = LIS_FALSE;

                memset(
                        F->ritz_value,
                        0,
                        F->scalar_bytes);

                memset(
                        F->selected_index,
                        0,
                        F->integer_bytes);

                memset(
                        F->selected_residual,
                        0,
                        F->real_bytes);

                  
                    
                                        
                                                   
                   
                F->pc=1;
K->filter.A = F->A;
K->filter.B = F->B;
K->filter.emin = F->esolver->feast_emin;
K->filter.emax = F->esolver->feast_emax;
K->filter.ratio = F->esolver->feast_ratio;
K->filter.npoints = F->esolver->feast_points;
K->filter.block_size = F->ss;
K->filter.X = F->X;
K->filter.Q = F->Q;
K->filter.BQ = F->B!=NULL ? F->BQ : NULL;
K->filter.rank_tol = F->rank_tol;
K->filter.active_rank = &F->active_rank;
K->filter.solver = F->shifted_solver;
K->filter.pc=0; K->filter.parent=9;
K->active=4; return 3;
resume_1: F->pc=0; F->err=K->completion;

                if( F->err )
                        goto cleanup;

                if( F->active_rank>0 )
                {
                          
                            
                                                              
                                                     
                           
                        F->pc=2;
K->rr.A = F->A;
K->rr.active_rank = F->active_rank;
K->rr.Q = F->Q;
K->rr.ritz_value = F->ritz_value;
K->rr.ritz_vector = F->Ritz;
K->rr.projected_niter = &F->projected_niter;
K->rr.projected_error = &F->projected_error;
K->rr.pc=0; K->rr.parent=9;
K->active=5; return 3;
resume_2: F->pc=0; F->err=K->completion;

                        if( F->err )
                                goto cleanup;

                          
                                                                     
                                                                     
                                                   
                           
                        F->pc=3;
K->cert.A = F->A;
K->cert.B = F->B;
K->cert.emin = F->esolver->feast_emin;
K->cert.emax = F->esolver->feast_emax;
K->cert.tol = F->tol;
K->cert.active_rank = F->active_rank;
K->cert.ritz_value = F->ritz_value;
K->cert.ritz_vector = F->Ritz;
K->cert.previous_count = F->previous_count;
K->cert.previous_trace = F->previous_trace;
K->cert.selected_index = F->selected_index;
K->cert.residual = F->selected_residual;
K->cert.selected_count = &F->selected_count;
K->cert.trace = &F->trace;
K->cert.trace_error = &F->trace_error;
K->cert.max_residual = &F->max_residual;
K->cert.converged = &F->converged;
K->cert.pc=0; K->cert.parent=9;
K->active=8; return 3;
resume_3: F->pc=0; F->err=K->completion;

                        if( F->err )
                                goto cleanup;
                }
                else
                {
                          
                                                                  
                                                                      
                                                                    
                                                                      
                           
                        F->selected_count = 0;
                        F->trace = 0.0;
                        F->max_residual = 0.0;

                        if( F->previous_count>=0 )
                        {
                                F->trace_error =
                                        lis_feast_real_abs((LIS_REAL)F->previous_trace);

                                if( F->trace_error<=F->tol &&
                                    F->previous_count==0 )
                                {
                                        F->converged =
                                                LIS_TRUE;
                                }
                        }
                }

                F->history_value =
                        F->trace_error;

                if( F->max_residual >
                    F->history_value )
                {
                        F->history_value =
                                F->max_residual;
                }

                F->esolver->rhistory[F->iter] =
                        F->history_value;

                  
                                                               
                                                                   
                                                             
                  
                                                                  
                                                                 
                                        
                   
                if( F->converged )
                {
                        for(F->i=0;
                            F->i<F->selected_count;
                            F->i++)
                        {
                                F->first =
                                        F->selected_index[F->i];

                                if( F->first<0 ||
                                    F->first>=F->active_rank )
                                {
                                        LIS_SETERR(
                                                LIS_BREAKDOWN,
                                                "Native FEAST selected "
                                                "an invalid Ritz index\n");

                                        F->err = LIS_BREAKDOWN;
                                        goto cleanup;
                                }

                                F->esolver->evalue[F->i] =
                                        F->ritz_value[F->first];

                                F->esolver->resid[F->i] =
                                        F->selected_residual[F->i];

                                F->esolver->iter[F->i] =
                                        F->iter;

                                F->esolver->iter2[F->i] =
                                        F->projected_niter;

                                F->err = lis_vector_duplicate(
                                        F->A,
                                        &F->esolver->evector[F->i]);

                                if( F->err )
                                        goto cleanup;

                                F->err = lis_vector_copy(
                                        F->Ritz[F->first],
                                        F->esolver->evector[F->i]);

                                if( F->err )
                                        goto cleanup;
                        }

                          
                                                                    
                                                                  
                                                                     
                                        
                           
                        if( F->selected_count>0 )
                        {
                                F->err = lis_vector_copy(
                                        F->esolver->evector[0],
                                        F->esolver->x);

                                if( F->err )
                                        goto cleanup;
                        }
                        else
                        {
                                F->err = lis_vector_set_all(
                                        0.0,
                                        F->esolver->x);

                                if( F->err )
                                        goto cleanup;
                        }

                          
                                                                    
                                                                     
                                                                   
                                                
                          
                                                                
                                                                        
                                                                     
                           
                        F->esolver->nev =
                                F->selected_count;

                        F->err = LIS_SUCCESS;
                        goto cleanup;
                }

                  
                                                                  
                                                                   
                               
                   
                if( F->selected_count>0 )
                {
                        F->first =
                                F->selected_index[0];

                        if( F->first<0 ||
                            F->first>=F->active_rank )
                        {
                                LIS_SETERR(
                                        LIS_BREAKDOWN,
                                        "Native FEAST selected an "
                                        "invalid Ritz index\n");

                                F->err = LIS_BREAKDOWN;
                                goto cleanup;
                        }

                        F->esolver->evalue[0] =
                                F->ritz_value[F->first];

                        F->esolver->resid[0] =
                                F->selected_residual[0];

                        F->err = lis_vector_copy(
                                F->Ritz[F->first],
                                F->esolver->x);

                        if( F->err )
                                goto cleanup;
                }
                else
                {
                        F->esolver->evalue[0] =
                                0.0;

                        F->esolver->resid[0] =
                                0.0;
                }

                F->esolver->iter[0] =
                        F->iter;

                F->esolver->iter2[0] =
                        F->projected_niter;

                  
                                            
                  
                                                                     
                                                                   
                                                                   
                   
                for(F->i=0;F->i<F->active_rank;F->i++)
                {
                        F->err = lis_vector_copy(
                                F->Ritz[F->i],
                                F->X[F->i]);

                        if( F->err )
                                goto cleanup;
                }

                F->previous_count =
                        F->selected_count;

                F->previous_trace =
                        F->trace;
        }

        LIS_SETERR(
                LIS_MAXITER,
                "Native FEAST reached the maximum "
                "number of outer iterations\n");

        F->err = LIS_MAXITER;

cleanup:

        if( F->err!=LIS_SUCCESS )
        {
                  
                                                                
                                                                       
                                           
                   
                F->esolver->nev =
                        0;
        }

        if( F->shifted_solver!=NULL )
        {
                lis_solver_destroy(
                        F->shifted_solver);
        }

        if( F->selected_residual!=NULL )
        {
                lis_free(
                        F->selected_residual);
        }

        if( F->selected_index!=NULL )
        {
                lis_free(
                        F->selected_index);
        }

        if( F->ritz_value!=NULL )
        {
                lis_free(
                        F->ritz_value);
        }

        F->esolver->retcode =
                F->err;

        LIS_DEBUG_FUNC_OUT;
        return nf_return(K,F->parent,F->err);

}
static int nf_step(lis_rci_native_feast *K){switch(K->active){
case 1:return nf_step_projector(K);
case 2:return nf_step_bgram(K);
case 3:return nf_step_brank(K);
case 4:return nf_step_filter(K);
case 5:return nf_step_rr(K);
case 6:return nf_step_standard(K);
case 7:return nf_step_generalized(K);
case 8:return nf_step_cert(K);
case 9:return nf_step_feast(K);
default:return nf_return(K,0,LIS_BREAKDOWN);}}
