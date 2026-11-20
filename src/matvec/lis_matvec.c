/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
	#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
	#include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <stdint.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

LIS_MATVEC_FUNC LIS_MATVEC  = lis_matvec;
LIS_MATVEC_FUNC LIS_MATVECH = lis_matvech;


/*
 * Generic serial USER kernels backed by LIS_MATRIX_USER_GET_ROW.
 *
 * Explicit application callbacks remain the fast path.  These helpers
 * let applications exposing complete rows avoid duplicating basic
 * matrix-vector kernels.
 *
 * LIS_PRECON_ROW_SOURCE owns only a reusable maximum-row-sized buffer;
 * it never materializes storage proportional to total matrix nnz.
 */
static LIS_INT
lis_matvec_user_rows(
	LIS_MATRIX A,
	const LIS_SCALAR *x,
	LIS_SCALAR *y)
{
#ifdef _OPENMP
	if( A->user_get_row_thread_safe && omp_get_max_threads()>1 )
	{
		LIS_INT parallel_err = LIS_SUCCESS;

#pragma omp parallel
		{
			LIS_PRECON_ROW_SOURCE row;
			LIS_INT local_err;
			LIS_INT row_ready;
			LIS_INT i,k;

			local_err = lis_precon_row_source_init_mode(A,&row,LIS_TRUE);
			row_ready = (local_err==LIS_SUCCESS);
			if( local_err )
			{
#pragma omp critical(lis_user_row_error)
				{
					if( parallel_err==LIS_SUCCESS ) parallel_err = local_err;
				}
			}

			/* Every thread in the team must encounter this worksharing region. */
#pragma omp for schedule(static)
			for(i=0;i<A->n;i++)
			{
				LIS_SCALAR sum = (LIS_SCALAR)0.0;
				if( !row_ready ) continue;
				local_err = lis_precon_row_source_get(&row,i);
				if( local_err )
				{
#pragma omp critical(lis_user_row_error)
					{
						if( parallel_err==LIS_SUCCESS ) parallel_err = local_err;
					}
					continue;
				}
				for(k=0;k<row.nnz;k++)
					sum += row.value[k] * x[row.index[k]];
				y[i] = sum;
			}
			if( row_ready ) lis_precon_row_source_destroy(&row);
		}
		return parallel_err;
	}
#endif

	{
		LIS_PRECON_ROW_SOURCE row;
		LIS_INT err;
		LIS_INT i,k;
		LIS_SCALAR sum;

		err = lis_precon_row_source_init(A,&row);
		if( err ) return err;
		for(i=0;i<A->n;i++)
		{
			err = lis_precon_row_source_get(&row,i);
			if( err )
			{
				lis_precon_row_source_destroy(&row);
				return err;
			}
			sum = (LIS_SCALAR)0.0;
			for(k=0;k<row.nnz;k++) sum += row.value[k] * x[row.index[k]];
			y[i] = sum;
		}
		lis_precon_row_source_destroy(&row);
		return LIS_SUCCESS;
	}
}

static LIS_INT
lis_matvech_user_rows(
	LIS_MATRIX A,
	const LIS_SCALAR *x,
	LIS_SCALAR *y)
{
	LIS_PRECON_ROW_SOURCE row;
	LIS_INT err;
	LIS_INT i,k;

	err = lis_precon_row_source_init(A,&row);
	if( err ) return err;

	for(i=0;i<A->n;i++)
	{
		y[i] = (LIS_SCALAR)0.0;
	}

	for(i=0;i<A->n;i++)
	{
		err = lis_precon_row_source_get(&row,i);
		if( err )
		{
			lis_precon_row_source_destroy(&row);
			return err;
		}

		for(k=0;k<row.nnz;k++)
		{
			y[row.index[k]] +=
				conj(row.value[k])
				*
				x[i];
		}
	}

	lis_precon_row_source_destroy(&row);

	return LIS_SUCCESS;
}


#ifdef USE_MPI
/* Correctness-first distributed USER row fallback.  Full global vectors
 * replace a sparse halo here; memory and communication are O(global_n)
 * per rank.  All ranks must enter a USER operation collectively. */
/*
 * Sparse, per-call row-backed USER exchange. We deliberately do not retain
 * a halo cache: rows can change without a structure-generation notification.
 *
 * Phase 1: collect/validate local rows and classify remote columns.
 * Phase 2: all-to-all request counts and global column indices.
 * Phase 3: forward Ax fetches only requested x entries; transpose A^H x
 *          returns contributions along those same reversed request edges.
 *
 * Repeated off-rank columns are not yet deduplicated, so communication is
 * bounded by off-rank NNZ rather than the number of unique remote columns.
 * Count metadata still costs O(nprocs) per call.
 */
static int
lis_user_sparse_owner(const LIS_INT *ranges,int np,LIS_INT col)
{
    int lo=0,hi=np;
    while(lo<hi)
    {
        int mid=lo+(hi-lo)/2;
        if(col>=ranges[mid+1])lo=mid+1;
        else hi=mid;
    }
    return lo<np && col>=ranges[lo] ? lo : -1;
}

static LIS_INT
lis_user_sparse_sync_error(LIS_MATRIX A,LIS_INT local)
{
    LIS_INT global=LIS_SUCCESS;
    int rc=MPI_Allreduce(&local,&global,1,LIS_MPI_INT,MPI_MAX,A->comm);
    return rc==MPI_SUCCESS ? global : LIS_ERR_ILL_ARG;
}

static LIS_INT lis_matvec_user_rows_mpi(
    LIS_MATRIX A,LIS_VECTOR X,LIS_VECTOR Y,LIS_INT adjoint)
{
    LIS_INT err=LIS_SUCCESS,n=0,nnz=0,i,j,k,required,got;
    LIS_INT *rowptr=NULL,*cols=NULL,*request_ids=NULL,*incoming_ids=NULL;
    LIS_INT *remote_slot=NULL;
    LIS_SCALAR *coeff=NULL,*local_x=NULL,*remote_vals=NULL,*incoming_vals=NULL;
    int *send_count=NULL,*recv_count=NULL,*send_at=NULL,*recv_at=NULL,*next=NULL;
    int np=0,rank=0,r,total_send=0,total_recv=0,rc;
    const size_t int_size=sizeof(LIS_INT), scalar_size=sizeof(LIS_SCALAR);

    if(A==NULL || X==NULL || Y==NULL)
        return LIS_ERR_ILL_ARG;
    if(A->nprocs<1 || A->nprocs>INT_MAX ||
       A->my_rank<0 || A->my_rank>=A->nprocs ||
       A->n<0 || A->n>INT_MAX || A->gn<0 || A->ranges==NULL ||
       A->is<0 || A->ie<A->is || A->ie>A->gn ||
       A->ie-A->is!=A->n || X->n<A->n || Y->n<A->n ||
       (A->n>0 && (X->value==NULL || Y->value==NULL)) ||
       A->user_get_row==NULL ||
       (uintmax_t)A->n+1>(uintmax_t)SIZE_MAX/int_size ||
       (uintmax_t)A->n>(uintmax_t)SIZE_MAX/scalar_size ||
       (uintmax_t)A->nprocs>(uintmax_t)SIZE_MAX/sizeof(int))
        err=LIS_ERR_ILL_ARG;

    if(!err)
    {
        np=(int)A->nprocs;
        rank=(int)A->my_rank;
        if(A->ranges[0]!=0 || A->ranges[np]!=A->gn ||
           A->ranges[rank]!=A->is || A->ranges[rank+1]!=A->ie)
            err=LIS_ERR_ILL_ARG;
        for(r=0;r<np && !err;r++)
            if(A->ranges[r]<0 || A->ranges[r+1]<A->ranges[r] ||
               A->ranges[r+1]>A->gn)
                err=LIS_ERR_ILL_ARG;
    }
    if(!err)
    {
        n=A->n;
        rowptr=(LIS_INT *)malloc(((size_t)n+1)*int_size);
        local_x=(LIS_SCALAR *)malloc(n ? (size_t)n*scalar_size : scalar_size);
        send_count=(int *)calloc((size_t)np,sizeof(int));
        recv_count=(int *)calloc((size_t)np,sizeof(int));
        send_at=(int *)calloc((size_t)np,sizeof(int));
        recv_at=(int *)calloc((size_t)np,sizeof(int));
        next=(int *)calloc((size_t)np,sizeof(int));
        if(!rowptr || !local_x || !send_count || !recv_count ||
           !send_at || !recv_at || !next)
            err=LIS_ERR_OUT_OF_MEMORY;
        else
        {
            if(n)memcpy(local_x,X->value,(size_t)n*scalar_size);
            rowptr[0]=0;
        }
    }
    /* All ranks stop together before entering application callbacks. */
    err=lis_user_sparse_sync_error(A,err);
    if(err)goto cleanup;

    for(i=0;i<n;i++)
    {
        required=0;
        err=A->user_get_row(A->user_data,i,&required,NULL,NULL);
        if(err)break;
        if(required<0 || required>INT_MAX-nnz)
        {
            err=LIS_ERR_ILL_ARG;
            break;
        }
        nnz+=required;
        rowptr[i+1]=nnz;
    }
    err=lis_user_sparse_sync_error(A,err);
    if(err)goto cleanup;

    if((uintmax_t)nnz>(uintmax_t)SIZE_MAX/int_size ||
       (uintmax_t)nnz>(uintmax_t)SIZE_MAX/scalar_size)
        err=LIS_ERR_ILL_ARG;
    if(!err)
    {
        cols=(LIS_INT *)malloc(nnz ? (size_t)nnz*int_size : int_size);
        coeff=(LIS_SCALAR *)malloc(
            nnz ? (size_t)nnz*scalar_size : scalar_size);
        remote_slot=(LIS_INT *)malloc(
            nnz ? (size_t)nnz*int_size : int_size);
        if(!cols || !coeff || !remote_slot)
            err=LIS_ERR_OUT_OF_MEMORY;
    }
    err=lis_user_sparse_sync_error(A,err);
    if(err)goto cleanup;

    for(i=0;i<n;i++)
    {
        required=rowptr[i+1]-rowptr[i];
        if(!required)continue;
        got=required;
        err=A->user_get_row(A->user_data,i,&got,
                             cols+rowptr[i],coeff+rowptr[i]);
        if(err)break;
        if(got!=required){err=LIS_ERR_ILL_ARG;break;}
        for(j=rowptr[i];j<rowptr[i+1];j++)
        {
            int owner;
            if(cols[j]<0 || cols[j]>=A->gn)
            {
                err=LIS_ERR_ILL_ARG;
                break;
            }
            owner=lis_user_sparse_owner(A->ranges,np,cols[j]);
            if(owner<0){err=LIS_ERR_ILL_ARG;break;}
            if(owner==rank)
                remote_slot[j]=-1;
            else
            {
                if(send_count[owner]==INT_MAX)
                {
                    err=LIS_ERR_ILL_ARG;
                    break;
                }
                send_count[owner]++;
            }
        }
        if(err)break;
    }
    err=lis_user_sparse_sync_error(A,err);
    if(err)goto cleanup;

    rc=MPI_Alltoall(send_count,1,MPI_INT,
                    recv_count,1,MPI_INT,A->comm);
    if(rc!=MPI_SUCCESS){err=LIS_ERR_ILL_ARG;goto cleanup;}
    for(r=0;r<np;r++)
    {
        if(send_count[r]<0 || recv_count[r]<0 ||
           send_count[r]>INT_MAX-total_send ||
           recv_count[r]>INT_MAX-total_recv)
        {
            err=LIS_ERR_ILL_ARG;
            break;
        }
        send_at[r]=total_send;
        recv_at[r]=total_recv;
        total_send+=send_count[r];
        total_recv+=recv_count[r];
    }
    if(!err &&
       ((uintmax_t)total_send>(uintmax_t)SIZE_MAX/int_size ||
        (uintmax_t)total_send>(uintmax_t)SIZE_MAX/scalar_size ||
        (uintmax_t)total_recv>(uintmax_t)SIZE_MAX/int_size ||
        (uintmax_t)total_recv>(uintmax_t)SIZE_MAX/scalar_size))
        err=LIS_ERR_ILL_ARG;
    if(!err)
    {
        request_ids=(LIS_INT *)malloc(
            total_send ? (size_t)total_send*int_size : int_size);
        incoming_ids=(LIS_INT *)malloc(
            total_recv ? (size_t)total_recv*int_size : int_size);
        remote_vals=(LIS_SCALAR *)malloc(
            total_send ? (size_t)total_send*scalar_size : scalar_size);
        incoming_vals=(LIS_SCALAR *)malloc(
            total_recv ? (size_t)total_recv*scalar_size : scalar_size);
        if(!request_ids || !incoming_ids || !remote_vals || !incoming_vals)
            err=LIS_ERR_OUT_OF_MEMORY;
    }
    err=lis_user_sparse_sync_error(A,err);
    if(err)goto cleanup;

    memcpy(next,send_at,(size_t)np*sizeof(int));
    for(j=0;j<nnz;j++)
    {
        if(remote_slot[j]==-1)continue;
        r=lis_user_sparse_owner(A->ranges,np,cols[j]);
        k=next[r]++;
        remote_slot[j]=k;
        request_ids[k]=cols[j];
    }
    rc=MPI_Alltoallv(request_ids,send_count,send_at,LIS_MPI_INT,
                     incoming_ids,recv_count,recv_at,LIS_MPI_INT,A->comm);
    if(rc!=MPI_SUCCESS){err=LIS_ERR_ILL_ARG;goto cleanup;}
    for(k=0;k<total_recv;k++)
        if(incoming_ids[k]<A->is || incoming_ids[k]>=A->ie)
        {
            err=LIS_ERR_ILL_ARG;
            break;
        }
    err=lis_user_sparse_sync_error(A,err);
    if(err)goto cleanup;

    if(adjoint)
    {
        if(n)memset(Y->value,0,(size_t)n*scalar_size);
        for(i=0;i<n;i++)
            for(j=rowptr[i];j<rowptr[i+1];j++)
            {
                LIS_SCALAR v=conj(coeff[j])*local_x[i];
                if(remote_slot[j]<0)
                    Y->value[cols[j]-A->is]+=v;
                else
                    remote_vals[remote_slot[j]]=v;
            }
        rc=MPI_Alltoallv(remote_vals,send_count,send_at,LIS_MPI_SCALAR,
                         incoming_vals,recv_count,recv_at,LIS_MPI_SCALAR,
                         A->comm);
        if(rc!=MPI_SUCCESS){err=LIS_ERR_ILL_ARG;goto cleanup;}
        for(k=0;k<total_recv;k++)
            Y->value[incoming_ids[k]-A->is]+=incoming_vals[k];
    }
    else
    {
        for(k=0;k<total_recv;k++)
            incoming_vals[k]=local_x[incoming_ids[k]-A->is];
        rc=MPI_Alltoallv(incoming_vals,recv_count,recv_at,LIS_MPI_SCALAR,
                         remote_vals,send_count,send_at,LIS_MPI_SCALAR,
                         A->comm);
        if(rc!=MPI_SUCCESS){err=LIS_ERR_ILL_ARG;goto cleanup;}
        for(i=0;i<n;i++)
        {
            LIS_SCALAR sum=(LIS_SCALAR)0.0;
            for(j=rowptr[i];j<rowptr[i+1];j++)
            {
                LIS_SCALAR v=remote_slot[j]<0 ?
                             local_x[cols[j]-A->is] :
                             remote_vals[remote_slot[j]];
                sum+=coeff[j]*v;
            }
            Y->value[i]=sum;
        }
    }
cleanup:
    free(rowptr);free(cols);free(coeff);free(remote_slot);free(local_x);
    free(send_count);free(recv_count);free(send_at);free(recv_at);
    free(next);free(request_ids);free(incoming_ids);
    free(remote_vals);free(incoming_vals);
    return err;
}

#endif

#undef __FUNC__
#define __FUNC__ "lis_matvec"
LIS_INT lis_matvec(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT err,i;
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	if( A->matrix_type==LIS_MATRIX_OPERATOR )
	{
		if( X->precision!=LIS_PRECISION_DEFAULT ||
		    Y->precision!=LIS_PRECISION_DEFAULT )
		{
			LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
			           "LIS_MATRIX_OPERATOR does not support quad precision yet\n");
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		if( A->operator_scale!=NULL )
		{
		        if( A->operator_A==NULL ||
		            A->operator_work==NULL )
		        {
		                LIS_SETERR(
		                        LIS_ERR_ILL_ARG,
		                        "scaled matrix operator is not initialized\n");
		                return LIS_ERR_ILL_ARG;
		        }

		        /* work = D*x */
		        for(i=0;i<A->n;i++)
		        {
		                A->operator_work->value[i] =
		                        A->operator_scale->value[i]
		                        * X->value[i];
		        }

		        /* y = A*work */
		        err = lis_matvec(
		                A->operator_A,
		                A->operator_work,
		                Y);
		        if( err ) return err;

		        /* y = D*y */
		        for(i=0;i<A->n;i++)
		        {
		                Y->value[i] *= A->operator_scale->value[i];
		        }

		        LIS_DEBUG_FUNC_OUT;
		        return LIS_SUCCESS;
		}

		if( A->operator_A==NULL || A->operator_B==NULL ||
		    A->operator_work==NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,
			           "LIS_MATRIX_OPERATOR is not initialized\n");
			return LIS_ERR_ILL_ARG;
		}

		err = lis_matvec(A->operator_A,X,A->operator_work);
		if( err ) return err;
		err = lis_matvec(A->operator_B,X,Y);
		if( err ) return err;
		err = lis_vector_scale(A->operator_alpha,A->operator_work);
		if( err ) return err;
		err = lis_vector_scale(A->operator_beta,Y);
		if( err ) return err;
		err = lis_vector_axpy((LIS_SCALAR)1.0,A->operator_work,Y);
		if( err ) return err;

		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}

	if( X->precision==LIS_PRECISION_DEFAULT )
	{
		x = X->value;
		y = Y->value;

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_csr(A, x, y);
			break;
		case LIS_MATRIX_BSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			if( A->bnr<=4 && A->bnc<=4 )
			{
				lis_matvec_bsr_xxx[A->bnr-1][A->bnc-1](A, x, y);
			}
			else
			{
				lis_matvec_bsr(A, x, y);
			}
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_csc(A, x, y);
			break;
		case LIS_MATRIX_BSC:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_bsc(A, x, y);
			break;
		case LIS_MATRIX_MSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_msr(A, x, y);
			break;
		case LIS_MATRIX_ELL:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_ell(A, x, y);
			break;
		case LIS_MATRIX_DIA:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_dia(A, x, y);
			break;
		case LIS_MATRIX_JAD:
			#ifdef USE_MPI
				#ifndef USE_OVERLAP
					LIS_MATVEC_SENDRECV;
				#else
					LIS_MATVEC_REALLOC;
				#endif
			#endif
			lis_matvec_jad(A, x, y);
			break;
		case LIS_MATRIX_VBR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_vbr(A, x, y);
			break;
		case LIS_MATRIX_DNS:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_dns(A, x, y);
			break;
		case LIS_MATRIX_COO:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_coo(A, x, y);
			break;
		case LIS_MATRIX_USER:
			if( A->user_matvec!=NULL )
			{
				err = A->user_matvec(A->user_data,x,y);
			}
			else if( A->user_get_row!=NULL )
			{
#ifdef USE_MPI
				if(A->nprocs>1)
					err=lis_matvec_user_rows_mpi(A,X,Y,LIS_FALSE);
				else
#endif
					err=lis_matvec_user_rows(A,x,y);
			}
			else
			{
				LIS_SETERR(
					LIS_ERR_NOT_IMPLEMENTED,
					"user matvec callback and row callback are not set\n");
				return LIS_ERR_NOT_IMPLEMENTED;
			}
			if( err ) return err;
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#ifdef USE_QUAD_PRECISION
	else
	{

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				lis_send_recv_mp(A->commtable,X);
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvec_csr_mp(A, X, Y);
			#else
				lis_matvec_csr_mp2(A, X, Y);
			#endif
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				lis_send_recv_mp(A->commtable,X);
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvec_csc_mp(A, X, Y);
			#else
				lis_matvec_csc_mp2(A, X, Y);
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matvech"
LIS_INT lis_matvech(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_INT err,i;
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	if( A->matrix_type==LIS_MATRIX_OPERATOR )
	{
		if( X->precision!=LIS_PRECISION_DEFAULT ||
		    Y->precision!=LIS_PRECISION_DEFAULT )
		{
			LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
			           "LIS_MATRIX_OPERATOR does not support quad precision yet\n");
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		if( A->operator_scale!=NULL )
		{
		        if( A->operator_A==NULL ||
		            A->operator_work==NULL )
		        {
		                LIS_SETERR(
		                        LIS_ERR_ILL_ARG,
		                        "scaled matrix operator is not initialized\n");
		                return LIS_ERR_ILL_ARG;
		        }

		        /* work = D*x */
		        for(i=0;i<A->n;i++)
		        {
		                A->operator_work->value[i] =
#ifdef _COMPLEX
		                        conj(A->operator_scale->value[i])
#else
		                        A->operator_scale->value[i]
#endif
		                        * X->value[i];
		        }

		        /* y = A*work */
		        err = lis_matvech(
		                A->operator_A,
		                A->operator_work,
		                Y);
		        if( err ) return err;

		        /* y = D*y */
		        for(i=0;i<A->n;i++)
		        {
#ifdef _COMPLEX
		                Y->value[i] *= conj(A->operator_scale->value[i]);
#else
		                Y->value[i] *= A->operator_scale->value[i];
#endif
		        }

		        LIS_DEBUG_FUNC_OUT;
		        return LIS_SUCCESS;
		}

		if( A->operator_A==NULL || A->operator_B==NULL ||
		    A->operator_work==NULL )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,
			           "LIS_MATRIX_OPERATOR is not initialized\n");
			return LIS_ERR_ILL_ARG;
		}

		err = lis_matvech(A->operator_A,X,A->operator_work);
		if( err ) return err;
		err = lis_matvech(A->operator_B,X,Y);
		if( err ) return err;
		err = lis_vector_scale(conj(A->operator_alpha),A->operator_work);
		if( err ) return err;
		err = lis_vector_scale(conj(A->operator_beta),Y);
		if( err ) return err;
		err = lis_vector_axpy((LIS_SCALAR)1.0,A->operator_work,Y);
		if( err ) return err;

		LIS_DEBUG_FUNC_OUT;
		return LIS_SUCCESS;
	}


	x = X->value;
	y = Y->value;

	if( X->precision==LIS_PRECISION_DEFAULT )
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_csr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_BSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_bsr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_csc(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_BSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_bsc(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_MSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_msr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_ELL:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_ell(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_JAD:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_jad(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_DIA:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_dia(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_VBR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_vbr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_DNS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_dns(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_COO:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvech_coo(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_USER:
			if( A->user_matvech!=NULL )
			{
				err = A->user_matvech(A->user_data,x,y);
			}
			else if( A->user_get_row!=NULL )
			{
#ifdef USE_MPI
				if(A->nprocs>1)
					err=lis_matvec_user_rows_mpi(A,X,Y,LIS_TRUE);
				else
#endif
					err=lis_matvech_user_rows(A,x,y);
			}
			else
			{
				LIS_SETERR(
					LIS_ERR_NOT_IMPLEMENTED,
					"user transpose/Hermitian matvec callback and row callback are not set\n");
				return LIS_ERR_NOT_IMPLEMENTED;
			}
			if( err ) return err;
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#ifdef USE_QUAD_PRECISION
	else
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvech_csr_mp(A, X, Y);
			#else
				lis_matvech_csr_mp2(A, X, Y);
			#endif
			#ifdef USE_MPI
				lis_reduce_mp(A->commtable,Y);
			#endif
			break;
		case LIS_MATRIX_CSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvech_csc_mp(A, X, Y);
			#else
				lis_matvech_csc_mp2(A, X, Y);
			#endif
			#ifdef USE_MPI
				lis_reduce_mp(A->commtable,Y);
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matvec_block"
LIS_INT lis_matvec_block(
        LIS_MATRIX A,
        LIS_INT nvec,
        LIS_VECTOR X[],
        LIS_VECTOR Y[])
{
        LIS_INT i,err;
        LIS_INT use_user_block;

        LIS_DEBUG_FUNC_IN;

        if( nvec<0 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "number of block vectors is negative\n");
                return LIS_ERR_ILL_ARG;
        }

        if( nvec==0 )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( A==NULL || X==NULL || Y==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "invalid block matvec argument\n");
                return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<nvec;i++)
        {
                if( X[i]==NULL || Y[i]==NULL )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "block matvec contains a null vector\n");
                        return LIS_ERR_ILL_ARG;
                }
        }

        /*
         * A USER matrix may provide a native block callback.
         *
         * The USER callback interface is currently default-precision.
         * If a vector requests another precision, preserve the existing
         * scalar lis_matvec() behavior through the fallback below.
         */
        use_user_block =
                A->matrix_type==LIS_MATRIX_USER &&
                A->user_matvec_block!=NULL;

        if( use_user_block )
        {
                for(i=0;i<nvec;i++)
                {
                        if( X[i]->precision!=LIS_PRECISION_DEFAULT )
                        {
                                use_user_block = LIS_FALSE;
                                break;
                        }
                }
        }

        if( use_user_block )
        {
                err = A->user_matvec_block(
                        A->user_data,
                        nvec,
                        X,
                        Y);

                if( err ) return err;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        /*
         * Compatibility fallback.
         *
         * This deliberately reuses lis_matvec() so all existing matrix
         * formats, operator semantics, MPI communication, and precision
         * handling remain unchanged.
         */
        for(i=0;i<nvec;i++)
        {
                err = lis_matvec(
                        A,
                        X[i],
                        Y[i]);

                if( err ) return err;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matvech_block"
LIS_INT lis_matvech_block(
        LIS_MATRIX A,
        LIS_INT nvec,
        LIS_VECTOR X[],
        LIS_VECTOR Y[])
{
        LIS_INT i,err;
        LIS_INT use_user_block;

        LIS_DEBUG_FUNC_IN;

        if( nvec<0 )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "number of block vectors is negative\n");
                return LIS_ERR_ILL_ARG;
        }

        if( nvec==0 )
        {
                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        if( A==NULL || X==NULL || Y==NULL )
        {
                LIS_SETERR(
                        LIS_ERR_ILL_ARG,
                        "invalid block matvech argument\n");
                return LIS_ERR_ILL_ARG;
        }

        for(i=0;i<nvec;i++)
        {
                if( X[i]==NULL || Y[i]==NULL )
                {
                        LIS_SETERR(
                                LIS_ERR_ILL_ARG,
                                "block matvech contains a null vector\n");
                        return LIS_ERR_ILL_ARG;
                }
        }

        use_user_block =
                A->matrix_type==LIS_MATRIX_USER &&
                A->user_matvech_block!=NULL;

        if( use_user_block )
        {
                for(i=0;i<nvec;i++)
                {
                        if( X[i]->precision!=LIS_PRECISION_DEFAULT )
                        {
                                use_user_block = LIS_FALSE;
                                break;
                        }
                }
        }

        if( use_user_block )
        {
                err = A->user_matvech_block(
                        A->user_data,
                        nvec,
                        X,
                        Y);

                if( err ) return err;

                LIS_DEBUG_FUNC_OUT;
                return LIS_SUCCESS;
        }

        for(i=0;i<nvec;i++)
        {
                err = lis_matvech(
                        A,
                        X[i],
                        Y[i]);

                if( err ) return err;
        }

        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matvec_optimize"
LIS_INT lis_matvec_optimize(LIS_MATRIX A, LIS_INT *matrix_type_maxperf)
{
  int    	   nprocs,my_rank;
  LIS_INT	   err,iter,i,matrix_type,ss,se;
  double	   time,time2,convtime;
  LIS_REAL         val;
  double	   commtime,comptime,flops,flops_maxperf;
  LIS_MATRIX       A1;
  LIS_VECTOR       X, Y;
  char             *lis_storagename2[]   = {"CSR", "CSC", "MSR", "DIA", "ELL", "JAD", "BSR", "BSC", "VBR", "COO", "DNS"};

	LIS_DEBUG_FUNC_IN;

#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
	nprocs  = 1;
	my_rank = 0;
#endif

	ss = 1;
	se = 11;

	lis_vector_duplicate(A,&X);
	lis_vector_duplicate(A,&Y);
	err = lis_vector_set_all(1.0,X);

	printf("\nmeasuring matvec performance...\n");
	iter = (int)(10000000 / A->nnz) + 1; 
	flops_maxperf = 0.0;
#ifdef _LONG__LONG
	printf("number of iterations = 1e7 / %lld + 1 = %lld\n", A->nnz, iter);
#else
	printf("number of iterations = 1e7 / %d + 1 = %d\n", A->nnz, iter);
#endif
	for (matrix_type=ss;matrix_type<se;matrix_type++)
	  {
	    if ( nprocs>1 && matrix_type==9 ) continue;
	    lis_matrix_duplicate(A,&A1);
	    lis_matrix_set_type(A1,matrix_type);
	    err = lis_matrix_convert(A,A1);
	    if( err ) CHKERR(err);
		    
	    comptime = 0.0;
	    commtime = 0.0;

	    for(i=0;i<iter;i++)
	      {
#ifdef USE_MPI
		MPI_Barrier(A1->comm);
		time = lis_wtime();
		lis_send_recv(A1->commtable,X->value);
		commtime += lis_wtime() - time;
#endif
		time2 = lis_wtime();
		lis_matvec(A1,X,Y);
		comptime += lis_wtime() - time2;
	      }
	    lis_vector_nrm2(Y,&val);

	    if( my_rank==0 )
	      {
		flops = 2.0*A->nnz*iter*1.0e-6 / comptime;
#ifdef USE_MPI
#ifdef _LONG__DOUBLE
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#endif
#else
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100);
#endif
#endif
#else
#ifdef _LONG__DOUBLE
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#endif
#else
#ifdef _LONG__LONG
		printf("matrix_type = %2lld (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#else
		printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops);
#endif
#endif
#endif
	      }
	    lis_matrix_destroy(A1);

	    if ( flops > flops_maxperf )
	      {
		*matrix_type_maxperf=matrix_type;
		flops_maxperf = flops;
	      }
	  }

	printf("matrix format is set to %s\n\n", lis_storagename2[*matrix_type_maxperf-1]);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
