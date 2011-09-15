#ifndef __VMML__VMMLIB_BLAS_DGEMM__HPP__
#define __VMML__VMMLIB_BLAS_DGEMM__HPP__


#include <vmmlib/matrix.hpp>
#include <vmmlib/exception.hpp>
#include <vmmlib/blas_includes.hpp>
#include <vmmlib/blas_types.hpp>

/** 
 *
 *   a wrapper for blas's DGEMM routine. 
 
 SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
 *     .. Scalar Arguments ..
 DOUBLE PRECISION ALPHA,BETA
 INTEGER K,LDA,LDB,LDC,M,N
 CHARACTER TRANSA,TRANSB
 *     ..
 *     .. Array Arguments ..
 DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGEMM  performs one of the matrix-matrix operations
 *
 *     C := alpha*op( A )*op( B ) + beta*C,
 *
 *  where  op( X ) is one of
 *
 *     op( X ) = X   or   op( X ) = X**T,
 *
 *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
 *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 *
 *
 *   more information in: http://www.netlib.org/blas/dgemm.f 
 *   or http://www.netlib.org/clapack/cblas/dgemm.c
 **
 */


namespace vmml
{
	
	namespace blas
	{
				
		
#if 0
		/* Subroutine */ 
		void cblas_dgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, 
						 blasint M, blasint N, blasint K,
						 double alpha, double *A, blasint lda, double *B, blasint ldb, double beta, double *C, blasint ldc);
		
#endif

		template< typename float_t >
		struct dgemm_params
		{
			CBLAS_ORDER     order;
			CBLAS_TRANSPOSE trans_a;
			CBLAS_TRANSPOSE trans_b;
			blas_int 		m;
			blas_int		n;
			blas_int		k;
			float_t			alpha;
			float_t*        a;
			blas_int        lda; //leading dimension of input array matrix left
			float_t*        b;
			blas_int        ldb; //leading dimension of input array matrix right
			float_t			beta;
			float_t*        c;
			blas_int        ldc; //leading dimension of output array matrix right
			
			friend std::ostream& operator << ( std::ostream& os, 
											  const dgemm_params< float_t >& p )
			{
				os 
				<< " (1)\torder "     << p.order << std::endl
				<< " (2)\ttrans_a "    << p.trans_a << std::endl
				<< " (3)\ttrans_b "     << p.trans_b << std::endl
				<< " (4)\tm "        << p.m << std::endl
				<< " (6)\tn "      << p.n << std::endl
				<< " (5)\tk "        << p.k << std::endl
				<< " (7)\talpha "       << p.alpha << std::endl 
				<< " (8)\ta "       << p.a << std::endl
				<< " (9)\tlda "       << p.lda << std::endl
				<< " (10)\tb "       << p.b << std::endl
				<< " (11)\tldb "   << p.ldb << std::endl
				<< " (12)\tbeta "        << p.beta << std::endl
				<< " (13)\tc "        << p.c << std::endl
				<< " (14)\tldc "        << p.ldc << std::endl
				<< std::endl;
				return os;
			}
			
		};
		
		
		
		template< typename float_t >
		inline void
		dgemm_call( dgemm_params< float_t >& p )
		{
			VMMLIB_ERROR( "not implemented for this type.", VMMLIB_HERE );
		}
		
		
		template<>
		inline void
		dgemm_call( dgemm_params< float >& p )
		{
			std::cout << "calling blas sgemm (single precision) " << std::endl;
			cblas_sgemm( 
					p.order,
					p.trans_a,
					p.trans_b,
					p.m,
					p.n,
					p.k,
					p.alpha,
					p.a,
					p.lda,
					p.b,
					p.ldb,
				    p.beta,
					p.c,
					p.ldc
					);
			
		}
		
		
		template<>
		inline void
		dgemm_call( dgemm_params< double >& p )
		{
			std::cout << "calling blas dgemm (double precision) " << std::endl;
			cblas_dgemm( 
				   p.order,
				   p.trans_a,
				   p.trans_b,
				   p.m,
				   p.n,
				   p.k,
				   p.alpha,
				   p.a,
				   p.lda,
				   p.b,
				   p.ldb,
				   p.beta,
				   p.c,
				   p.ldc
				   );
		}
		
	} // namespace blas
	
	
	
	template< size_t M, size_t K, size_t N, typename float_t >
	struct blas_dgemm
	{
		
		typedef matrix< M, K, float_t > matrix_left_t;
		typedef matrix< K, N, float_t > matrix_right_t;
		typedef matrix< M, N, float_t > matrix_out_t;
		//typedef typename evalues_type::const_iterator evalue_const_iterator;
		
		
		blas_dgemm();
		~blas_dgemm() {};
		
		//bool compute( const matrix_left_t& A, const matrix_right_t& B, matrix_out_t& C );
		bool compute( const matrix_left_t& ABT, matrix_out_t& C );
				
		//inline bool test_success( blas::lapack_int info );
		
		blas::dgemm_params< float_t > p;
		
		const blas::dgemm_params< float_t >& get_params(){ return p; };
		
		
	}; // struct blas_dgemm
	
	
	template< size_t M, size_t K, size_t N, typename float_t >
	blas_dgemm< M, K, N, float_t >::blas_dgemm()
	{
		p.order      = CblasColMajor; //
		p.trans_a    = CblasNoTrans;
		p.trans_b    = CblasNoTrans;
		p.m          = M;
		p.n          = N;
		p.k          = K;
		p.alpha      = 1;
		p.a          = 0;
		p.lda        = M;
		p.b          = 0;
		p.ldb        = K; //no transpose
		p.beta       = 1;
		p.c          = 0;
		p.ldc        = M;
	}
	
	
		
#if 0
	template< size_t M, size_t K, size_t N, typename float_t >
	bool
	blas_dgemm< M, K, N, float_t >::compute( const matrix_left_t& A, 
									  const matrix_right_t& B,
									  matrix_out_t& C )
	{
		// lapack destroys the contents of the input matrix
		matrix_left_t AA( A );
		matrix_right_t BB( B );
		
		p.a         = AA.array;
		p.b         = BB.array;
		p.c         = C.array;
		
		blas::dgemm_call< float_t >( p );
		
		return true;
				
	}	
#endif	
	template< size_t M, size_t K, size_t N, typename float_t >
	bool
	blas_dgemm< M, K, N, float_t >::compute( const matrix_left_t& ABT, matrix_out_t& C )
	{
		// lapack destroys the contents of the input matrix
		matrix_left_t* AA = new matrix_left_t( ABT );
		
		p.trans_b   = CblasTrans;
		p.a         = AA->array;
		p.b         = AA->array;
		p.ldb       = N; 
		p.c         = C.array;
		
		blas::dgemm_call< float_t >( p );
		
		//std::cout << p << std::endl; //debug

		delete AA;
		
		return true;
	}	
	
	
	
	
} // namespace vmml

#endif	

