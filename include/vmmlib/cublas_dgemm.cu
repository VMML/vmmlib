#ifndef __VMML__VMMLIB_CUBLAS_DGEMM__HPP__
#define __VMML__VMMLIB_CUBLAS_DGEMM__HPP__


#include <vmmlib/matrix.hpp>
#include <vmmlib/exception.hpp>
#include <vmmlib/cublas_includes.hpp>
#include <vmmlib/cublas_types.hpp>

/** 
 *
 *   a wrapper for CUBLAS DGEMM routine. 
 
 cublasStatus_t cublasDgemm(
			cublasHandle_t handle,
			cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k,
			const double const double const double const double double
			*alpha, *A, int lda, *B, int ldb, *beta,
			*C, int ldc
			)
 
 *
 *  Purpose
 *  =======
 *
 *  CUBLAS DGEMM  is a CUDA implementation of the level3 BLAS DGEMM implementation
 *  performs one of the matrix-matrix operations
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
	
	namespace cublas
	{
		
#if 0
		/* CUBLAS DGEMM Subroutine */ 
		cublasStatus_t cublasDgemm(
								   cublasHandle_t handle,
								   cublasOperation_t transa, cublasOperation_t transb, 
								   int m, int n, int k,
								   const double *alpha, 
								   const double *A, int lda,
								   const double *B, int ldb,
								   const *beta,
								   double *C, int ldc
								   )

		/* DGEMM Subroutine */ 
		void cblas_dgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, 
						 blasint M, blasint N, blasint K,
						 double alpha, double *A, blasint lda, double *B, blasint ldb, double beta, double *C, blasint ldc);
		
#endif
		
		template< typename float_t >
		struct dgemm_params
		{
			//FIXME maybe add const 
			cublasHandle_t		handle;
			cublasOperation_t	trans_a;
			cublasOperation_t	trans_b;
			cublas_int			m;
			cublas_int			n;
			cublas_int			k;
			float_t				alpha;
			float_t*			h_a; //host
			float_t*			d_a; //device
			cublas_int			lda; //leading dimension of input array matrix left
			float_t*			h_b;
			float_t*			d_b;
			cublas_int			ldb; //leading dimension of input array matrix right
			float_t				beta;
			float_t*			h_c;
			float_t*			d_c;
			cublas_int			ldc; //leading dimension of output array matrix right
			
			friend std::ostream& operator << ( std::ostream& os, 
											  const dgemm_params< float_t >& p )
			{
				os 
				<< " (1)\thandle "		<< p.handle << std::endl
				<< " (2)\ttrans_a "		<< p.trans_a << std::endl
				<< " (3)\ttrans_b "     << p.trans_b << std::endl
				<< " (4)\tm "			<< p.m << std::endl
				<< " (6)\tn "			<< p.n << std::endl
				<< " (5)\tk "			<< p.k << std::endl
				<< " (7)\talpha "       << p.alpha << std::endl 
				<< " (8)\th_a "			<< p.h_a << std::endl
				<< " (9)\tlda "			<< p.lda << std::endl
				<< " (10)\th_b "		<< p.h_b << std::endl
				<< " (11)\tldb "		<< p.ldb << std::endl
				<< " (12)\tbeta "       << p.beta << std::endl
				<< " (13)\th_c "		<< p.h_c << std::endl
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
			//std::cout << "calling cublas sgemm (single precision) " << std::endl;
			cublasStatus_t stat  = cublasSgemm( 
						p.handle,
						p.trans_a,
						p.trans_b,
						p.m,
						p.n,
						p.k,
						&p.alpha,
						p.d_a,
						p.lda,
						p.d_b,
						p.ldb,
						&p.beta,
						p.d_c,
						p.ldc
						);
			
		}
		
		template<>
		inline void
		dgemm_call( dgemm_params< double >& p )
		{
			//std::cout << "calling cublas dgemm (double precision) " << std::endl;

			cublasStatus_t stat  = cublasDgemm( 
						p.handle,
						p.trans_a,
						p.trans_b,
						p.m,
						p.n,
						p.k,
						&p.alpha,
						p.d_a,
						p.lda,
						p.d_b,
						p.ldb,
						&p.beta,
						p.d_c,
						p.ldc
						);
		}
		
	} // namespace cublas
	
	
	
	template< size_t M, size_t K, size_t N, typename float_t >
	struct cublas_dgemm
	{
		
		typedef matrix< M, K, float_t > matrix_left_t;
		typedef matrix< K, M, float_t > matrix_left_t_t;
		typedef matrix< K, N, float_t > matrix_right_t;
		typedef matrix< N, K, float_t > matrix_right_t_t;
		typedef matrix< M, N, float_t > matrix_out_t;
		typedef vector< M, float_t > vector_left_t;
		typedef vector< N, float_t > vector_right_t;
		
		cublas_dgemm();
		~cublas_dgemm();
		
		bool compute( const matrix_left_t& A_, const matrix_right_t& B_, matrix_out_t& C_ );
		bool compute( const matrix_left_t& A_, matrix_out_t& C_ );
		
		
		cublas::dgemm_params< float_t > p;
		
		const cublas::dgemm_params< float_t >& get_params(){ return p; };
		
		
	}; // struct cublas_dgemm
	
	
	template< size_t M, size_t K, size_t N, typename float_t >
	cublas_dgemm< M, K, N, float_t >::cublas_dgemm()
	{
        cublasStatus_t cstat = cublasCreate( &p.handle ); if ( cstat > 0 ) { printf( "cublasCreate: status error=%d\n", cstat ); }
		p.trans_a    = CUBLAS_OP_N;
		p.trans_b    = CUBLAS_OP_N;
		p.m          = M;
		p.n          = N;
		p.k          = K;
		p.alpha      = 1.0f;
		p.h_a        = 0;
		p.d_a        = 0;
		p.lda        = M;
		p.h_b        = 0;
		p.d_b        = 0;
		p.ldb        = K; //no transpose, use N for transpose
		p.beta       = 0.0;
		p.h_c        = 0;
		p.d_c        = 0;
		p.ldc        = M;
	}
	
	template< size_t M, size_t K, size_t N, typename float_t >
	cublas_dgemm< M, K, N, float_t >::~cublas_dgemm()
	{
		/*cublasStatus_t cuerr = cublasDestroy( p.handle );
		if ( cuerr > 0 ) 
		{ printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }*/
	}

	
	
	
	template< size_t M, size_t K, size_t N, typename float_t >
	bool
	cublas_dgemm< M, K, N, float_t >::compute( 
											const matrix_left_t& A_, 
											const matrix_right_t& B_,
											matrix_out_t& C_ 
											)
	{
		// cublas needs non-const data
		matrix_left_t* AA = new matrix_left_t( A_ );
		matrix_right_t* BB = new matrix_right_t( B_ );
		C_.zero();
		
		p.h_a         = AA->array;
		p.h_b         = BB->array;
		p.h_c         = C_.array;
				
		// memory sizes of matrices
		size_t mem_size_A = sizeof(float_t) * M * K;
		size_t mem_size_B = sizeof(float_t) * K * N;
		size_t mem_size_C = sizeof(float_t) * M * N;
		
		// allocate device memory
		cudaError_t cuerr = cudaMalloc( (void**) &p.d_a, mem_size_A ); if ( cuerr > 0 ) { printf( "cudaMalloc: cublas error=%d\n", cuerr ); }
		cuerr = cudaMalloc( (void**) &p.d_b, mem_size_B ); if ( cuerr > 0 ) { printf( "cudaMalloc: cublas error=%d\n", cuerr ); } 
		cuerr = cudaMalloc( (void**) &p.d_c, mem_size_C ); if ( cuerr > 0 ) { printf( "cudaMalloc: cublas error=%d\n", cuerr ); }  
		
		// copy host memory to device
		cuerr = cudaMemcpy( p.d_a, p.h_a, mem_size_A, cudaMemcpyHostToDevice); if ( cuerr > 0 ) { printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }
		cuerr = cudaMemcpy( p.d_b, p.h_b, mem_size_B, cudaMemcpyHostToDevice); if ( cuerr > 0 ) { printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }

		// call CUBLAS V2
		cublas::dgemm_call< float_t >( p );
		//std::cout << p << std::endl; //debug
	
		cudaDeviceSynchronize();

		// copy result from device to host
		cuerr = cudaMemcpy( p.h_c, p.d_c, mem_size_C, cudaMemcpyDeviceToHost); if ( cuerr > 0 ) { printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }
		
		// clean up memory
		cudaFree( p.d_a );
		cudaFree( p.d_b );
		cudaFree( p.d_c );
		cudaDeviceReset();
		
		delete AA;
		delete BB;
		
		return true;
	}	
	
	
	template< size_t M, size_t K, size_t N, typename float_t >
	bool
	cublas_dgemm< M, K, N, float_t >::compute( 
											const matrix_left_t& A_, 
											matrix_out_t& C_ 
											)
	{
		// cublas needs non-const data
		matrix_left_t* AA = new matrix_left_t( A_ );
		C_.zero();
		
		p.h_a         = AA->array;
		p.h_b         = AA->array;
		p.h_c         = C_.array;
		p.trans_b     = CUBLAS_OP_T;
		p.ldb         = N; 
		
		// memory sizes of matrices
		size_t mem_size_A = sizeof(float_t) * M * K;
		size_t mem_size_B = sizeof(float_t) * K * N;
		size_t mem_size_C = sizeof(float_t) * M * N;
		
		// allocate device memory
		cudaError_t cuerr = cudaMalloc( (void**) &p.d_a, mem_size_A ); if ( cuerr > 0 ) { printf( "cudaMalloc: cublas error=%d\n", cuerr ); }
		cuerr = cudaMalloc( (void**) &p.d_b, mem_size_B ); if ( cuerr > 0 ) { printf( "cudaMalloc: cublas error=%d\n", cuerr ); } 
		cuerr = cudaMalloc( (void**) &p.d_c, mem_size_C ); if ( cuerr > 0 ) { printf( "cudaMalloc: cublas error=%d\n", cuerr ); }  
		
		// copy host memory to device
		cuerr = cudaMemcpy( p.d_a, p.h_a, mem_size_A, cudaMemcpyHostToDevice); if ( cuerr > 0 ) { printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }
		cuerr = cudaMemcpy( p.d_b, p.h_b, mem_size_B, cudaMemcpyHostToDevice); if ( cuerr > 0 ) { printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }

		// call CUBLAS V2
		cublas::dgemm_call< float_t >( p );
		//std::cout << p << std::endl; //debug
	
		cudaDeviceSynchronize();

		// copy result from device to host
		cuerr = cudaMemcpy( p.h_c, p.d_c, mem_size_C, cudaMemcpyDeviceToHost); if ( cuerr > 0 ) { printf( "cudaMemcpy: cublas error=%d\n", cuerr ); }
		
		// clean up memory
		cudaFree( p.d_a );
		cudaFree( p.d_b );
		cudaFree( p.d_c );
		cudaDeviceReset();

		delete AA;
		
		return true;
	}	

	
	
} // namespace vmml

#endif	

