#include "cublas_dgemm_test.hpp"

#ifdef VMMLIB_USE_CUDA

#include <vmmlib/cublas_dgemm.cu>

#define TEST( x ) \
{ \
    ok = x; \
    global_ok &= ok; \
}

namespace vmml
{
	
	bool
	cublas_dgemm_test::run()
	{
        bool global_ok = true;
		bool ok = false;

		//cublas compute
		//A*B = D (MxK, KxN, MxN)
		matrix< 3, 6, float > A;
		matrix< 6, 2, float > B;
		matrix< 3, 2, float > D;
		matrix< 3, 2, float > D_check;
		
		float AData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
		A = AData;
		float BData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		B = BData;
		
		cublas_dgemm< 3, 6, 2, float > cublas_multiplier;
		cublas_multiplier.compute( A, B, D );
		
		float DData[] = { 161, 182, 377, 434, 593, 686 };
		D_check = DData;
		
		TEST(D == D_check);
		
		log( "matrix-matrix multiplication (MxK) x (KxN) = (MxN)", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix A\n" << A << "\n"
            << "input matrix B\n" << B << "\n"
            << "matrix C should be\n" << D_check << "\n"
            << "matrix C is\n" << D << "\n"
            << std::endl;
			log_error( ss.str() );
		}
		
		//covariance computation
		
		matrix< 3, 3, float > C;
		matrix< 3, 3, float > C_check;
		
		cublas_dgemm< 3, 6, 3, float > blas_cov;
		blas_cov.compute( A, C );
		
		float CData[] = { 91, 217, 343, 217, 559, 901, 343, 901, 1459 };
		C_check = CData;
		
		TEST(C == C_check);
		
		log( "symmetric matrix-matrix multiplication (input left matrix) (MxK) x (KxN) = (MxN), while M=N", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix (left matrix)\n" << A << "\n"
            << "covariance matrix should be\n" << C_check << "\n"
            << "covariance matrix is\n" << C << "\n"
            << std::endl;
			log_error( ss.str() );
		}
		
		
		return global_ok;
	}
	
	
	
} // namespace vmml
#endif
