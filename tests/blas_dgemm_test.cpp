#include "blas_dgemm_test.hpp"

#include <vmmlib/blas_dgemm.hpp>

namespace vmml
{
	
	bool
	blas_dgemm_test::run()
	{
		bool ok = false;
		
		matrix< 3, 6, double > A;
		matrix< 3, 3, double > C;
		matrix< 3, 3, double > C_check;
		
		double AData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
		A = AData;
		
		blas_dgemm< 3, 6, 3, double > blas_cov;
		blas_cov.compute( A, C );
		
		double CData[] = { 91, 217, 343, 217, 559, 901, 343, 901, 1459 };
		C_check = CData;
		
		ok = C == C_check;
		
		log( "symmetric matrix matrix multiplication (MxK) x (KxN) = (MxN), while M=N", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix\n" << A << "\n"
            << "covariance matrix should be\n" << C_check << "\n"
            << "covariance matrixis\n" << C << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		
		//A*B = D (MxK, KxN, MxN)
		matrix< 6, 2, double > B;
		matrix< 3, 2, double > D;
		matrix< 3, 2, double > D_check;
		
		double BData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		B = BData;

		blas_dgemm< 3, 6, 2, double > blas_v_dgemm;
		blas_v_dgemm.compute( A, B, D );

		double DData[] = { 161, 182, 377, 434, 593, 686 };
		D_check = DData;
		
		ok = D == D_check;
		
		log( "matrix matrix multiplication (MxK) x (KxN) = (MxN)", ok );
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
		
		
		
		return ok;
		return true;
	}
	
	
	
} // namespace vmml

