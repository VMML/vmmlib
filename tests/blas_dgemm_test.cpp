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
		
		log( "symmetric matrix matrix multiplication", ok );
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
		
		return ok;
		return true;
	}
	
	
	
} // namespace vmml

