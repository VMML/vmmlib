#include "cublas_dgemm_test.hpp"

#ifdef  VMMLIB_CUBLAS

#include <vmmlib/cublas_dgemm.cu>

namespace vmml
{
	
	bool
	cublas_dgemm_test::run()
	{
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
		
		ok = D == D_check;
		
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
		
		
		
		return ok;
	}
	
	
	
} // namespace vmml
#endif
