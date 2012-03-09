#include "blas_daxpy_test.hpp"

#include <vmmlib/blas_daxpy.hpp>

namespace vmml
{
	
	bool
	blas_daxpy_test::run()
	{
		bool ok = false;
		
		double a  = 0.8147;
		vector< 2, double > B;
		vector< 2, double > C;
		vector< 2, double > C_check;
		
		double BData[] = {  0.9649, 0.1576 };
		B = BData;
		double CData[] = {  0.786104, 0.128397 };
		C_check = CData;
		
		
		blas_daxpy< 2, double > blas_daxpy1;
		blas_daxpy1.compute( a, B, C );
		
		ok = C.equals( C_check, 0.0001 );
		
		log( "compute single daxpy", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "daxpy product of \n" << a << "\n"
			<< "and \n" << B << "\n"
			<< "should be\n" << C_check << "\n"
            << "is\n" << C << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		
		
		//use multiple daxpys for final 
		
		matrix< 3,3, float > right_m;
		matrix< 2,3, float > res_m;
		matrix< 2,3, float > res_m_check;
		matrix< 2,3, float > left_m;
		
		float lData[] = {  0.9649  ,  0.9706 ,   0.4854, 0.1576  ,  0.9572 ,   0.8003 };
		left_m = lData;
		float rData[] = {  0.8147,    0.9134 ,   0.2785, 0.9058 ,   0.6324  ,  0.5469, 0.1270 ,   0.0975  ,  0.9575 };
		right_m = rData;
		float res_data[] = {  1.7269,1.5425,1.2643, 	1.0971,0.8273,1.3337 };
		res_m_check = res_data;
		
		blas_daxpy< 2, float > blas_daxpy2;
		blas_daxpy2.compute_mmm( left_m, right_m, res_m );
		
		ok = res_m.equals( res_m_check, 0.0001 );
		
		log( "compute matrix-matrix multiplication by multiple daxpy's", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "left matrix \n" << left_m << "\n"
			<< "right matrix \n" << right_m << "\n"
			<< "result matrix \n"
			<< "should be\n" << res_m_check << "\n"
            << "is\n" << res_m << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		
		
	
		return ok;
	}
	
	
	
} // namespace vmml

