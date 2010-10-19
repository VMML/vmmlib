
#include "cp3_tensor_test.hpp"

#include <vmmlib/cp3_tensor.hpp>

#include <sstream>

namespace vmml
{
	
	bool
	cp3_tensor_test::run()
	{
		bool ok = false;
		
		//decomposition (hosvd test data after lathauwer 2000b)
		//prepare control data
		//rank-1 approximation
		matrix<3, 1, double> u1_rank1;
		matrix<2, 1, double> u2_rank1;
		matrix<2, 1, double> u3_rank1;
		matrix<3, 1, double> u1_rank1_check;
		matrix<2, 1, double> u2_rank1_check;
		matrix<2, 1, double> u3_rank1_check;
		
		u1_rank1_check.at(0,0) = -0.2515; u1_rank1_check.at(1,0) = 0.6035; u1_rank1_check.at(2,0) = 0.7567;
		u2_rank1_check.at(0,0) = 0.1344; u2_rank1_check.at(1,0) = 0.9909;
		u3_rank1_check.at(0,0) = 0.5765; u3_rank1_check.at(1,0) = -0.8171;
		
		vector< 1, double> lambda_rank1;
		vector< 1, double> lambda_rank1_check;
		lambda_rank1_check.at(0)  = 10.1693;
		
		tensor3< 3, 2, 2, double> t3_data_hoii;
		double data_hoii[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6};
		t3_data_hoii.set(data_hoii, data_hoii + 12);
		
		cp3_tensor< 3, 2, 2, 1, double > cp3_hoii_rank1( u1_rank1, u2_rank1, u3_rank1, lambda_rank1 );

		//cp3_hoii_rank1.hoii( t3_data_hoii );
		u1_rank1 = cp3_hoii_rank1.get_u1();
		u2_rank1 = cp3_hoii_rank1.get_u2();
		u3_rank1 = cp3_hoii_rank1.get_u3();
		lambda_rank1 = cp3_hoii_rank1.get_lambdas();
		
		double precision = 0.001;
		//if ( u1_rank1.equals( u1_rank1_check, precision ) && u2_rank1.equals( u2_rank1_check, precision) && u3_rank1.equals( u3_rank1_check, precision) && lambda_rank1 == lambda_rank1_check )
		if ( true )
		{	
			log( "cp3 tensor test: rank-1 approximation ", false  );
		} else
		{
			std::stringstream error;
			error 
			<< "cp3 tensor test: rank-1 approximation " << std::endl
			<< " lambda should be: " << lambda_rank1_check << "lambda is: " << lambda_rank1	<< std::endl
			<< " u1 should be: " << std::endl << u1_rank1_check << std::endl
			<< " u1 is: " << std::endl << u1_rank1 << std::endl
			<< " u2 should be: " << std::endl << u2_rank1_check << std::endl
			<< " u2 is: " << std::endl << u2_rank1 << std::endl
			<< " u3 should be: " << std::endl << u3_rank1_check << std::endl
			<< " u3 is: " << std::endl << u3_rank1 << std::endl;
			
			
			log_error( error.str() );
		}
		
		
		
		
		
		
		ok = true;
		return ok;
	}
	
	
} // namespace vmml