#include "t3_ttm_test.hpp"

#include <vmmlib/t3_ttm.hpp>

namespace vmml
{
	
	bool
	t3_ttm_test::run()
	{
		bool ok = false;
		
		tensor3< 2, 3, 4, int >  t3;
		//test tensor3 matrix multiplication
		t3.fill_increasing_values();
		matrix<5, 4, int> u3;
		u3.fill(1);
		tensor3<2, 3, 5, int> t3_jji;
		t3_ttm::multiply_horizontal_bwd(t3, u3, t3_jji );
		
		tensor3<6, 3, 5, int> t3_iji;
		matrix<6, 2, int> u1;
		u1.fill(2);
		t3_ttm::multiply_lateral_bwd(t3_jji, u1, t3_iji );
		
		tensor3<6, 7, 5, int> t3_iii;
		matrix<7, 3, int> u2;
		u2.fill(3);
		t3_ttm::multiply_frontal_bwd(t3_iji, u2, t3_iii );
		
		
		t3.fill_increasing_values();
		tensor3<6, 7, 5, int> t3_reco;
		t3_ttm::full_tensor3_matrix_multiplication( t3, u1, u2, u3, t3_reco ); 
		tensor3<6, 7, 5, int> t3_reco2;
		t3_ttm::full_tensor3_matrix_kronecker_mult( t3, u1, u2, u3, t3_reco2 ); 
		ok = t3_reco == t3_reco2;
		
		tensor3<6, 7, 5, int> t3_iii_test;
		t3_iii_test.fill(1656);
		if ( t3_iii_test == t3_reco && t3_iii_test == t3_iii && ok )
		{	
			log( "tensor3 matrix multiplication along all three modes", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "T3_result (all values should be 1656): " << std::endl << t3_iii
			<< std::endl;
			log_error( error.str() );
		}
		
		
		return ok;
	}
	
	
	
} // namespace vmml

