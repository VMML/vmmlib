#include "tucker3_tensor_test.hpp"

#include <vmmlib/tucker3_tensor.hpp>

#include <sstream>

namespace vmml
{
	
	bool
	tucker3_tensor_test::run()
	{
		bool ok = false;
		
		tensor3< 2, 3, 4, uint16_t >  core;
		core.fill_increasing_values();
		tensor3< 6, 7, 5, uint16_t >  t3_reco;
		matrix<6, 2, uint16_t> u1;
		matrix<7, 3, uint16_t> u2;
		matrix<5, 4, uint16_t> u3;
		u1.fill(2);
		u2.fill(3);
		u3.fill(1);
		
		tucker3_tensor<2, 3, 4, 6, 7, 5, uint16_t > tuck3( core, u1, u2, u3 );
		
		//tucker3 reconstruction
		tuck3.reconstruction( t3_reco );
		
		tensor3<6, 7, 5, uint16_t> t3_reco_test;
		t3_reco_test.fill(1656);
		//std::cout << "Tucker3 reconstruction (all values should be 1656): " << std::endl << t3_reco << std::endl;
		//std::cout << "Tucker3 core : " << std::endl << tuck3.get_core() << std::endl;
		
		if ( t3_reco_test == t3_reco)
		{	
			log( "tucker3 reconstruction", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker3 reconstruction (all values should be 1656): " << std::endl << t3_reco
			<< std::endl;
			log_error( error.str() );
		}
		

		//decomposition
		//tucker method1
		tensor3<6, 7, 5, uint16_t> t3_2;
		t3_2.fill_increasing_values();
		matrix<6, 6, uint16_t> u1_2;
		matrix<7, 7, uint16_t> u2_2;
		matrix<5, 5, uint16_t> u3_2;
		tucker3_tensor<6, 7, 5, 6, 7, 5, uint16_t > tuck3_deco( t3_2, u1_2, u2_2, u3_2 );
		tensor3<6, 7, 5, uint16_t> t3_deco;
		tuck3_deco.decomposition( t3_deco );
		

		
		//rank reduction
		tensor3< 1, 2, 3, uint16_t >  core_red;
		matrix<6, 1, uint16_t> u1_red;
		matrix<7, 2, uint16_t> u2_red;
		matrix<5, 3, uint16_t> u3_red;
		
		tucker3_tensor< 1, 2, 3, 6, 7, 5, uint16_t > tuck3_red( core_red, u1_red, u2_red, u3_red );
		tuck3_red.rank_reduction( tuck3 );
		
		u1_red.fill(2);
		u2_red.fill(3);
		u3_red.fill(1);
		double data[] = { 0, 1, 6, 7, 12, 13 };
		core_red.set(data, data+6);
		
		if ( tuck3_red.get_u1() == u1_red && tuck3_red.get_u2() == u2_red && tuck3_red.get_u3() == u3_red && tuck3_red.get_core() == core_red)
		{	
			log( "tucker3 rank reduction", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker3 rank reduction: " << std::endl
			<< "u1 should be: " << u1_red << std::endl
			<< "u1 is: " << tuck3_red.get_u1() << std::endl
			<< "u2 should be: " <<  u2_red << std::endl
			<< "u2 is: " <<  tuck3_red.get_u2() << std::endl
			<< "u3 should be: " << u3_red << std::endl
			<< "u3 is: " << tuck3_red.get_u3() << std::endl
			<< "core should be: " << core_red << std::endl
			<< "core is: " << tuck3_red.get_core() << std::endl;

			log_error( error.str() );
		}
		
		
		//basis matrices subsampling
		tensor3< 3, 4, 3, uint16_t > t3_sub;
		tensor3< 2, 3, 4, uint16_t > core_sub;
		matrix< 3, 2, uint16_t > u1_sub;
		matrix< 4, 3, uint16_t > u2_sub;
		matrix< 3, 4, uint16_t > u3_sub;
		tucker3_tensor< 2, 3, 4, 3, 4, 3, uint16_t > tuck3_sub( core_sub, u1_sub, u2_sub, u3_sub );
		
		tuck3_sub.subsampling( tuck3, 2);
		tuck3_sub.reconstruction( t3_sub );
		
		tensor3< 3, 4, 3, uint16_t > t3_sub_test;
		t3_sub_test.fill(1656);
		if ( t3_sub_test == t3_sub )
		{	
			log( "basis matrices subsampling", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "basis matrices subsampling with factor 2: " << std::endl << t3_sub
			<< std::endl;
			log_error( error.str() );
		}
		
		//basis matrices region of interest selection
		tensor3< 1, 1, 3, uint16_t > t3_roi;
		tensor3< 2, 3, 4, uint16_t > core_roi;
		matrix< 1, 2, uint16_t > u1_roi;
		matrix< 1, 3, uint16_t > u2_roi;
		matrix< 3, 4, uint16_t > u3_roi;
		tucker3_tensor< 2, 3, 4, 1, 1, 3, uint16_t > tuck3_roi( core_roi, u1_roi, u2_roi, u3_roi );
		
		tuck3_roi.region_of_interest( tuck3, 0, 1, 1, 2, 1, 4);
		tuck3_roi.reconstruction( t3_roi );

		tensor3< 1, 1, 3, uint16_t > t3_roi_test;
		t3_roi_test.fill(1656);
		if ( t3_roi_test == t3_roi)
		{	
			log( "basis matrices region of interest selection", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "basis matrices region of interest selection: " << std::endl << t3_roi
			<< std::endl;
			log_error( error.str() );
		}
		
		
	
		ok = true;
		return ok;
	}
	
	
	
} // namespace vmml

