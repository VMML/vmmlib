#include "t3_hosvd_test.hpp"
#include "vmmlib/t3_hosvd.hpp"

#include <sstream>

namespace vmml
{
	
	bool
	t3_hosvd_test::run()
	{
		bool ok = false;
		
		double precision = 0.001;


		// hosvd test data after lathauwer et al. 2000a
		tensor3< 3, 3, 3, double > t3_data_hosvd;
		double data_hosvd[] = { 
			0.9073, 1.7842, 2.1236, 0.8924, 1.7753, -0.6631, 2.1488, 4.2495, 1.8260, 
			0.7158, 1.6970, -0.0704, -0.4898, -1.5077, 1.9103, 0.3054, 0.3207, 2.1335, 
			-0.3698, 0.0151, 1.4429, 2.4288, 4.0337, -1.7495, 2.3753, 4.7146, -0.2716
		};
		t3_data_hosvd.set(data_hosvd, data_hosvd + 27);

		//prepare control data
		matrix< 3, 3, double > u1_hosvd;
		matrix< 3, 3, double > u2_hosvd;
		matrix< 3, 3, double > u3_hosvd;
		matrix< 3, 3, double > u1_hosvd_check;
		matrix< 3, 3, double > u2_hosvd_check;
		matrix< 3, 3, double > u3_hosvd_check;
		double data_u1_hosvd[] = { -0.1122204303093513, 0.7738508263292329, -0.6233469929427331, 
			-0.5770542482332139, -0.5614452683822977, -0.5931168562830853, 
			-0.8089591890924934, 0.293145201565186, 0.5095595364450776};
		u1_hosvd_check.set( data_u1_hosvd, data_u1_hosvd + 9);
		double data_u2_hosvd[] = { -0.4624061014026708, 0.01022823125054092, 0.8866092604248326,
			-0.8866395226347692, -0.01338403108646822, -0.4622674816746957,
			0.007138227200902682, -0.9998581154330659, 0.01525761170912346 };
		u2_hosvd_check.set( data_u2_hosvd, data_u2_hosvd + 9);
		double data_u3_hosvd[] = { -0.6208224807897551, -0.498575577624314, -0.6049808598151153,
			0.05741499780009954, -0.7985513947597073, 0.5991820991609654,
			-0.7818458705996144, 0.3372507425105363, 0.5243843736253236 };
		u3_hosvd_check.set( data_u3_hosvd, data_u3_hosvd + 9);
		
		
		///hosvd
		t3_hosvd< 3, 3, 3, 3, 3, 3, double >::hosvd( t3_data_hosvd, u1_hosvd, u2_hosvd, u3_hosvd );
		
		if ( u1_hosvd.equals( u1_hosvd_check, precision ) && u2_hosvd.equals( u2_hosvd_check, precision ) && u3_hosvd.equals( u3_hosvd_check, precision ))
		{	
			log( "HOSVD (3x3x3) compute factor matrices U1, U2, U3", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOSVD (3x3x3): " << std::endl
			<< "U1 should be: " << std::endl << u1_hosvd_check << std::endl
			<< "U1 is: " << std::endl << u1_hosvd << std::endl
			<< "U2 should be: " << std::endl << u2_hosvd_check << std::endl
			<< "U2 is: " << std::endl << u2_hosvd << std::endl
			<< "U3 should be: " << std::endl << u3_hosvd_check << std::endl
			<< "U3 is: " << std::endl << u3_hosvd << std::endl;
			
			
			log_error( error.str() );
		}
		
		
		//hosvd with test data from Lathauwer et al. 2000b (non-cubic)
		matrix< 3, 3, double > u1_2;
		matrix< 2, 2, double > u2_2;
		matrix< 2, 2, double > u3_2;
		matrix< 3, 3, double > u1_2_check;
		matrix< 2, 2, double > u2_2_check;
		matrix< 2, 2, double > u3_2_check;
		double data_u1_2[] = {-0.246452, 0.499323, 0.830625, 0.521727, -0.653918, 0.547898, 0.816738, 0.56839, -0.0993512};
		u1_2_check.set( data_u1_2, data_u1_2 + 9);
		double data_u2_2[] = {0.171472, 0.985189, 0.985189, -0.171472};
		u2_2_check.set( data_u2_2, data_u2_2 + 4);
		double data_u3_2[] = {-0.510464, 0.859899, 0.859899, 0.510464};
		u3_2_check.set( data_u3_2, data_u3_2 + 4);
		
		double data_2[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6};
		tensor3< 3, 2, 2, double > t3_data_2;
		t3_data_2.set(data_2, data_2 + 12);
		
		///hosvd
		t3_hosvd< 3, 2, 2, 3, 2, 2, double >::hosvd( t3_data_2, u1_2, u2_2, u3_2 );
		
		ok = u1_2.equals( u1_2_check, precision );
		ok = ok && (u2_2.equals( u2_2_check, precision ));
		ok = ok && (u3_2.equals( u3_2_check, precision ));
		
		if ( ok )
		{	
			log( "HOSVD (3x2x2) compute factor matrices U1, U2, U3", ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOSVD (3x2x2): " << std::endl
			<< "U1 should be: " << std::endl << u1_2_check << std::endl
			<< "U1 is: " << std::endl << u1_2 << std::endl
			<< "U2 should be: " << std::endl << u2_2_check << std::endl
			<< "U2 is: " << std::endl << u2_2 << std::endl
			<< "U3 should be: " << std::endl << u3_2_check << std::endl
			<< "U3 is: " << std::endl << u3_2 << std::endl;
			
			
			log_error( error.str() );
		}		
		
		
		
		//hoeigs test
		u1_hosvd.zero(); u2_hosvd.zero(); u3_hosvd.zero();
		matrix<3, 3, double> u1_hoeigs_check;
		matrix<3, 3, double> u2_hoeigs_check;
		matrix<3, 3, double> u3_hoeigs_check;
		double data_u1_hoeigs[] = { 0.1122204363, 0.7738507986, -0.623347044, 
			0.5770543218, -0.5614452958, -0.5931168199,
			0.8089591861, 0.2931452692, 0.5095595717};
		u1_hoeigs_check.set( data_u1_hoeigs, data_u1_hoeigs + 9);
		double data_u2_hoeigs[] = { -0.4624060392, -0.01022823341, 0.8866092563,
			-0.886639595, 0.01338403113, -0.4622673988,
			0.007138226647, 0.9998582006, 0.01525761187 };
		u2_hoeigs_check.set( data_u2_hoeigs, data_u2_hoeigs + 9);
		double data_u3_hoeigs[] = { 0.6208223701, 0.4985756576, 0.6049809456,
			-0.05741500854, 0.7985514402, -0.5991820693,
			0.7818458676, -0.3372507095, -0.524384439 };
		u3_hoeigs_check.set( data_u3_hoeigs, data_u3_hoeigs + 9);
		
		t3_hosvd< 3, 3, 3, 3, 3, 3, double >::hoeigs( t3_data_hosvd, u1_hosvd, u2_hosvd, u3_hosvd );
		
		if ( u1_hosvd.equals( u1_hoeigs_check, precision ) && u2_hosvd.equals( u2_hoeigs_check, precision ) && u3_hosvd.equals( u3_hoeigs_check, precision ))
		{	
			log( "HOEIGS compute factor matrices U1, U2, U3", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOEIGS: " << std::endl
			<< "U1 should be: " << std::endl << u1_hoeigs_check << std::endl
			<< "U1 is: " << std::endl << u1_hosvd << std::endl
			<< "U2 should be: " << std::endl << u2_hoeigs_check << std::endl
			<< "U2 is: " << std::endl << u2_hosvd << std::endl
			<< "U3 should be: " << std::endl << u3_hoeigs_check << std::endl
			<< "U3 is: " << std::endl << u3_hosvd << std::endl;
			
			log_error( error.str() );
		}
		

		return ok;
	}

} //end vmml namespace