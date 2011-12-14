#include "t3_hooi_test.hpp"
#include "vmmlib/t3_hooi.hpp"

#include <sstream>

namespace vmml
{
	
	bool
	t3_hooi_test::run()
	{
		bool ok = false;
		
		double precision = 0.001;

		//test data from lathauwer et al. 2000b paper
		tensor3< 3, 2, 2, double> t3_data;
		double data[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6 };
		t3_data.set( data, data + 12 );

		matrix< 3, 2, double > u1; u1.zero();
		matrix< 2, 2, double > u2; u2.zero();
		matrix< 2, 1, double > u3; u3.zero();
		matrix< 3, 2, double > u1_check;
		matrix< 2, 2, double > u2_check;
		matrix< 2, 1, double > u3_check;
		double data_u1[] = { 
			-0.2789474111071824, 0.4141266306147135, 
			0.5983607967045262, 0.7806355076145295, 
			0.7511009910815754, -0.4680890279285661 }; //original from paper (u1): {-0.2789, -0.4141, 0.5984, -0.7806, 0.7511, 0.4681};
		u1_check.set( data_u1, data_u1 + 6);
		double data_u2[] = { 
			0.09816424894941811, 0.9951702267593202, 
			0.9951702267593202, -0.098164248949418 }; //original in paper (u2): 0.0982, -0.9952, 0.9952, 0.0982};
		u2_check.set( data_u2, data_u2 + 4);
		double data_u3[] = { -0.5104644303570166, 0.8598988692516616 };//original in paper (u3): {0.5105, -0.8599};
		u3_check.set( data_u3, data_u3 + 2);
		
		tensor3< 2, 2, 1, double > core;
		tensor3< 2, 2, 1, double > core_check;
		double data_core[] = { -10.14733447424582, 0.0, 0.0, -2.760705584847321 };
		core_check.set( data_core, data_core + 4);
		
		typedef t3_hooi< 2, 2, 1, 3, 2, 2, double > hooi_type;
		hooi_type::als( t3_data, u1, u2, u3, core, hooi_type::init_hosvd() );
		
		ok = u1.equals( u1_check, precision );
		ok = ok && ( u2.equals( u2_check, precision ) );
		ok = ok && ( u3.equals( u3_check, precision ) );
		ok = ok && ( u3.equals( u3_check, precision ) );
		ok = ok && ( core.equals( core_check, precision ) );
		
		if ( ok )
		{	
			log( "HOOI rank-(2,2,1) approximation" , ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOOI rank-(2,2,1) approximation: " << std::setprecision(16) << std::endl
			<< "U1 should be: " << std::endl << u1_check << std::endl
			<< "U1 is: " << std::endl << u1 << std::endl
			<< "U2 should be: " << std::endl << u2_check << std::endl
			<< "U2 is: " << std::endl << u2 << std::endl
			<< "U3 should be: " << std::endl << u3_check << std::endl
			<< "U3 is: " << std::endl << u3 << std::endl
			<< "core should be: " << std::endl << core_check << std::endl
			<< "core is: " << std::endl << core << std::endl;
			
			
			log_error( error.str() );
		}
		
	
		//(1a) derive core tensor (with pseudo inverse)
		core.zero();
		t3_hooi< 2, 2, 1, 3, 2, 2, double >::derive_core( t3_data, u1_check, u2_check, u3_check, core );
		
		if ( core.equals( core_check, precision ))
		{	
			log( "derive core tensor", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "derive core: " << std::endl
			<< "core should be: " << std::endl << core_check << std::endl
			<< "core is: " << std::setprecision(16) << std::endl << core << std::endl;
			
			log_error( error.str() );
		}
		
		
		//(1b) derive core tensor with orthogonal basis
		core.zero();
		t3_hooi< 2, 2, 1, 3, 2, 2, double >::derive_core_orthogonal_bases( t3_data, u1_check, u2_check, u3_check, core );
		
		if ( core.equals( core_check, precision ))
		{	
			log( "derive core tensor (orthogonal bases)", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "derive core (orthogonal bases): " << std::endl
			<< "core should be: " << std::endl << core_check << std::endl
			<< "core is: " << std::setprecision(16) << std::endl << core << std::endl;
			
			log_error( error.str() );
		}
		
		return ok;
	}
	
} //end vmml namespace

