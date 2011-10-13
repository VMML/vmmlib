#include "qtucker3_tensor_test.hpp"
#include <vmmlib/qtucker3_tensor.hpp>
#include <sstream>

namespace vmml
{
	
	bool
	qtucker3_tensor_test::run()
	{
		bool ok = false;
		double precision = 0.001;
		
		typedef double T_value;
		typedef double T_coeff;
		typedef tensor3< 3, 2, 2, T_value> t3_type;
		typedef matrix< 3, 2, T_coeff > u1_type;
		typedef matrix< 2, 2, T_coeff > u2_type;
		typedef matrix< 2, 1, T_coeff > u3_type;
		typedef tensor3< 2, 2, 1, T_coeff > core_type;
		typedef t3_hooi< 2, 2, 1, 3, 2, 2, float > hooi_type;
		typedef qtucker3_tensor< 2, 2, 1, 3, 2, 2, T_value, T_coeff > tuck3_type;
		
		
		//test data from lathauwer et al. 2000b paper (same test as in t3_hooi_test
		T_value data_als[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6 };
		t3_type t3_data_als;
		t3_data_als.set( data_als, data_als + 12 );
		
		T_coeff data_u1_hooi[] = { 
			-0.2789474111071824, 0.4141266306147135, 
			0.5983607967045262, 0.7806355076145295, 
			0.7511009910815754, -0.4680890279285661 }; //original from paper (u1): {-0.2789, -0.4141, 0.5984, -0.7806, 0.7511, 0.4681};
		u1_type u1_hooi_check;
		u1_hooi_check.set( data_u1_hooi, data_u1_hooi + 6);
		
		T_coeff data_u2_hooi[] = { 
			0.09816424894941811, 0.9951702267593202, 
			0.9951702267593202, -0.098164248949418 }; //original in paper (u2): 0.0982, -0.9952, 0.9952, 0.0982};
		u2_type u2_hooi_check;
		u2_hooi_check.set( data_u2_hooi, data_u2_hooi + 4);
		
		T_coeff data_u3_hooi[] = {-0.5104644303570166, 0.8598988692516616};//original in paper (u3): {0.5105, -0.8599};
		u3_type u3_hooi_check;
		u3_hooi_check.set( data_u3_hooi, data_u3_hooi + 2);
		
		T_coeff data_core_hooi[] = { -10.14733447424582, 0.0, 0.0, -2.760705584847321 };
		core_type core_hooi_check;
		core_hooi_check.set( data_core_hooi, data_core_hooi + 4);
		//end fill test data
		
		u1_type u1_hooi; u1_hooi.zero();
		u2_type u2_hooi; u2_hooi.zero();
		u3_type u3_hooi; u3_hooi.zero();
		core_type core_hooi;
		
		tuck3_type tuck3_hooi;
		tuck3_hooi.tucker_als( t3_data_als, hooi_type::init_hosvd() );
		tuck3_hooi.get_u1( u1_hooi );
		tuck3_hooi.get_u2( u2_hooi );
		tuck3_hooi.get_u3( u3_hooi );
		tuck3_hooi.get_core( core_hooi );
		
		ok = u1_hooi.equals( u1_hooi_check, precision );
		ok = ok && u2_hooi.equals( u2_hooi_check, precision );
		ok = ok && u3_hooi.equals( u3_hooi_check, precision );
		ok = ok && core_hooi.equals( core_hooi_check, precision);
		
		if ( ok )
		{	
			log( "Tucker ALS: rank-(2,2,1) approximation (same test as T3_HOOI)" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker ALS: rank-(2,2,1) approximation: " << std::setprecision(16) << std::endl
			<< "U1 should be: " << std::endl << u1_hooi_check << std::endl
			<< "U1 is: " << std::endl << u1_hooi << std::endl
			<< "U2 should be: " << std::endl << u2_hooi_check << std::endl
			<< "U2 is: " << std::endl << u2_hooi << std::endl
			<< "U3 should be: " << std::endl << u3_hooi_check << std::endl
			<< "U3 is: " << std::endl << u3_hooi << std::endl
			<< "core should be: " << std::endl << core_hooi_check << std::endl
			<< "core is: " << std::endl << core_hooi << std::endl;
			
			
			log_error( error.str() );
		}
		
		//number of nonzeros
		
		size_t number_nonzeros = tuck3_hooi.nnz( );
		size_t number_nonzeros2 = tuck3_hooi.nnz( 0.1 );		
		ok = ( number_nonzeros == 16 ) && (number_nonzeros2 == 12);
		log( "get number of nonzeros" , ok  );
		
		//quantization
		typedef unsigned char T_value_2;
		typedef unsigned short T_coeff_2;
		typedef tensor3< 3, 2, 2, T_value_2 > t3q_type;
		typedef qtucker3_tensor< 2, 2, 2, 3, 2, 2, T_value_2, T_coeff_2 > tuck3q_type;
		typedef t3_hooi< 2, 2, 2, 3, 2, 2, float > hooi_type1;
		
		
		//fill test data
		T_value_2 data_hooi_3[] = { 0, 13, 122, 123, 124, 95, 10, 40, 25, 54, 33, 76};
		t3q_type t3_data_hooi_3;
		t3_data_hooi_3.set(data_hooi_3, data_hooi_3 + 12);
		
		t3q_type t3_data_hooi_3_reco;
		float u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max;

		tuck3q_type tuck3_hooi_3;
		tuck3_hooi_3.enable_quantify_linear();
		tuck3_hooi_3.decompose( t3_data_hooi_3, u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max, hooi_type1::init_hosvd() );
		tuck3_hooi_3.reconstruct( t3_data_hooi_3_reco, u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max );
		double rmse = t3_data_hooi_3_reco.rmse( t3_data_hooi_3 );
		double rmse_check = 5.392896562454479; 
		
		float u_min, u_max;
		tuck3_hooi_3.decompose( t3_data_hooi_3, u_min, u_max, core_min, core_max, hooi_type1::init_hosvd() );
		tuck3_hooi_3.reconstruct( t3_data_hooi_3_reco, u_min, u_max, core_min, core_max );
		double rmse2 = t3_data_hooi_3_reco.rmse( t3_data_hooi_3 );
		
		
		if ( (rmse == rmse_check) && (rmse != 0) && (rmse2 == rmse_check))
		{	
			log( "quantized Tucker ALS ank-(2,2,1) approximation" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "quantized Tucker ALS rank-(2,2,1) approximation: " << std::setprecision(16) << std::endl
			<< "RMSE should be: " << rmse_check << ", is: " << rmse << std::endl
			<< "Tucker3 is : " << std::endl << tuck3_hooi_3 << std::endl;
			
			
			log_error( error.str() );
		}		
		
		
		//tucker3 reconstruction
		typedef int T_value_3;
		typedef int T_coeff_3;
		typedef tensor3< 6, 7, 5, T_value_3 > t3r_type;
		typedef tensor3< 2, 3, 4, T_coeff_3 > t3r_core_type;
		typedef matrix< 6, 2, T_coeff_3 > u1r_type;
		typedef matrix< 7, 3, T_coeff_3 > u2r_type;
		typedef matrix< 5, 4, T_coeff_3 > u3r_type;
		typedef qtucker3_tensor< 2, 3, 4, 6, 7, 5, T_value_3, T_coeff_3 > tuck3r_type;
		
		
		t3r_core_type core; core.fill_increasing_values();
		u1r_type u1; u1.fill(2);
		u2r_type u2; u2.fill(3);
		u3r_type u3; u3.fill(1);
		t3r_type t3_reco;
		
		tuck3r_type tuck3( core, u1, u2, u3 );
		tuck3.reconstruct( t3_reco, 2, 2, 3, 3, 1, 1, 0, 23 );
		
		t3r_type t3_reco_check;
		t3_reco_check.fill(1656);
		
	
		ok = t3_reco_check == t3_reco ;
		log( "tucker3 reconstruction", ok );
		
		
		
		return ok;
	}
	
	
	
} // namespace vmml

