#include "tucker3_tensor_test.hpp"

#include <vmmlib/tucker3_tensor.hpp>
#include <sstream>

namespace vmml
{
	
	bool
	tucker3_tensor_test::run()
	{
		bool ok = false;
		double precision = 0.001;
		
		//test data from lathauwer et al. 2000b paper (same test as in t3_hooi_test
		tensor3< 3, 2, 2, double> t3_data_als;
		double data_als[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6 };
		t3_data_als.set( data_als, data_als + 12 );
		
		matrix<3, 2, double> u1_hooi; u1_hooi.zero();
		matrix<2, 2, double> u2_hooi; u2_hooi.zero();
		matrix<2, 1, double> u3_hooi; u3_hooi.zero();
		matrix<3, 2, double> u1_hooi_check;
		matrix<2, 2, double> u2_hooi_check;
		matrix<2, 1, double> u3_hooi_check;
		double data_u1_hooi[] = { -0.2789474111071824, -0.4141266306147135, 0.5983607967045262, -0.7806355076145295, 0.7511009910815754, 0.4680890279285661}; //original from paper (u1): {-0.2789, -0.4141, 0.5984, -0.7806, 0.7511, 0.4681};
		u1_hooi_check.set( data_u1_hooi, data_u1_hooi + 6);
		double data_u2_hooi[] = { -0.09816424894941811, -0.9951702267593202, -0.9951702267593202, 0.098164248949418}; //original in paper (u2): 0.0982, -0.9952, 0.9952, 0.0982};
		u2_hooi_check.set( data_u2_hooi, data_u2_hooi + 4);
		double data_u3_hooi[] = {-0.5104644303570166, 0.8598988692516616};//original in paper (u3): {0.5105, -0.8599};
		u3_hooi_check.set( data_u3_hooi, data_u3_hooi + 2);
		
		tensor3< 2, 2, 1, double > core_hooi;
		tensor3< 2, 2, 1, double > core_hooi_check;
		double data_core_hooi[] = { 10.14733447424582, 0.0, 0.0, -2.760705584847321 };
		core_hooi_check.set( data_core_hooi, data_core_hooi + 4);
		
		tucker3_tensor< 2, 2, 1, 3, 2, 2, double, double > tuck3_hooi;
		
		tuck3_hooi.tucker_als( t3_data_als );
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
		tensor3< 3, 2, 2, unsigned char> t3_data_hooi_3;
		unsigned char data_hooi_3[] = { 0, 13, 122, 123, 124, 95, 10, 40, 25, 54, 33, 76};
		t3_data_hooi_3.set(data_hooi_3, data_hooi_3 + 12);
		
		tucker3_tensor< 2, 2, 2, 3, 2, 2, unsigned char, unsigned short > tuck3_hooi_3;
		
		float u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max, u_min, u_max;
		
		//std::cout << "start quant" << std::endl;
		tuck3_hooi_3.enable_quantify_linear();
		tuck3_hooi_3.decompose( t3_data_hooi_3, u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max );
		//std::cout << "Tucker3 is : " << std::endl << tuck3_hooi_3 << std::endl;
		
		/*std::cout
		<< "u1_min: " << u1_min << ", u1_max: " << u1_max << std::endl
		<< "u2_min: " << u2_min << ", u2_max: " << u2_max << std::endl
		<< "u3_min: " << u3_min << ", u3_max: " << u3_max << std::endl;*/
		
		tensor3< 3, 2, 2, unsigned char > t3_data_hooi_3_reco;
		tuck3_hooi_3.reconstruct( t3_data_hooi_3_reco, u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max );
		double rmse = t3_data_hooi_3_reco.rmse( t3_data_hooi_3 );
		double rmse_check = 5.392896562454479; 
		//std::cout << "reco : " << std::endl << t3_data_hooi_3_reco << std::endl;
		
		tuck3_hooi_3.decompose( t3_data_hooi_3, u_min, u_max, core_min, core_max );
		//std::cout << "Tucker3 is (1 u min/max) : " << std::endl << tuck3_hooi_3 << std::endl;
		//std::cout
		//<< "u_min: " << u_min << ", u_max: " << u_max << std::endl;

		tuck3_hooi_3.reconstruct( t3_data_hooi_3_reco, u_min, u_max, core_min, core_max );
		//std::cout << "reco : " << std::endl << t3_data_hooi_3_reco << std::endl;
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
		
		
		
		//export
		std::vector< float > export_data;
		tuck3_hooi.export_to( export_data );
		
		double export_data_check[] = {
			-0.2789473533630371, 0.5983607769012451, 0.751101016998291, -0.4141266345977783, -0.7806354761123657, 0.4680890440940857, 
			-0.09816429018974304, -0.9951702356338501, -0.9951702356338501, 0.09816429018974304,
			-0.5104644894599915, 0.8598988652229309,
			10.14733409881592, 1.639127731323242e-07, -2.384185791015625e-07, -2.76070547103881};
		
		ok = true;
		for (int i = 0; (i < 16) && ok; ++i )
		{
			//std::cout << std::setprecision( 16) << export_data[i] << std::endl;
			ok = (abs(export_data_check[i] - export_data[i])) < 1.0e-6;
		}

		log( "export tucker3", ok  );
		
		
		//import tucker3 from vector
		std::vector< float > in_data = export_data;
				
		matrix<3, 2, double> u1_imported;
		matrix<2, 2, double> u2_imported;
		matrix<2, 1, double> u3_imported;
		tensor3<2, 2, 1, double> core_imported;
		tucker3_tensor< 2, 2, 1, 3, 2, 2, double, double > tuck3_import;
		
		tuck3_import.import_from( in_data );
		
		tuck3_import.get_core( core_imported ); tuck3_import.get_u1( u1_imported  ); tuck3_import.get_u2( u2_imported ); tuck3_import.get_u3( u3_imported );
		
		if ( u1_imported.equals( u1_hooi_check, precision ) && u2_imported.equals( u2_hooi_check, precision ) && u3_imported.equals( u3_hooi_check, precision ) && core_imported.equals( core_hooi_check, precision))
		{	
			log( "import tucker3" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "import tucker3: " << std::endl
			<< "U1 should be: " << std::endl << u1_hooi_check << std::endl
			<< "U1 is: " << std::endl << u1_imported << std::endl
			<< "U2 should be: " << std::endl << u2_hooi_check << std::endl
			<< "U2 is: " << std::endl << u2_imported << std::endl
			<< "U3 should be: " << std::endl << u3_hooi_check << std::endl
			<< "U3 is: " << std::endl << u3_imported << std::endl
			<< "core should be: " << std::endl << core_hooi_check << std::endl
			<< "core is: " << std::endl << core_imported << std::endl;
			
			
			log_error( error.str() );
		}		
		
		//export bytes (for one min/max value for all Us)
		std::vector<unsigned char> export_data_vec;
		tuck3_hooi_3.export_quantized_to( export_data_vec );
		
		//check exported data
		double export_data2_check[] = {
			-0.921997, 0.921997, -251.345, 16.3135,  /* all min/max values -0.715833, 0.865947, -0.721095, 0.692817, -0.921997, 0.921997,*/
			29346, 7327, 8188, 63543, 18451, 43301,
			8145, 7140, 7140, 57390,
			0, 19007, 19007, 65535,
			0, 60923, 61001, 59355, 61571, 52960, 51805, 65535 };
		
		ok = true;
		
		float * float_ptr = (float*)&(export_data_vec[0]);
		//check u2 min/max
		float u1_min_e = *float_ptr; float_ptr++;
		float u1_max_e = *float_ptr; float_ptr++;
		ok = ((fabs(float(export_data2_check[0]) - u1_min_e)) < 1.0e-4) && ((fabs(float(export_data2_check[1]) - u1_max_e)) < 1.0e-4 );
		//std::cout<<"#### U1 min value === " << u1_min_e<<", U1 max value === " << u1_max_e << std::endl;

		float core_min_e = *float_ptr; float_ptr++;
		float core_max_e = *float_ptr; float_ptr++;
		ok = ok ? ((fabs(float(export_data2_check[2]) - core_min_e)) < 1.0e-3) && ((fabs(float(export_data2_check[3]) - core_max_e)) < 1.0e-3 ) : ok;
		//std::cout<<"#### core min value === " << core_min_e <<", core max value === " << core_max_e << std::endl;
		//std::cout<<"#### shold be core min value === " << export_data2_check[2] <<", core max value === " << export_data2_check[3] << std::endl;
		
		unsigned short* value_ptr = (unsigned short*)float_ptr;
		unsigned short value;
		size_t index = 4;
		size_t end_index = index + 6;
		//check u1
		for ( ; (index < end_index) && ok; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### U1 value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < 1.0e-6;
		}
		//check u2
		end_index += 4;
		for ( ; (index < end_index) && ok ; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### U2 value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < 1.0e-6;
		}
		//check u3		
		end_index += 4;
		for ( ; (index < end_index) ; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### U3 value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < 1.0e-6;
		}
		
		//check core values
		end_index += 8;
		for ( ; (index < end_index) && ok; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### core value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < 1.0e-6;
		}
		
		log( "export tucker3 (bytes) ", ok  );
		
		
		//import from bytes
		tucker3_tensor< 2, 2, 2, 3, 2, 2, unsigned char, unsigned short > tuck3_import2;
		tuck3_import2.enable_quantify_coeff();
		tuck3_import2.import_quantized_from( export_data_vec );
		
		tensor3< 3, 2, 2, unsigned char > import_reco;
		tuck3_import2.reconstruct( import_reco );
		rmse = import_reco.rmse( t3_data_hooi_3 );
		rmse_check = 5.392896562454479; 

		if ( (rmse <= rmse_check+0.5)&& rmse > 0 )
		{
			log( "import tucker3 (bytes)" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "import tucker3 (bytes): " << std::endl
			<< "RMSE should be: " << rmse_check << ", is: " << rmse << std::endl
			<< "Tucker3 is: " << std::endl << tuck3_import2 << std::endl;

			log_error( error.str() );
		}		
		

		//export bytes (with optimized quantization of core) 
		std::vector<unsigned char> export_data_vec2;
		tuck3_hooi_3.export_hot_quantized_to( export_data_vec2 );
		
		tucker3_tensor< 2, 2, 2, 3, 2, 2, unsigned char, unsigned short > tuck3_import3;
		tuck3_import3.import_hot_quantized_from( export_data_vec2 );
		
		tensor3< 3, 2, 2, unsigned char > import_reco3;
		tuck3_import3.reconstruct( import_reco3 );
		rmse = import_reco3.rmse( t3_data_hooi_3 );
		rmse_check = 5.400617248673217; 
		
		ok = (rmse <= rmse_check+0.5)&& rmse > 0;
		log( "export/import tucker3 (bytes) with hot core quantization" , ok  );
			
		//export bytes for TtM 
		std::vector<unsigned char> export_data_vec4;
		tuck3_hooi_3.enable_quantify_log();
		tuck3_hooi_3.decompose( t3_data_hooi_3, u_min, u_max, core_min, core_max );
		tuck3_hooi_3.export_ttm_quantized_to( export_data_vec4 );
				
		tucker3_tensor< 2, 2, 2, 3, 2, 2, unsigned char, unsigned short > tuck3_import4;
		tuck3_import4.import_ttm_quantized_from( export_data_vec4 );
				
		tensor3< 3, 2, 2, unsigned char > import_reco4;
		tuck3_import4.reconstruct( import_reco4 );
		rmse = import_reco4.rmse( t3_data_hooi_3 );
		rmse_check = 5.400617248673217; 
				
		ok = (rmse <= rmse_check+0.5)&& rmse > 0;
		log( "export/import TtM tucker3 (bytes)" , ok  );

		
		//tucker3 reconstruction
		tensor3< 2, 3, 4, int >  core;
		core.fill_increasing_values();
		tensor3< 6, 7, 5, int >  t3_reco;
		matrix<6, 2, int> u1;
		matrix<7, 3, int> u2;
		matrix<5, 4, int> u3;
		u1.fill(2);
		u2.fill(3);
		u3.fill(1);
		
		tucker3_tensor<2, 3, 4, 6, 7, 5, int, int > tuck3( core, u1, u2, u3 );

		tuck3.reconstruct( t3_reco );
		
		tensor3<6, 7, 5, int> t3_reco_test;
		t3_reco_test.fill(1656);
		
		//thresholded tucker
		tucker3_tensor<2, 3, 4, 6, 7, 5, int, int > tuck3_copy( tuck3 );
		//std::cout << "tucker3: " << std::endl << tuck3_copy << std::endl;

		tensor3< 6, 7, 5, int >  t3_reco_thresh1;
		tensor3< 6, 7, 5, int >  t3_reco_thresh2;
		size_t nnz_core = 0;
		tuck3_copy.threshold_core( 4, nnz_core );
		tuck3_copy.reconstruct( t3_reco_thresh1 );
		tuck3_copy.threshold_core( 12, nnz_core );
		tuck3_copy.reconstruct( t3_reco_thresh2 );
		
		tensor3<6, 7, 5, int> t3_reco_thresh1_check;
		t3_reco_thresh1_check.fill(1596);
		tensor3<6, 7, 5, int> t3_reco_thresh2_check;
		t3_reco_thresh2_check.fill(1188);
		
		ok = ( t3_reco_test == t3_reco ) && ( t3_reco_thresh1 == t3_reco_thresh1_check ) && ( t3_reco_thresh2 == t3_reco_thresh2_check);
		log( "tucker3 reconstruction (incl. core thresholding)", ok );
		

		//rank reduction
		tensor3< 1, 2, 3, int >  core_red;
		matrix<6, 1, int> u1_red;
		matrix<7, 2, int> u2_red;
		matrix<5, 3, int> u3_red;
		tensor3< 1, 2, 3, int >  core_red2;
		matrix<6, 1, int> u1_red2;
		matrix<7, 2, int> u2_red2;
		matrix<5, 3, int> u3_red2;
		
		u1_red.fill(2);
		u2_red.fill(3);
		u3_red.fill(1);
		double data[] = { 0, 1, 6, 7, 12, 13 };
		core_red.set(data, data+6);		
		
		tucker3_tensor< 1, 2, 3, 6, 7, 5, int, int > tuck3_red( core_red, u1_red, u2_red, u3_red );
		tuck3_red.reduce_ranks( tuck3 );
				
		tuck3_red.get_u1( u1_red2 ); tuck3_red.get_u2( u2_red2 ); tuck3_red.get_u3( u3_red2 ); tuck3_red.get_core( core_red2 );
		
		if (  u1_red2 == u1_red && u2_red2 == u2_red && u3_red2 == u3_red && core_red2 == core_red)
		{	
			log( "tucker3 reduce ranks", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker3 reduce ranks: " << std::endl
			<< "u1 should be: " << u1_red << std::endl
			<< "u1 is: " << u1_red2 << std::endl
			<< "u2 should be: " <<  u2_red << std::endl
			<< "u2 is: " <<  u2_red2 << std::endl
			<< "u3 should be: " << u3_red << std::endl
			<< "u3 is: " << u3_red2 << std::endl
			<< "core should be: " << core_red << std::endl
			<< "core is: " << core_red2 << std::endl;

			log_error( error.str() );
		}
		
		
		//factor matrices subsampling
		tensor3< 3, 4, 3, int > t3_sub;
		tensor3< 2, 3, 4, int > core_sub;
		matrix< 3, 2, int > u1_sub;
		matrix< 4, 3, int > u2_sub;
		matrix< 3, 4, int > u3_sub;
		tucker3_tensor< 2, 3, 4, 3, 4, 3, int, int > tuck3_sub( core_sub, u1_sub, u2_sub, u3_sub );
		
		tuck3_sub.subsampling( tuck3, 2);
		tuck3_sub.reconstruct( t3_sub );
		
		tensor3< 3, 4, 3, int > t3_sub_test;
		t3_sub_test.fill(1656);
		if ( t3_sub_test == t3_sub )
		{	
			log( "factor matrices subsampling", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "factor matrices subsampling with factor 2: " << std::endl << t3_sub
			<< std::endl;
			log_error( error.str() );
		}
		
		//factor matrices subsampling, average data
		tensor3< 3, 4, 3, int > t3_sub_avg;
		tensor3< 2, 3, 4, int > core_sub_avg;
		matrix< 3, 2, int > u1_sub_avg;
		matrix< 4, 3, int > u2_sub_avg;
		matrix< 3, 4, int > u3_sub_avg;
		tucker3_tensor< 2, 3, 4, 3, 4, 3, int, int > tuck3_sub_avg( core_sub_avg, u1_sub_avg, u2_sub_avg, u3_sub_avg );
		
		/*matrix<6, 2, int> u1_new;
		u1_new.fill(5);
		u1_new.at(1,1) = 2;
		u1_new.at(4,0) = 3;
		u1_new.at(3,1) = 8;
		tuck3.set_u1(u1_new);*/
		
		tuck3_sub_avg.subsampling_on_average( tuck3, 2);
		tuck3_sub_avg.reconstruct( t3_sub_avg );
		
		t3_sub_test.fill(1656);
		if ( t3_sub_test == t3_sub_avg )
		{	
			log( "factor matrices subsampling on average", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "factor matrices subsampling on average with factor 2: " << std::endl << t3_sub_avg
			<< std::endl;
			log_error( error.str() );
		}
		
		
		//factor matrices region of interest selection
		tensor3< 1, 1, 3, int > t3_roi;
		tensor3< 2, 3, 4, int > core_roi;
		matrix< 1, 2, int > u1_roi;
		matrix< 1, 3, int > u2_roi;
		matrix< 3, 4, int > u3_roi;
		tucker3_tensor< 2, 3, 4, 1, 1, 3, int, int > tuck3_roi( core_roi, u1_roi, u2_roi, u3_roi );
		
		tuck3_roi.region_of_interest( tuck3, 0, 1, 1, 2, 1, 4);
		tuck3_roi.reconstruct( t3_roi );

		tensor3< 1, 1, 3, int > t3_roi_test;
		t3_roi_test.fill(1656);
		if ( t3_roi_test == t3_roi)
		{	
			log( "factor matrices region of interest selection", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "factor matrices region of interest selection: " << std::endl << t3_roi
			<< std::endl;
			log_error( error.str() );
		}
		
		{
			
			
			
		}
		
		
	
		ok = true;
		return ok;
	}
	
	
	
} // namespace vmml

