#include "tucker3_exporter_importer_test.hpp"
#include <vmmlib/tucker3_exporter.hpp>
#include <vmmlib/tucker3_importer.hpp>
#include <sstream>

//@SUS: check hot core quantization, is it correctly implementented

namespace vmml
{
	
	bool
	tucker3_exporter_importer_test::run()
	{
		bool ok = false;
		
		double precision = 1.0e-6;

		typedef tensor3< 3, 2, 2, double> t3_type;
		typedef t3_hooi< 2, 2, 1, 3, 2, 2, float > hooi_type;
		typedef tucker3_tensor< 2, 2, 1, 3, 2, 2, double, double > tuck3_type;
		typedef tucker3_exporter<2, 2, 1, 3, 2, 2, double, double> tuck3_exporter_t;
		typedef tucker3_importer< 2, 2, 1, 3, 2, 2, double, double > tuck3_importer_t;
		typedef matrix< 3, 2, double > u1_type;
		typedef matrix< 2, 2, double > u2_type;
		typedef matrix< 2, 1, double > u3_type;
		typedef tensor3<2, 2, 1, double> core_type;
		typedef std::vector< float > out_vec_type;
		
		//prepare check data
		double data_u1[] = { 
			-0.2789474111071824, 0.4141266306147135, 
			0.5983607967045262, 0.7806355076145295, 
			0.7511009910815754, -0.4680890279285661 }; //original from paper (u1): {-0.2789, -0.4141, 0.5984, -0.7806, 0.7511, 0.4681};
		u1_type u1_check;
		u1_check.set( data_u1, data_u1 + 6);
		
		double data_u2[] = { 
			0.09816424894941811, 0.9951702267593202, 
			0.9951702267593202, -0.098164248949418 }; //original in paper (u2): 0.0982, -0.9952, 0.9952, 0.0982};
		u2_type u2_check;
		u2_check.set( data_u2, data_u2 + 4);
		
		double data_u3[] = {-0.5104644303570166, 0.8598988692516616};//original in paper (u3): {0.5105, -0.8599};
		u3_type u3_check;
		u3_check.set( data_u3, data_u3 + 2);
		
		double data_core[] = { -10.14733447424582, 0.0, 0.0, -2.760705584847321 };
		core_type core_check;
		core_check.set( data_core, data_core + 4);
				
		//test data from lathauwer et al. 2000b paper (same test as in t3_hooi_test
		t3_type t3_data_als;
		double data_als[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6 };
		t3_data_als.set( data_als, data_als + 12 );
		
		tuck3_type tuck3;
		tuck3.tucker_als( t3_data_als, hooi_type::init_hosvd() );
		
		//export
		out_vec_type export_data;
		tuck3_exporter_t::export_to( export_data, tuck3 );
		
		double export_data_check[] = {
			-0.2789473533630371, 0.5983607769012451, 0.751101016998291, 0.4141266345977783, 0.7806354761123657, -0.4680890440940857, 
			0.09816429018974304, 0.9951702356338501, 0.9951702356338501, -0.09816429018974304,
			-0.5104644894599915, 0.8598988652229309,
			-10.14733409881592, -3.874301910400391e-07, 5.960464477539062e-07, -2.76070547103881};
		
		ok = true;
		for (int i = 0; (i < 16) && ok; ++i )
		{
			ok = (abs(export_data_check[i] - export_data[i])) < precision;
		}
		
		log( "export tucker3", ok  );
		
		//import tucker3 from vector
		
		tuck3_type tuck3i;
		tuck3_importer_t::import_from( export_data, tuck3i );
		
		u1_type u1_imported; u2_type u2_imported; u3_type u3_imported;
		core_type core_imported;
		tuck3i.get_core( core_imported ); tuck3i.get_u1( u1_imported  ); tuck3i.get_u2( u2_imported ); tuck3i.get_u3( u3_imported );
		
		if ( u1_imported.equals( u1_check, precision ) && u2_imported.equals( u2_check, precision ) && u3_imported.equals( u3_check, precision ) && core_imported.equals( core_check, precision))
		{	
			log( "import tucker3" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "import tucker3: " << std::endl
			<< "U1 should be: " << std::endl << u1_check << std::endl
			<< "U1 is: " << std::endl << u1_imported << std::endl
			<< "U2 should be: " << std::endl << u2_check << std::endl
			<< "U2 is: " << std::endl << u2_imported << std::endl
			<< "U3 should be: " << std::endl << u3_check << std::endl
			<< "U3 is: " << std::endl << u3_imported << std::endl
			<< "core should be: " << std::endl << core_check << std::endl
			<< "core is: " << std::endl << core_imported << std::endl;
			
			
			log_error( error.str() );
		}		
		
		//tests for quantized exports
		typedef unsigned char T_value;
		typedef unsigned short T_coeff;
		typedef tensor3< 3, 2, 2, T_value> t3q_type; 
		typedef qtucker3_tensor< 2, 2, 2, 3, 2, 2, T_value, T_coeff > tuck3q_type; 
		typedef t3_hooi< 2, 2, 2, 3, 2, 2, float > hooi_type1;
		typedef tucker3_exporter< 2, 2, 2, 3, 2, 2, T_value, T_coeff > tuck3q_exporter_t;
		typedef tucker3_importer< 2, 2, 2, 3, 2, 2, T_value, T_coeff > tuck3q_importer_t;
		typedef std::vector<T_value> out_vecq_type;
		typedef matrix< 3, 2, T_coeff > u1q_type;
		typedef matrix< 2, 2, T_coeff > u2q_type;
		typedef matrix< 2, 2, T_coeff > u3q_type;
		typedef tensor3< 2, 2, 2, T_coeff > coreq_type;
		
		
		unsigned char data_2[] = { 0, 13, 122, 123, 124, 95, 10, 40, 25, 54, 33, 76};
		t3q_type t3_data2;
		t3_data2.set(data_2, data_2 + 12);
		
		
		float core_min, core_max, u_min, u_max;
		
		tuck3q_type tuck3q;
		tuck3q.enable_quantify_linear();
		tuck3q.decompose( t3_data2, u_min, u_max, core_min, core_max, hooi_type1::init_hosvd() );
		
		t3q_type t3_reco2;
		tuck3q.reconstruct( t3_reco2, u_min, u_max, core_min, core_max );
		
		//export bytes (for one min/max value for all Us)
		out_vecq_type export_vec;
		tuck3q_exporter_t::export_quantized_to( export_vec, tuck3q );
		
		//check exported data
		double export_data2_check[] = {
			-0.715831, 0.921996, -251.346, 39.7633,  /* all min/max values -0.715833, 0.865947, -0.721095, 0.692817, -0.921997, 0.921997,*/
			24791, 0, 969, 63293, 12525, 40502,
			56365, 57497, 57497, 921,
			65535, 44136, 13150, 65535,
			0, 56016, 56087, 54573, 56556, 64473, 65535, 52911 };
		
		ok = true;
		
		float * float_ptr = (float*)&(export_vec[0]);
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
			ok = (fabs(float(export_data2_check[index]) - float(value))) < precision;
		}
		//check u2
		end_index += 4;
		for ( ; (index < end_index) && ok ; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### U2 value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < precision;
		}
		//check u3		
		end_index += 4;
		for ( ; (index < end_index) && ok; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### U3 value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < precision;
		}
		
		//check core values
		end_index += 8;
		for ( ; (index < end_index) && ok; ++index ) {
			value = *value_ptr;
			//std::cout<<"#### core value === " << value << ", should be " << export_data2_check[index] << std::endl;
			value_ptr++;
			ok = (fabs(float(export_data2_check[index]) - float(value))) < precision;
		}
		
		log( "export tucker3 (bytes) ", ok  );
		
		
		
		//prepare check data
		T_coeff du1[] = { 24791, 63293, 0, 12525, 969, 40502 };
		u1q_type u1q_check; 
		u1q_check.set(du1, du1 + 6);
		
		T_coeff du2[] = { 56365, 57497, 57497, 921 };
		u2q_type u2q_check;
		u2q_check.set(du2, du2 + 4);
		
		T_coeff du3[] = { 65535, 13150, 44136, 65535 };
		u3q_type u3q_check;
		u3q_check.set(du3, du3 + 4);
		
		T_coeff dco[] = { 0, 56087, 56016, 54573, 56556, 65535, 64473, 52911 };
		coreq_type coreq_check;
		coreq_check.set(dco, dco + 8);
		
		u1q_type u1q; u2q_type u2q; u3q_type u3q;
		coreq_type coreq;
		
		//do import from bytes
		tuck3q_type tuck3qi;
		tuck3q_importer_t::import_quantized_from( export_vec, tuck3qi );
		
		tuck3qi.get_u1( u1q ); tuck3qi.get_u2( u2q ); tuck3qi.get_u3( u3q );
		tuck3qi.get_core( coreq );
		
		ok = u1q == u1q_check;
		ok = ok && (u2q == u2q_check);
		ok = ok && (u3q == u3q_check);
		ok = ok && (coreq == coreq_check);
		
		if (ok)
		{
			log( "import tucker3 (bytes)" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "import tucker3 (bytes): " << std::endl
			<< "Tucker3 is: " << std::endl << tuck3q << std::endl;
			
			log_error( error.str() );
		}	
				
		
		//prepare check data
	
		T_coeff dcho[] = { 65285, 40, 43, 79, 4, 127, 123, 98 };
		coreq_check.zero();
		coreq_check.set(dcho, dcho + 8);
		

		//export bytes (with optimized quantization of core) 
		out_vecq_type export_vec_2;
		tuck3q_exporter_t::export_hot_quantized_to( export_vec_2, tuck3q );
		
		tuck3q_type tuck3qi_2;
		tuck3q_importer_t::import_hot_quantized_from( export_vec_2, tuck3qi_2 );
		
		tuck3qi_2.get_u1( u1q ); tuck3qi_2.get_u2( u2q ); tuck3qi_2.get_u3( u3q );
		tuck3qi_2.get_core( coreq );
		
		ok = u1q == u1q_check;
		ok = ok && (u2q == u2q_check);
		ok = ok && (u3q == u3q_check);
		ok = ok && (coreq == coreq_check);

		if (ok)
		{
			log( "export/import tucker3 (bytes) with hot core quantization" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "export/import tucker3 (bytes) with hot core quantization: " << std::endl
			<< "Tucker3 is: " << std::endl << tuck3qi_2 << std::endl;
			
			log_error( error.str() );
		}	
		
		//prepare check data
		T_coeff dchlo[] = { 127, 26, 28, 52, 0, 85, 82, 65 };
		coreq_check.zero();
		coreq_check.set(dchlo, dchlo + 8);
		
		
		//export bytes for TTM 
		out_vecq_type export_vec_3;
		tuck3q.enable_quantify_log();
		tuck3q.decompose( t3_data2, u_min, u_max, core_min, core_max, hooi_type1::init_hosvd() );
		tuck3q_exporter_t::export_ttm_quantized_to( export_vec_3, tuck3q );
		
		tuck3q_type tuck3qi_3;
		tuck3q_importer_t::import_ttm_quantized_from( export_vec_3, tuck3qi_3 );
		
		tuck3qi_3.get_u1( u1q ); tuck3qi_3.get_u2( u2q ); tuck3qi_3.get_u3( u3q );
		tuck3qi_3.get_core( coreq );
		
		ok = u1q == u1q_check;
		ok = ok && (u2q == u2q_check);
		ok = ok && (u3q == u3q_check);
		ok = ok && (coreq == coreq_check);
		
		if (ok)
		{
			log( "export/import TTM tucker3 (bytes)" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "export/import TTM tucker3 (bytes): " << std::endl
			<< "Tucker3 is: " << std::endl << tuck3qi_3 << std::endl;
			
			log_error( error.str() );
		}	
	
		
		return ok;
	}
	
} //end vmml namespace

