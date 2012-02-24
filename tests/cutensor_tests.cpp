#include "cutensor_tests.hpp"

#ifdef VMMLIB_CUDA

#include <vmmlib/t3_cublas_ttm.hpp>
#include <vmmlib/t3_cublas_hosvd.hpp>
#include <vmmlib/t3_cublas_hooi.hpp>

namespace vmml
{
	
	bool
	cutensor_tests::run()
	{
		bool ok = false;
		
		tensor3< 2, 3, 4, int >  t3;
		
		//CUDA TTM for tensor3
		t3.fill_increasing_values();
		matrix<5, 4, int> u3;
		u3.fill(1);
		tensor3<2, 3, 5, int> t3_jji; t3_jji.zero();
		t3_cublas_ttm::multiply_horizontal_bwd(t3, u3, t3_jji );
				
		tensor3<6, 3, 5, int> t3_iji; t3_iji.zero();
		matrix<6, 2, int> u1;
		u1.fill(2);
		t3_cublas_ttm::multiply_lateral_bwd(t3_jji, u1, t3_iji );
		
		tensor3<6, 7, 5, int> t3_iii; t3_iii.zero();
		matrix<7, 3, int> u2;
		u2.fill(3);
		t3_cublas_ttm::multiply_frontal_bwd(t3_iji, u2, t3_iii );		
		
		t3.fill_increasing_values();
		tensor3<6, 7, 5, int> t3_reco;
		t3_cublas_ttm::full_tensor3_matrix_multiplication( t3, u1, u2, u3, t3_reco ); 
		tensor3<6, 7, 5, int> t3_reco2;
		t3_cublas_ttm::full_tensor3_matrix_kronecker_mult( t3, u1, u2, u3, t3_reco2 ); 

		ok = t3_reco == t3_reco2;
		
		tensor3<6, 7, 5, int> t3_iii_test;
		t3_iii_test.fill(1656);
		if ( t3_iii_test == t3_reco && t3_iii_test == t3_iii && ok )
		{	
			log( "tensor3 CUBLAS matrix multiplication along all three modes", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "T3_result (all values should be 1656): " << std::endl << t3_iii
			<< std::endl;
			log_error( error.str() );
		}
		
		

		//CUDA HOEIGS for tensor3
		double precision = 0.001;
		
		//hoeigs test
		// hosvd test data after lathauwer et al. 2000a
		tensor3< 3, 3, 3, double > t3_data_hoeigs;
		double data_lathw[] = { 
			0.9073, 1.7842, 2.1236, 0.8924, 1.7753, -0.6631, 2.1488, 4.2495, 1.8260, 
			0.7158, 1.6970, -0.0704, -0.4898, -1.5077, 1.9103, 0.3054, 0.3207, 2.1335, 
			-0.3698, 0.0151, 1.4429, 2.4288, 4.0337, -1.7495, 2.3753, 4.7146, -0.2716
		};
		t3_data_hoeigs.set(data_lathw, data_lathw + 27);

		matrix< 3, 3, double > u1_hoeigs;
		matrix< 3, 3, double > u2_hoeigs;
		matrix< 3, 3, double > u3_hoeigs;
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
		
		typedef t3_cublas_hosvd< 3, 3, 3, 3, 3, 3, double > hosvd_t;
		hosvd_t::apply_all( t3_data_hoeigs, u1_hoeigs, u2_hoeigs, u3_hoeigs );
		
		if ( u1_hoeigs.equals( u1_hoeigs_check, precision ) && u2_hoeigs.equals( u2_hoeigs_check, precision ) && u3_hoeigs.equals( u3_hoeigs_check, precision ))
		{	
			log( "CUBLAS HOEIGS compute factor matrices U1, U2, U3", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOEIGS: " << std::endl
			<< "U1 should be: " << std::endl << u1_hoeigs_check << std::endl
			<< "U1 is: " << std::endl << u1_hoeigs << std::endl
			<< "U2 should be: " << std::endl << u2_hoeigs_check << std::endl
			<< "U2 is: " << std::endl << u2_hoeigs << std::endl
			<< "U3 should be: " << std::endl << u3_hoeigs_check << std::endl
			<< "U3 is: " << std::endl << u3_hoeigs << std::endl;
			
			log_error( error.str() );
		}
		
		//CUDA HOOI for tensor3
		
		//test data from lathauwer et al. 2000b paper
		tensor3< 3, 2, 2, double> t3_data;
		double data[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6 };
		t3_data.set( data, data + 12 );
		
		matrix< 3, 2, double > u1_hooi; u1_hooi.zero();
		matrix< 2, 2, double > u2_hooi; u2_hooi.zero();
		matrix< 2, 1, double > u3_hooi; u3_hooi.zero();
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
		
		tensor3< 2, 2, 1, double > core_hooi;
		tensor3< 2, 2, 1, double > core_check;
		double data_core[] = { -10.14733447424582, 0.0, 0.0, -2.760705584847321 };
		core_check.set( data_core, data_core + 4);
		
		typedef t3_cublas_hooi< 2, 2, 1, 3, 2, 2, double > hooi_cublas_t;
		hooi_cublas_t::als( t3_data, u1_hooi, u2_hooi, u3_hooi, core_hooi, hooi_cublas_t::init_hosvd() );
		
		ok = u1_hooi.equals( u1_check, precision );
		ok = ok && ( u2_hooi.equals( u2_check, precision ) );
		ok = ok && ( u3_hooi.equals( u3_check, precision ) );
		ok = ok && ( u3_hooi.equals( u3_check, precision ) );
		ok = ok && ( core_hooi.equals( core_check, precision ) );
		
		if ( ok )
		{	
			log( "HOOI rank-(2,2,1) approximation" , ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOOI rank-(2,2,1) approximation: " << std::setprecision(16) << std::endl
			<< "U1 should be: " << std::endl << u1_check << std::endl
			<< "U1 is: " << std::endl << u1_hooi << std::endl
			<< "U2 should be: " << std::endl << u2_check << std::endl
			<< "U2 is: " << std::endl << u2_hooi << std::endl
			<< "U3 should be: " << std::endl << u3_check << std::endl
			<< "U3 is: " << std::endl << u3_hooi << std::endl
			<< "core should be: " << std::endl << core_check << std::endl
			<< "core is: " << std::endl << core_hooi << std::endl;
			
			
			log_error( error.str() );
		}
				
		return ok;
	}
	
	
	
} // namespace vmml

#endif
