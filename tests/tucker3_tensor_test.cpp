#include "tucker3_tensor_test.hpp"

#include <vmmlib/tucker3_tensor.hpp>

#include <sstream>

namespace vmml
{
	
	bool
	tucker3_tensor_test::run()
	{
		bool ok = false;
		
		//decomposition (hosvd test data after lathauwer 2000a)
		//prepare control data
		matrix<3, 3, double> u1_hosvd;
		matrix<3, 3, double> u2_hosvd;
		matrix<3, 3, double> u3_hosvd;
		matrix<3, 3, double> u1_hosvd_check;
		matrix<3, 3, double> u2_hosvd_check;
		matrix<3, 3, double> u3_hosvd_check;
		double data_u1_hosvd[] = { 0.1121, -0.7739, -0.6233, 0.5771, 0.5613, -0.5932, 0.8090, -0.2932, 0.5095 };
		u1_hosvd_check.set( data_u1_hosvd, data_u1_hosvd + 9);
		double data_u2_hosvd[] = { 0.4624, 0.0102, 0.8866, 0.8866, -0.0135, -0.4623, -0.0072, -0.9999, -0.0152 };
		u2_hosvd_check.set( data_u2_hosvd, data_u2_hosvd + 9);
		double data_u3_hosvd[] = { 0.6208, -0.4986, 0.6050, -0.0575, -0.7986, -0.5992, 0.7819, 0.3371, -0.5244 };
		u3_hosvd_check.set( data_u3_hosvd, data_u3_hosvd + 9);
		
		tensor3<3, 3, 3, double> core_hosvd;
		tensor3<3, 3, 3, double> core_hosvd_check;
		double data_core_hosvd[] = { 
			8.7088, 0.0489, -0.2797, -0.0256, 3.2546, -0.2853, 0.0, 0.0, 0.0, 
			0.1066, 3.2737, 0.3223, 3.1965, -0.2130, 0.7829, 0.0, 0.0, 0.0,
			-0.0033, -0.1797, -0.2222, 0.2948, -0.0378, -0.3704, 0.0, 0.0, 0.0
		};
		core_hosvd_check.set( data_core_hosvd, data_core_hosvd + 27);

		
		//do decomposition
		tensor3< 3, 3, 3, double> t3_data_hosvd;
		double data_hosvd[] = { 
			0.9073, 0.7158, -0.3698, 0.8924, -0.4989, 2.4288, 2.1488, 0.3054, 2.3753, 
			1.7842, 1.6970, 0.0151, 1.7753, -1.5077, 4.0337, 4.2495, 0.3207, 4.7146,
			2.1236, -0.704, 1.4429, -0.6631, 1.9103, -1.7495, 1.8260, 2.1335, -0.2716
		};
		t3_data_hosvd.set(data_hosvd, data_hosvd + 27);
		

		tucker3_tensor< 3, 3, 3, 3, 3, 3, double > tuck3_hosvd( t3_data_hosvd, u1_hosvd, u2_hosvd, u3_hosvd );
		//tuck3_hosvd.decomposition( t3_data_hosvd );		
		
		//std::cout << u1_hosvd_check << std::endl << u2_hosvd_check << std::endl << u3_hosvd_check << std::endl;
		
		//derive core tensor
		/*tuck3_hosvd.derive_core( t3_data_hosvd, core_hosvd, u1_hosvd_check, u2_hosvd_check, u3_hosvd_check );
		
		
		if ( core_hosvd.equals( core_hosvd_check, 0.1 ))
		{	
			log( "HOSVD derive core tensor", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOSVD derive core: " << std::endl
			<< "core should be: " << core_hosvd_check << std::endl
			<< "core is: " << core_hosvd << std::endl;
			
			log_error( error.str() );
		}*/
		
		
		//decomposition
		//tucker method1
		//prepare control data
		matrix<6, 6, double> u1_deco;
		matrix<7, 7, double> u2_deco;
		matrix<5, 5, double> u3_deco;
		tensor3<6, 7, 5, double> core_deco;
		double data_u1_deco[] = { 0.39409, 0.31438, 0.79477, -0.038406, -0.31008, 0.12873, 0.4169, -0.53849, -0.11145, -0.011199, -0.51248, -0.51092,
			0.4411, 0.62016, -0.48386, -0.40984, -0.052559, -0.12644, 0.35799, -0.38185, -0.19023, -0.20749, -0.061116, 0.80191,
		0.40864, -0.24275, 0.22764, -0.2222, 0.78387, -0.24176, 0.42565, 0.14788, -0.18382, 0.85907, 0.14237, 0.069925 };
		u1_deco.set( data_u1_deco, data_u1_deco + 36);
		double data_u2_deco[] = { 0.37977, -0.27662, -0.37921, -0.43239, -0.60914, -0.27689, 0.027818, 0.41347, 0.59511, -0.47933, -0.13254, 0.32123, 0.095224, -0.33956,
			0.36785, 0.39565, 0.051661, 0.57533, -0.32655, -0.0903, 0.50959, 0.32284, -0.52527, -0.16687, 0.61675, 0.10942, -0.052308, -0.44378,
			0.39252, -0.11557, 0.30416, -0.14084, -0.17109, 0.83077, -0.027864, 0.36717, 0.14822, 0.70948, -0.15908, 0.020459, -0.44672, -0.33849,
		0.39545, -0.31869, -0.0071048, -0.19736, 0.61437, -0.11677, 0.55853 };
		u2_deco.set( data_u2_deco, data_u2_deco + 49);
		double data_u3_deco[] = { 0.4429, 0.20937, -0.34207, 0.61521, 0.5143, 0.37964, 0.69004, -0.30785, -0.47388, -0.24573, 0.51481, -0.30504, 0.33105, 
		-0.51488, 0.51694, 0.43335, 0.1715, 0.69319, 0.36143, -0.4143, 0.45486, -0.59795, -0.4451, 0.034881, -0.48604 };
		u3_deco.set( data_u3_deco, data_u3_deco + 25);
		
		double data_core_deco[] = { 2007.7, -2.3252, -9.037, -8.3985, 15.156, -3.2019, -1.5649, 3.2852, 8.9891, 181.9, 166.99, -20.471, 18.155, 58.109, 13.389, -10.986, -40.448, 9.5224, -113.94, 1.6497, 42.873, -6.9224, -40.183, -34.849, 33.371, -59.766, -142.27, 8.0361, 
			-11.799, 74.117, 60.689, 17.149, -24.997, 49.599, 7.1353, -6.4525, -72.047, 1.687, 106.31, 4.7981, -71.856, -1.8353, 16.323, 22.606, 29.91, -17.912, -63.184, -40.601, -53.833, 81.257, -56.533, -67.03, 32.031, 44.125, 157.82, 81.359, 
			-131.17, -210.35, -20.246, -60.333, -87.536, 12.731, 43.177, -85.653, 88.362, 83.687, -105.87, 35.914, 132.23, 12.536, 21.216, -65.823, 124.14, -134.18, 67.559, -68.979, 10.943, 113.23, -2.5747, 56.162, -59.419, -136.28, 22.089, -52.339, 
			-2.3786, 23.036, -131.44, -9.4512, 23.962, 23.488, -64.454, -53.058, 158.79, 5.0995, -86.745, 47.356, -140.88, 0.81563, 28.832, -3.4195, 98.289, -110.88, -74.096, 0.13307, 126.49, -117.83, 26.728, -82.889, -68.421, -26.545, 14.857, -49.305, 
			-78.954, 26.993, 58.295, 0.72088, 47.81, 36.51, -87.717, 97.861, -35.38, 37.791, -1.4775, 34.25, -25.063, -8.1026, 2.1341, 190.25, -80.073, 45.774, 4.0446, 17.17, 70.644, -124.02, -39.105, -40.683, -35.555, 128.88, -2.9354, 56.124, 
			-50.355, 33.157, 58.033, 117.73, -47.401, 2.9577, -71.132, 1.2545, 40.366, 35.36, 98.407, -9.2858, 4.1012, -2.6178, 48.598, 101.1, 34.827, -55.407, -39.182, 26.673, -51.725, 20.17, -37.401, -48.994, 2.3187, -68.968, 47.912, -8.6009, 
			7.6446, 73.173, -51.285, 1.8827, -73.729, 55.471, -10.906, 40.436, 125.49, 30.189, -7.0578, -106.45, -33.484, 87.351, -76.129, 55.187, -33.385, -11.311, 38.171, -57.902, 56.593, 
			12.727, -42.201, 134, -32.284, -13.556, 23.213, -59.794, 7.9093, -15.15, -90.559, -24.622, -118.46, -27.751, 1.521, -30.563, 14.512, -8.2528, -73.867, -42.656, 36.996, -33.781		
		};
		
		core_deco.set( data_core_deco, data_core_deco + 210);

		//do decomposition
		tensor3<6, 7, 5, double> t3_2;
		double data_decomposition[] = { 63.0, 97, 202, 133,  84,  31,  51, 99, 236, 115, 110,  87, 159, 171, 110, 189, 217,  53, 208,  88, 230, 212, 188, 100,  82, 136,  85,  51,
		210, 241, 188,  28, 133, 147,  76, 116, 130, 249,  96, 197, 220, 127, 227,   4,  60, 198, 196,  36, 104, 128, 179,  52,  42,  18, 131,  61,
		71, 129, 208, 233, 242, 184, 133, 136,  97, 100,  81,  40, 237,  56, 146,  17,  14,  84,  73, 187, 215, 105,  91,  96,  52, 175, 191, 169,
		208, 196, 221,  27, 147, 211, 212, 202, 212,  17, 185, 241,  82, 235, 120, 180, 247, 156, 222, 249,  83, 79, 152,  25, 200, 129,  71, 205,
		175, 192, 139, 144, 201,  19, 137, 252, 127, 103, 207, 121, 192, 118, 209,  68, 175,  94, 146, 112, 150, 243, 212, 145, 119, 124,  21, 117,
		19, 254,  97, 128,  67, 144, 220, 181, 166, 162, 232, 148, 138, 169, 60, 180,  93,  53, 224, 196,  90, 102, 238, 104,  86,  16,  59,  89,
		65, 180, 152, 220, 140, 190, 135, 243, 179,  76,  21, 191, 186, 219, 76,   2,  32,  86, 215, 183, 173, 40,  95,  99,  60,  43,  34, 205,
		92, 230, 209,  81, 230, 114, 135, 189,  81, 250, 251,  27, 130, 244};
		
		tensor3<6, 7, 5, double> core_tmp;
		t3_2.set(data_decomposition, data_decomposition + 210);
		
		//part 1: derive core
		matrix<6, 6, double> u1_2;
		matrix<7, 7, double> u2_2;
		matrix<5, 5, double> u3_2;
		tucker3_tensor<6, 7, 5, 6, 7, 5, double > tuck3_deco( t3_2, u1_2, u2_2, u3_2 );
		tuck3_deco.derive_core(t3_2, core_tmp, u1_deco, u2_deco, u3_deco);
		if ( core_tmp.equals( core_deco, 0.1 ))
		{	
			log( "tucker3 derive core tensor", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker3 derive core: " << std::endl
			<< "core should be: " << core_deco << std::endl
			<< "core is: " << core_tmp << std::endl;
			
			log_error( error.str() );
		}
		
		//part 2: compute basis matrices
		//tuck3_deco.decomposition( t3_2 );		
		//TODO
		ok = true;
		matrix<6, 6, double> u1_tmp;
		u1_tmp = tuck3_deco.get_u1();
		matrix<7, 7, double> u2_tmp;
		u2_tmp = tuck3_deco.get_u2();
		matrix<5, 5, double> u3_tmp;
		u3_tmp = tuck3_deco.get_u3();
		
		//if ( u1_tmp == u1_deco && u2_tmp == u2_deco && u3_tmp == u3_deco )
		if ( ok)
		{	
			log( "tucker3 compute basis matrices", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker3 compute basis matrices: " << std::endl
			<< "u1 should be: " << u1_deco << std::endl
			<< "u1 is: " << u1_tmp << std::endl
			<< "u2 should be: " <<  u2_deco << std::endl
			<< "u2 is: " <<  u2_tmp << std::endl
			<< "u3 should be: " << u3_deco << std::endl
			<< "u3 is: " << u3_tmp << std::endl;
			
			log_error( error.str() );
		}
		
		
		
		//tucker3 reconstruction
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
		

		//rank reduction
		tensor3< 1, 2, 3, uint16_t >  core_red;
		matrix<6, 1, uint16_t> u1_red;
		matrix<7, 2, uint16_t> u2_red;
		matrix<5, 3, uint16_t> u3_red;
		
		tucker3_tensor< 1, 2, 3, 6, 7, 5, uint16_t > tuck3_red( core_red, u1_red, u2_red, u3_red );
		tuck3_red.progressive_rank_reduction( tuck3 );
		
		u1_red.fill(2);
		u2_red.fill(3);
		u3_red.fill(1);
		double data[] = { 0, 1, 6, 7, 12, 13 };
		core_red.set(data, data+6);
		
		if ( tuck3_red.get_u1() == u1_red && tuck3_red.get_u2() == u2_red && tuck3_red.get_u3() == u3_red && tuck3_red.get_core() == core_red)
		{	
			log( "tucker3 progressive rank reduction", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "Tucker3 progressive rank reduction: " << std::endl
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
		
		//basis matrices subsampling, average data
		tensor3< 3, 4, 3, uint16_t > t3_sub_avg;
		tensor3< 2, 3, 4, uint16_t > core_sub_avg;
		matrix< 3, 2, uint16_t > u1_sub_avg;
		matrix< 4, 3, uint16_t > u2_sub_avg;
		matrix< 3, 4, uint16_t > u3_sub_avg;
		tucker3_tensor< 2, 3, 4, 3, 4, 3, uint16_t > tuck3_sub_avg( core_sub_avg, u1_sub_avg, u2_sub_avg, u3_sub_avg );
		
		/*matrix<6, 2, uint16_t> u1_new;
		u1_new.fill(5);
		u1_new.at(1,1) = 2;
		u1_new.at(4,0) = 3;
		u1_new.at(3,1) = 8;
		tuck3.set_u1(u1_new);*/
		
		tuck3_sub_avg.subsampling_on_average( tuck3, 2);
		tuck3_sub_avg.reconstruction( t3_sub_avg );
		
		t3_sub_test.fill(1656);
		if ( t3_sub_test == t3_sub_avg )
		{	
			log( "basis matrices subsampling on average", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "basis matrices subsampling on average with factor 2: " << std::endl << t3_sub_avg
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

