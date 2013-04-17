#define CP_LOG 1
#define TUCKER_LOG 1

#include <iostream>

#include "include/vmmlib/tensor3.hpp"
#include "include/vmmlib/tensor4.hpp"
#include "include/vmmlib/t3_converter.hpp"
#include "include/vmmlib/tensor_mmapper.hpp"
#include "include/vmmlib/t3_ttm.hpp"
#include "include/vmmlib/tucker3_tensor.hpp"
#include "include/vmmlib/cp3_tensor.hpp"

using namespace vmml;


int main (int argc, char * const argv[]) {

    std::cout.setf(std::ios::fixed);
    std::cout.precision(16);
    std::cerr.precision(16);
	//**** create a tensor of size axbxc or 4x2x3 and with unsigned char values
	const size_t a = 4;
	const size_t b = 3;
	const size_t c = 2;
	typedef unsigned char values_t;
 	typedef tensor3< a, b, c, values_t > t3_t;
	t3_t t3; //after initializing a tensor3, the tensor is still empty
	t3.fill_increasing_values(); //fills the empty tensor with the values 0,1,2,3...
	
	//the operators << are overwritten such that you can directly output to the console
//	std::cout << "test tensor3:\n" << t3 << std::endl; 
	
 	typedef tensor4< a, b, c, c, values_t > t4_t;
	t4_t t4; //after initializing a tensor4, the tensor is still empty
	t4.fill_increasing_values(); //fills the empty tensor with the values 0,1,2,3...
	
	//the operators << are overwritten such that you can directly output to the console
//	std::cout << "test tensor4:\n" << t4 << std::endl; 
	
	
	//**** load a tensor3 from a file and memory map it
	const size_t d = 512;
	typedef tensor3< d,d,d, values_t > t3_512u_t;
	typedef t3_converter< d,d,d, values_t > t3_conv_t;
	typedef tensor_mmapper< t3_512u_t, t3_conv_t > t3map_t;
	
	std::string in_dir = "/home/rballester/datasets/survey_testdata";
	std::string file_name = "hnut512_uint.raw";
	t3_512u_t t3_hazelnut;
	t3_conv_t t3_conv;

	t3map_t t3_mmap( in_dir, file_name, true, t3_conv ); //true means that the file is read-only
	t3_mmap.get_tensor( t3_hazelnut );
	
	//**** get slices: this can be done as a forward or a backward procedure (Kiers, 2000 or Lathauwer et al., 2000). both types are implemented
	//the tensor3 is organized as I3 frontal matrices (slices)
	
	typedef matrix< d, d, values_t > matrix_t;
	//get frontal slice
	matrix_t slice;
	t3_hazelnut.get_frontal_slice_fwd( d/2, slice );
	slice.write_to_raw( in_dir, "frontal_slice.raw");

	//get horizontal slice
	t3_hazelnut.get_horizontal_slice_fwd( d/2, slice );
	slice.write_to_raw( in_dir, "horizontal_slice.raw");

	//get lateral slice
	t3_hazelnut.get_lateral_slice_fwd( d/2, slice );
	slice.write_to_raw( in_dir, "lateral_slice.raw");
	
	//**** unfoldings
	
	//frontal unfoldings after Kiers, 2000
	typedef matrix< a, b*c, values_t > fwd_front_unfolding_t;
	typedef matrix< b, c*a, values_t > fwd_horiz_unfolding_t;
	typedef matrix< c, a*b, values_t > fwd_lat_unfolding_t;
	
	fwd_front_unfolding_t unf_front_fwd;
	fwd_horiz_unfolding_t unf_horiz_fwd;
	fwd_lat_unfolding_t unf_lat_fwd;

	t3.frontal_unfolding_fwd( unf_front_fwd );     //I1x(I2*I3)
	t3.horizontal_unfolding_fwd( unf_horiz_fwd );  //I2x(I3*I1)
	t3.lateral_unfolding_fwd( unf_lat_fwd );       //I3x(I1*I2)
	
//	std::cout << "forward unfolded tensor (frontal)\n" << unf_front_fwd << std::endl;
//	std::cout << "forward unfolded tensor (horizontal)\n" << unf_horiz_fwd << std::endl;
//	std::cout << "forward unfolded tensor (lateral)\n" << unf_lat_fwd << std::endl;
	

	//backward unfoldings after Lathauwer et al., 2000
	typedef matrix< a, b*c, values_t > bwd_lat_unfolding_t;
	typedef matrix< b, a*c, values_t > bwd_front_unfolding_t;
	typedef matrix< c, a*b, values_t > bwd_horiz_unfolding_t;
	
	bwd_front_unfolding_t unf_front_bwd;
	bwd_horiz_unfolding_t unf_horiz_bwd;
	bwd_lat_unfolding_t unf_lat_bwd;

	t3.lateral_unfolding_bwd( unf_lat_bwd );     //I1x(I2*I3)
	t3.frontal_unfolding_bwd( unf_front_bwd );     //I2x(I1*I3)
	t3.horizontal_unfolding_bwd( unf_horiz_bwd );  //I3x(I1*I2)
		
//	std::cout << "backward unfolded tensor (lateral)\n" << unf_lat_bwd << std::endl;
//	std::cout << "backward unfolded tensor (frontal)\n" << unf_front_bwd << std::endl;
//	std::cout << "backward unfolded tensor (horizontal)\n" << unf_horiz_bwd << std::endl;
	
	//**** tensor times matrix multiplications (example with fwd multiplications, bwd is also possible)
	const size_t f = 2;
	matrix< f, a, values_t > m1;
	matrix< f, b, values_t > m2;
	matrix< f, c, values_t > m3;
	
	m1.fill( 1 );
	m2.fill( 2 );
	m3.fill( 3 );
	
	//resulting tensors after ttm multiplications
 	tensor3< f, b, c, values_t > t3_out_1;
 	tensor3< a, f, c, values_t > t3_out_2;
 	tensor3< a, b, f, values_t > t3_out_3;
 	tensor3< f, f, f, values_t > t3_out;
	
	t3_ttm::multiply_frontal_fwd(    t3, m1, t3_out_1 ); 
	t3_ttm::multiply_horizontal_fwd( t3, m2, t3_out_2 ); 
	t3_ttm::multiply_lateral_fwd(    t3, m3, t3_out_3 ); 
	
	//does a forward cyclic TTM along mode 1, 2, 3
	t3_ttm::full_tensor3_matrix_multiplication( t3, m1, m2, m3, t3_out );
	
	std::cout << "t3_out after TTM1, TTM2, and TTM3:\n" << t3_out << std::endl;
	
	
	//**** Tucker3 tensor 
	
	const size_t r = 2;
	typedef tucker3_tensor< r, r, r, a, b, c, values_t, float > tucker3_t;
	//Various constructors available: For example only core, core and the three factor matrices, the input data and the factor matrices or another tucker3 tensor
	tucker3_t tuck3_dec; //empty tucker3 tensor
	typedef t3_hooi< r, r, r, a, b, c, float > hooi_t;
	

	//Decomposition or Tucker ALS
	
	//choose initialization of Tucker ALS (init_hosvd, init_random, init_dct)

	//Example for initialization with init_rand
	tuck3_dec.tucker_als( t3, hooi_t::init_random());
//	std::cout << "Tucker3 decomposition (init_random):\n" << tuck3_dec << std::endl;
//
//	//Example for initialization with init_hosvd
//	tuck3_dec.tucker_als( t3, hooi_t::init_hosvd());
//	std::cout << "Tucker3 decomposition (init_hosvd):\n" << tuck3_dec << std::endl;
//	
//	//Reconstruction
//	t3_t t3_reco;
//	tuck3_dec.reconstruct( t3_reco );
//	std::cout << "reconstructed tensor:\n" << t3_reco << std::endl;
//	double rms_err = t3.rmse( t3_reco ) ;
//	std::cout << "RMSE original and Tucker3 tensor approximation:\t" << rms_err << std::endl;
//
//	//NOTE: set in t3_hooi.cpp the variable '#define TUCKER_LOG 1' if you like to see the number of iterations and the fit as console output
//
//
////	
////	//**** CP3 tensor
////	typedef cp3_tensor< r, a, b, c, values_t, float > cp3_t;
////	typedef t3_hopm< r, a, b, c, float > t3_hopm_t;
////	
////	cp3_t cp3_dec;
////	
////	//Decomposition or CP ALS
////	//choose initialization of Tucker ALS (init_hosvd, init_random)
////
////	int max_cp_iter = 20;
////	cp3_dec.cp_als( t3, t3_hopm_t::init_random(), max_cp_iter );
////	std::cout << "CP3 decomposition (init_random):\n" << cp3_dec << std::endl;
////	
////	t3_t t3_cp_reco;
////	cp3_dec.reconstruct( t3_cp_reco );
////	std::cout << "reconstructed tensor:\n" << t3_reco << std::endl;
////	rms_err = t3.rmse( t3_cp_reco ) ;
////	std::cout << "RMSE original and CP3 tensor approximation:\t" << rms_err << std::endl;
////	
////	
////	//NOTE: set in t3_hopm.cpp the variable '#define CP_LOG 1' if you like to see the number of iterations and the fit as console output
////
////	
////	
//    std::cout << "\nEnded vmmlib tensor classes demo!\n";
//    return 0;
}
