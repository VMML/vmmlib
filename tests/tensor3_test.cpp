#include "tensor3_test.hpp"

#include <vmmlib/t3_converter.hpp>
#include <vmmlib/tensor3.hpp>
#include <sstream>

namespace vmml
{

bool
tensor3_test::run()
{
    bool ok = false;
        
	
	tensor3< 2, 3, 4, int >  t3;
    tensor3< 2, 3, 4, int >  t3_tmp;

	//test size
	if (t3.size() == 24)
	{	
		log( "size()", true  );
	} else
	{
		std::stringstream error;
		error << "size should be 24, but size is: " << t3.size() << std::endl;
		log_error( error.str() );
	}

	
	//test at()
	//TODO
    t3.at( 0,0,0 ) = 255.0;
    t3( 0,2,0 ) = 128.0;
	ok = true;
	if (ok)
	{	
		log( "at() ", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3, 255 @ (0, 0, 0) and 128 @ (0, 2, 0): " << std::endl << t3 << std::endl;
		log_error( error.str() );
	}
	
	
	//test set tensor from input
	int data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };
	t3.set(data, data + 24);
	ok = false;
	int check_value = 1;
	for ( size_t i3 = 0; i3 < 4; ++i3 )
	{
		for( size_t i1 = 0; i1 < 2; ++i1 )
		{
			for( size_t i2 = 0; i2 < 3 && ok; ++i2, ++check_value )
			{
				if( t3.at(i1, i2, i3) != check_value)
				{
					std::stringstream error;
					error << "T3 set from input, values from 1... 24: " << std::endl << t3 << std::endl;
					log_error( error.str() );
					ok = false;
				} else {
					ok = true;
				}
				
			}
		}
	}
	if (ok)
	{	
		log( "set() ", true  );
	}
	
	
	//test fil()
 	t3.fill( 4.0 );
	ok = false;
	check_value = 4;
	for ( size_t i3 = 0; i3 < 4; ++i3 )
	{
		for( size_t i1 = 0; i1 < 2; ++i1 )
		{
			for( size_t i2 = 0; i2 < 3 && ok; ++i2 )
			{
				if( t3.at(i1, i2, i3) != check_value )
				{
					std::stringstream error;
					error << "T3 with all values to 4: " << std::endl << t3 << std::endl;
					log_error( error.str() );
					ok = false;
				} else {
					ok = true;
				}
				
			}
		}
	}
	if(ok) 
	{
		log( "fill()", true  );
	}
	
	
	//operations
    tensor3< 2, 3, 4, int >  t3_2;
    tensor3< 2, 3, 4, int >  t3_3;
    tensor3< 2, 3, 4, int >  t3_result;
	t3.fill( 2 );
	t3_2.fill( 5 );
	t3_3.fill( 5 );
	
	//test equals operator
	if ( !(t3==t3_2) && (t3_2==t3_3) && t3 != t3_2 && !(t3_2!=t3_3))
	{	
		log( "operator== and operator!= ", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3: " << std::endl << t3
		<< "T3_2: " << std::endl << t3_2
		<< "T3_3: " << std::endl << t3_3
		<< std::endl;
		log_error( error.str() );
	}
	//test sum with other t3
	t3_result = t3 + t3_2;
	t3_3.fill( 7 );
	if ( t3_result == t3_3)
	{	
		log( "operator+ and operator += with other tensor3", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3_result: " << std::endl << t3_result
		<< "T3_3: " << std::endl << t3_3
		<< std::endl;
		log_error( error.str() );
	}
	
	//test subtraction
	t3_result = t3_2 - t3;
	t3_3.fill( 3 );
	if ( t3_result == t3_3)
	{	
		log( "operator- and operator -= with other tensor3 ", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3_result: " << std::endl << t3_result
		<< "T3_3: " << std::endl << t3_3
		<< std::endl;
		log_error( error.str() );
	}
	
	
	//test sum with scalar (shift)
	t3_result = t3 + 4;
	t3_3.fill( 6 );
	if ( t3_result == t3_3)
	{	
		log( "operator+ and operator += with scalar", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3_result: " << std::endl << t3_result
		<< "T3_3: " << std::endl << t3_3
		<< std::endl;
		log_error( error.str() );
	}
	
	//test subtraction with scalar (negative shift)
	t3_result = t3 - 2;
	t3_3.zero();
	if ( t3_result == t3_3 )
	{	
		log( "operator- and operator -= with scalar", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3_result: " << std::endl << t3_result
		<< "T3_3: " << std::endl << t3_3
		<< std::endl;
		log_error( error.str() );
	}
	
	//test spherical weighting
	tensor3< 3,3,3, double > t3_result_sw_check;
	double data_sw[] = { 0.176921194126, 0.243116716007, 0.176921194126, 0.243116716007, 0.367879404384, 0.243116716007, 0.176921194126, 0.243116716007, 0.176921194126,
		0.243116716007, 0.367879404384, 0.243116716007, 0.367879404384, 0.9999999, 0.367879404384, 0.243116716007, 0.367879404384, 0.243116716007,
		0.176921194126, 0.243116716007, 0.176921194126, 0.243116716007, 0.367879404384, 0.243116716007, 0.176921194126, 0.243116716007, 0.176921194126};
	t3_result_sw_check.set(data_sw, data_sw + 27);
	
	tensor3< 3,3,3, unsigned int > t3_sw;
	tensor3< 3,3,3, double > t3_result_sw;
	t3_sw.fill( 1 );
	t3_sw.apply_spherical_weights( t3_result_sw );
	double precision = 1.0e-5;
	if ( t3_result_sw.equals( t3_result_sw_check, precision ) )
	{	
		log( "apply spherical weights", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "Apply spherical weights " << std::setprecision(12)
		<< "is: " << std::endl << t3_result_sw
		<< "should be: " << std::endl << t3_result_sw_check
		<< std::endl;
		log_error( error.str() );
	}
	

	
	
	//test fill increasing values
	t3.fill_increasing_values();
	
	ok = false;
	check_value = 0;
	for ( size_t i3 = 0; i3 < 4; ++i3 )
	{
		for( size_t i1 = 0; i1 < 2; ++i1 )
		{
			for( size_t i2 = 0; i2 < 3 && ok; ++i2, ++check_value )
			{
				if( t3.at(i1, i2, i3) != check_value)
				{
					std::stringstream error;
					error << "T3 with all values from 0 to 23:: " << std::endl << t3 << std::endl;
					log_error( error.str() );
					ok = false;
				} else {
					ok = true;
				}
				
			}
		}
	}
	if (ok)
	{	
		log( "fill_increasing_values()", true  );
	} 
	
	
	//test zero()
	t3.zero();
	ok = false;
	check_value = 0;
	for ( size_t i3 = 0; i3 < 4; ++i3 )
	{
		for( size_t i1 = 0; i1 < 2; ++i1 )
		{
			for( size_t i2 = 0; i2 < 3 && ok; ++i2 )
			{
				if( t3.at(i1, i2, i3) != check_value )
				{
					std::stringstream error;
					error << "T3 with all values to 0: " << std::endl << t3 << std::endl;
					log_error( error.str() );
					ok = false;
				} else {
					ok = true;
				}
				
			}
		}
	}
	if(ok) 
	{
		log( "zero()", true  );
	}
	
	//test fill_random and fill_random_signed -> not checked, since every time new values -> cannot be tested with values
	t3.fill_random();

    tensor3< 2, 3, 4, int16_t >  t3s;
	t3s.fill_random_signed();

	
	//test get_n_vector functions
	int i1 = 1;
	int i2 = 2;
	int i3 = 3;
	t3.fill_increasing_values();
	t3 += 1;

	vmml::vector< 3, int > test_I2_data ; 
	vmml::vector< 3, int > I2_data ; 
	t3.get_row( i1, i3, I2_data );
	test_I2_data.set(22, 23, 24);
	
	vmml::vector< 2, int > test_I1_data ; 
	vmml::vector< 2, int > I1_data ; 
	t3.get_column( i2, i3, I1_data );
	test_I1_data.set(21, 24);
	
	vmml::vector< 4, int > test_I3_data ; 
	vmml::vector< 4, int > I3_data ; 
	t3.get_I3_vector( i1, i2, I3_data );
	test_I3_data.set( 6, 12, 18, 24);
	
	
	if (I2_data == test_I2_data && I1_data == test_I1_data && I3_data == test_I3_data)
	{	
		log( "get_n_vector/get_row/get_column/get_tube()", true  );
	} else
	{
		std::stringstream error;
		error 
		     << "I2_vector (22, 23, 24): " << I2_data << std::endl
		     << "I1_vector (21, 24): " << I1_data << std::endl
			 << "I3_vector (6, 12, 18, 24): " << I3_data << std::endl;
		log_error( error.str() );
	}

	I2_data.set( 1 );
	test_I2_data.set( 1, 1, 1 );
	t3.set_row( i1, i3, I2_data );
	vmml::vector< 3, int > I2_data_2 ; 
	t3.get_row( i1, i3, I2_data_2 );
	
	I1_data.set( 2 );
	test_I1_data.set( 2, 2 );
	t3.set_column( i2, i3, I1_data );
	vmml::vector< 2, int > I1_data_2 ; 
	t3.get_column( i2, i3, I1_data_2 );
	
	I3_data.set( 3 );
	test_I3_data.set( 3, 3, 3, 3 );
	t3.set_tube( i1, i2, I3_data );
	vmml::vector< 4, int > I3_data_2 ; 
	t3.get_tube( i1, i2, I3_data_2 );

	if (I2_data_2 == test_I2_data && I1_data_2 == test_I1_data && I3_data_2 == test_I3_data)
	{	
		log( "set_n_vector/set_row/set_column/set_tube()", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "I2_vector (1, 1, 1): " << I2_data << std::endl
		<< "I1_vector (2, 2): " << I1_data << std::endl
		<< "I3_vector (3, 3, 3, 3): " << I3_data << std::endl;
		log_error( error.str() );
	}
	
	
	//test get_slice() functions
	//frontal slice
	matrix< 2, 3, int > mat_frontal;
	matrix< 2, 3, int > test_mat_frontal;
	int data2[] = { 13, 14, 15, 16, 17, 3 };
	test_mat_frontal.set(data2, data2 + 6);
	t3.get_frontal_slice_fwd( 2, mat_frontal );

	matrix< 2, 3, int > mat_frontal_2;
	matrix< 2, 3, int > test_mat_frontal_2;
	test_mat_frontal_2.fill(7);
	t3.set_frontal_slice_fwd( 2, test_mat_frontal_2);
	t3.get_frontal_slice_fwd( 2, mat_frontal_2);

	if (mat_frontal == test_mat_frontal && mat_frontal_2 == test_mat_frontal_2)
	{	
		log( "get/set_frontal_slice_fwd() (I2xI1)", true  );
	} else
	{
		std::stringstream error;
		error 
		     << "after get_frontal_slice_fwd at i3 = 2: " << mat_frontal << std::endl
		     << "after set_frontal_slice_fwd after i3 = 2: " << mat_frontal_2 << std::endl;
		log_error( error.str() );
	}

	
	t3.set_frontal_slice_fwd( 2, test_mat_frontal);
	matrix< 3, 2, int > mat_frontal_bwd = transpose( mat_frontal );
	matrix< 3, 2, int > test_mat_frontal_bwd  = transpose( test_mat_frontal );
	t3.get_frontal_slice_bwd( 2, mat_frontal_bwd );
	
	matrix< 3, 2, int > mat_frontal_2_bwd  = transpose ( mat_frontal );;
	matrix< 3, 2, int > test_mat_frontal_2_bwd = transpose( test_mat_frontal_2 );
	t3.set_frontal_slice_bwd( 2, test_mat_frontal_2_bwd);
	t3.get_frontal_slice_bwd( 2, mat_frontal_2_bwd);
	
	
	if (mat_frontal_bwd == test_mat_frontal_bwd && mat_frontal_2_bwd == test_mat_frontal_2_bwd)
	{	
		log( "get/set_frontal_slice_bwd() (I1xI2)", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_frontal_slice_bwd at i3 = 2 is: " << mat_frontal_bwd << std::endl
		<< "after get_frontal_slice_bwd at i3 = 2 should be: " << test_mat_frontal_bwd << std::endl
		<< "after set_frontal_slice_bwd after i3 = 2 is: " << mat_frontal_2_bwd << std::endl
		<< "after set_frontal_slice_bwd after i3 = 2 should be: " << test_mat_frontal_2_bwd << std::endl;
		log_error( error.str() );
	}
	
	
	
	//lateral slice
	matrix< 2, 4, int > mat_lateral;
	matrix< 2, 4, int > test_mat_lateral;
	int data3[] = { 1, 7, 7, 19, 4, 10, 7, 1 };
	test_mat_lateral.set(data3, data3 + 8);
	t3.get_lateral_slice_bwd( 0, mat_lateral );
	
	matrix< 2, 4, int > mat_lateral_2;
	matrix< 2, 4, int > test_mat_lateral_2;
	test_mat_lateral_2.fill(6);
	t3.set_lateral_slice_bwd( 0, test_mat_lateral_2 );
	t3.get_lateral_slice_bwd( 0, mat_lateral_2 );
	
	if (mat_lateral == test_mat_lateral && mat_lateral_2 == test_mat_lateral_2)
	{	
		log( "get/set_lateral_slice_bwd() (I1xI3)", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_lateral_slice_bwd i2 = 0: " << mat_lateral << std::endl
		<< "after set_lateral_slice_bwd i2 = 0: " << mat_lateral_2 << std::endl;
		log_error( error.str() );
	}
	
	
	t3.set_lateral_slice_bwd( 0, test_mat_lateral );
	matrix< 4, 2, int > mat_lateral_fwd = transpose( mat_lateral );
	matrix< 4, 2, int > test_mat_lateral_fwd  = transpose( test_mat_lateral );
	t3.get_lateral_slice_fwd( 0, mat_lateral_fwd );
	
	matrix< 4, 2, int > mat_lateral_2_fwd = transpose( mat_lateral_2 );
	matrix< 4, 2, int > test_mat_lateral_2_fwd  = transpose( test_mat_lateral_2 );
	t3.set_lateral_slice_fwd( 0, test_mat_lateral_2_fwd );
	t3.get_lateral_slice_fwd( 0, mat_lateral_2_fwd );

	if (mat_lateral_fwd == test_mat_lateral_fwd && mat_lateral_2_fwd == test_mat_lateral_2_fwd)
	{	
		log( "get/set_lateral_slice_fwd() (I3xI1)", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_lateral_slice_fwd i2 = 0 is: " << mat_lateral_fwd << std::endl
		<< "after get_lateral_slice_fwd i2 = 0 should be: " << test_mat_lateral_fwd << std::endl
		<< "after set_lateral_slice_fwd i2 = 0 is: " << mat_lateral_2_fwd << std::endl
		<< "after set_lateral_slice_fwd i2 = 0 should be: " << test_mat_lateral_2_fwd << std::endl;
		log_error( error.str() );
	}

	
	//horizontal slice
	matrix< 3, 4, int > mat_horizontal;
	matrix< 3, 4, int > test_mat_horizontal;
	//int data4[] = { 6, 5, 3, 6, 11, 3, 6, 7, 7, 6, 1, 3 };
	int data4[] = { 6, 6, 6, 6, 5, 11, 7, 1, 3, 3, 7, 3 };
	test_mat_horizontal.set(data4, data4 + 12);
	t3.get_horizontal_slice_fwd( 1, mat_horizontal );
	
	matrix< 3, 4, int > mat_horizontal_2;
	matrix< 3, 4, int > test_mat_horizontal_2;
	test_mat_lateral_2.fill(5);
	t3.set_horizontal_slice_fwd( 1, test_mat_horizontal_2 );
	t3.get_horizontal_slice_fwd( 1, mat_horizontal_2 );
	
	if (mat_horizontal == test_mat_horizontal &&  mat_horizontal_2 == test_mat_horizontal_2 )
	{	
		log( "get/set_horizontal_slice_fwd() (I2xI3)", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_horizontal_slice_fwd at i1 = 0: " << mat_horizontal << std::endl
		<< "after set_horizontal_slice_fwd at i1 = 0: " << mat_horizontal_2 << std::endl;
		log_error( error.str() );
	}
	
	t3.get_horizontal_slice_fwd( 1, test_mat_horizontal );
	matrix< 4, 3, int > mat_horizontal_bwd  = transpose( mat_horizontal );
	matrix< 4, 3, int > test_mat_horizontal_bwd = transpose( test_mat_horizontal );
	t3.get_horizontal_slice_bwd( 1, mat_horizontal_bwd );
	
	matrix< 4, 3, int > mat_horizontal_2_bwd  = transpose( mat_horizontal_2 );
	matrix< 4, 3, int > test_mat_horizontal_2_bwd = transpose( test_mat_horizontal_2 );
	t3.set_horizontal_slice_bwd( 1, test_mat_horizontal_2_bwd );
	t3.get_horizontal_slice_bwd( 1, mat_horizontal_2_bwd );

	if (mat_horizontal_bwd == test_mat_horizontal_bwd &&  mat_horizontal_2_bwd == test_mat_horizontal_2_bwd )
	{	
		log( "get/set_horizontal_slice_bwd() (I3xI2)", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_horizontal_slice_bwd at i1 = 0 is: " << mat_horizontal_bwd << std::endl
		<< "after get_horizontal_slice_bwd at i1 = 0 should be: " << test_mat_horizontal_bwd << std::endl
		<< "after set_horizontal_slice_bwd at i1 = 0 is: " << mat_horizontal_2_bwd << std::endl
		<< "after set_horizontal_slice_bwd at i1 = 0 should be: " << test_mat_horizontal_2_bwd << std::endl;
		log_error( error.str() );
	}
	
	
	
	//CP reconstruction
	
	tensor3< 4,4,4, double > t3_cp_input;
	double data_in_cp[] = { 
		0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858, 
		0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
		0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
		0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071 };
	t3_cp_input.set(data_in_cp, data_in_cp + 64);
	
	vector< 4, double> lambda; lambda.at(0) = 3.4996; lambda.at(1) = 1.6186; lambda.at(2) = 1.0947; lambda.at(3) = 0.80494;
	double data_u1_cp[] = { 
		0.414514672613941, -0.273871195562966, 0.622972704020820, -0.106663870725802, 
		0.373815968628613, -0.303860590002700, 0.583235093038218, 0.038436291238145,
		0.307288348024630, -0.134400219672316, 0.462215931211866, -0.019751415378651, 
		0.770722439633298, 0.902551877199397, 0.241035825610295, 0.993355601880911};
	matrix< 4, 4, double > u1_cp; u1_cp.set( data_u1_cp, data_u1_cp + 16);
	double data_u2_cp[] = { 
		0.363947543028453, 0.298190617522100, 0.309440303599311, 0.906460391523818, 
		0.582584049017999, -0.442905614843157, -0.385671235746994, 0.178817966899174, 
		0.700932982871312, -0.827455897125261, -0.577509589347220, -0.345891330880255, 
		0.191913952792576, -0.173878435335793, 0.649605319126744, -0.163440755426602 };
	matrix< 4, 4, double > u2_cp; u2_cp.set( data_u2_cp, data_u2_cp + 16);
	double data_u3_cp[] = { 
		0.528860463703445, -0.453374574314148, 0.081548921839262, 0.026850375225737, 
		0.421141692804753, 0.103741153724736, 0.379185048676421, 0.600983237030590,
		0.588063451217321, 0.731736212173060, 0.724001032600722, -0.545272439076217, 
		0.443990610098116, 0.498248315786364, 0.570430518994385, -0.583760373220334};
	matrix< 4, 4, double > u3_cp; u3_cp.set( data_u3_cp, data_u3_cp + 16);
	
	tensor3< 4,4,4, double > t3_cp_reco_check;
	double data_out_cp[] = { 
		0.354263061197741, 0.336071353932582, 0.340120155317465, 0.148789058227363,
		0.335156282115524, 0.284369809600727, 0.270075906182746, 0.127689570299207,
		0.248778281916907, 0.271657788743462, 0.293346312353042, 0.118870461584074,
		0.347767704804178, 1.119919940137974, 1.528048459613556, 0.399390301753092,

		0.241873571540868, 0.267325817660143, 0.334780629072782, 0.301651586726927,
		0.277064331887503, 0.253524404359242, 0.282150329524853, 0.268828551075853,
		0.208805038309086, 0.198138779535580, 0.228624285742602, 0.217032724171575,
		0.925152576629409, 0.641979031363691, 0.446794524932619, 0.178096907232692,

		0.408959960261278, 0.458603477148098, 0.565021279854640, 0.533195254992912,
		0.300410894904167, 0.426297991808946, 0.575920757418752, 0.513247731493541,
		0.303903577435438, 0.299195909145241, 0.360430286717156, 0.385593403984964,
		0.559935740903585, 0.298948990808052, 0.267704058359116, 0.313877703199338,

		0.334349191835064, 0.331983603647059, 0.392220714079066, 0.406517347684063,
		0.234640418580267, 0.303232696231040, 0.405821982955168, 0.393611395358567,
		0.239172809919676, 0.216514108351268, 0.254464520462408, 0.296451592996035,
		0.276358969906296, 0.233761715518065, 0.311620024949421, 0.277319482627673 };
				
	t3_cp_reco_check.set(data_out_cp, data_out_cp + 64);
	
	matrix< 4, 4, double > v1; 
	matrix< 4, 4, double > v2; 
	matrix< 4, 4, double > v3; 
	u1_cp.transpose_to( v1 );
	u2_cp.transpose_to( v2 );
	u3_cp.transpose_to( v3 );
			
	tensor3< 4, 4, 4, double > t3_cp_reco;
	matrix< 4, 16, double > temp;
	t3_cp_reco.reconstruct_CP( lambda,  v1, v2, v3, temp);
													
	precision = 1.0e-4;
	if ( t3_cp_reco.equals( t3_cp_reco_check, precision) )
	{	
		log( "tensor3 CP reconstruction", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "tensor3 CP reconstruction: is " << std::endl << t3_cp_reco
		<< std::endl << "should be: " << t3_cp_reco_check
		<< std::endl;
		log_error( error.str() );
	}
	
	//innerproduct between a tensor and a cp decomposition
	
	double innerp = t3_cp_input.tensor_inner_product( lambda, u1_cp, u2_cp, u3_cp);
	
	double innerp_check = 11.5230085;
	if ( (innerp - innerp_check) < precision )
	{	
		log( "tensor3 inner product t3 and cp of t3", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "tensor3 inner product t3 and cp of t3: is " << std::endl << innerp
		<< std::endl << "should be: " << innerp_check
		<< std::endl;
		log_error( error.str() );
	}
	
	//unfolding
	//matrix< I3, I1*I2, T> m_horizontal;
	//matrix< I1, I2*I3, T> m_lateral;
	//matrix< I2, I1*I3, T> m_frontal;
	matrix< 4, 6, int> m_horizontal;
	matrix< 2, 12, int> m_lateral;
	matrix< 3, 8, int> m_frontal;
	
	t3.fill_increasing_values();
	t3.horizontal_unfolding_bwd( m_horizontal);
	t3.lateral_unfolding_bwd( m_lateral);
	t3.frontal_unfolding_bwd(  m_frontal);
	
	matrix< 4, 6, int> m_horizontal_test;
	int data5[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
	m_horizontal_test.set(data5, data5 + 24);
	matrix< 2, 12, int> m_lateral_test;
	int data6[] = { 0, 6, 12, 18, 1, 7, 13, 19, 2, 8, 14, 20, 3, 9, 15, 21, 4, 10, 16, 22, 5, 11, 17, 23 };
	m_lateral_test.set(data6, data6 + 24);
	matrix< 3, 8, int> m_frontal_test;
	int data7[] = { 0, 3, 6, 9, 12, 15, 18, 21, 1, 4, 7, 10, 13, 16, 19, 22, 2, 5, 8, 11, 14, 17, 20, 23 };
	m_frontal_test.set(data7, data7 + 24);
	
	if ( m_horizontal_test == m_horizontal && m_lateral_test == m_lateral && m_frontal_test == m_frontal )
	{	
		log( "backward unfolding along all modes", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "backward unfolding: " << std::endl
		<< "unfolding_horizontal should be: " << std::endl << m_horizontal_test << std::endl
		<< "unfolding_horizontal is: " << std::endl << m_horizontal << std::endl
		<< "unfolding_lateral should be: " << std::endl << m_lateral_test << std::endl
		<< "unfolding_lateral is: " << std::endl << m_lateral << std::endl
		<< "unfolding_frontal should be: " << std::endl << m_frontal_test << std::endl
		<< "unfolding_frontal is: " << std::endl << m_frontal << std::endl
		<< std::endl;
		log_error( error.str() );
	}
	
	
	//compute frobenius norm of a tensor3
	double f_norm_check = 65.75712889109438;
	
	t3.fill_increasing_values();
	double f_norm = t3.frobenius_norm();
	
	ok = f_norm == f_norm_check;
	
	//compute frobenius norm of the difference of a tensor3 and another tensor3
	f_norm_check = 10;
	
	t3.fill_increasing_values();
	t3_tmp.fill_increasing_values();
	t3_tmp.at(1,1,1) = 0;
	
	f_norm = t3.frobenius_norm( t3_tmp );
	
	ok = (f_norm == f_norm_check ) && ok;

	if ( ok )
	{	
		log( "compute frobenius norm", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "compute frobenius norm: should be: " << f_norm_check << " is: " << f_norm << std::endl;
		log_error( error.str() );
	}

	//compute rmse between two tensor3
	double rmse_check = 2.041241452319315;
	
	t3.fill_increasing_values();
	t3_tmp.fill_increasing_values();
	t3_tmp.at(1,1,1) = 0;
	
	double rmse = t3.rmse( t3_tmp );
	ok = rmse == rmse_check;
	
	//psnr
	double psnr_check = 48.43872076218154;
	double psnr = t3.compute_psnr( t3_tmp, t3.get_max() );
	ok = ok && (psnr == psnr_check );
	
	//avg. frob. norm difference
	double afn_check = 13.42261772780059;
	double afn = t3.avg_frobenius_norm();
	ok = ok && ( afn == afn_check );
	
	
	if ( ok )
	{	
		log( "compute RMSE, PSNR, averaged frob.norm difference ", ok  );
	} else
	{
		std::stringstream error;
		error 
		<< "compute RMSE: should be: " << rmse_check << " is: " << std::setprecision(16) << rmse << std::endl
		<< "compute PSNR: should be: " << psnr_check << " is: " << std::setprecision(16) << psnr << std::endl
		//<< "compute frob.norm is: " << t3.frobenius_norm(); << ", size() " << t3.size() << std::endl
		<< "compute avg.frob.norm : should be: " << afn_check << " is: " << std::setprecision(16) << afn << std::endl;
		log_error( error.str() );
	}
	
	
	//set diagonal values in a cubic tensor3, i.e., R=I1, R=I2, R=I3
	tensor3< 3, 3, 3, int >  t3_diag;
	tensor3< 3, 3, 3, int >  t3_diag_check;
	t3_diag_check.zero(); t3_diag_check.at(1, 1, 1) = 1; t3_diag_check.at(2, 2, 2) = 2;
	vector< 3, int > diag_values;
	diag_values.at(0) = 0; diag_values.at(1) = 1; diag_values.at(2) = 2;
	t3_diag.diag( diag_values );

	
	if ( t3_diag == t3_diag_check )
	{	
		log( "fill diagonal values", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "fill diagonal values: should be: " << std::endl << t3_diag_check << std::endl 
		<< " is: " << std::endl << t3_diag << std::endl;
		log_error( error.str() );
	}
	
	{	//tensor3 type cast
		tensor3< 2, 3, 4, double >  t3_type_a;
		tensor3< 2, 3, 4, int >  t3_type_b;
		tensor3< 2, 3, 4, int >  t3_type_b_check;
		
		t3_type_a.fill(2.4); 
		t3_type_b_check.fill(2);
		t3_type_b.cast_from( t3_type_a );
		
		if (t3_type_b_check == t3_type_b )
		{	
			log( "type cast ", true  );
		} else
		{
			std::stringstream error;
			error << "type cast - tensor3 type double: " << std::endl << t3_type_a << std::endl
			<< " tensor3 type int should be: " << std::endl << t3_type_b_check << std::endl
			<< " is: " << t3_type_b << std::endl;
			log_error( error.str() );
		}
	}
	{
		//tensor3: from float_t to uint_t
		tensor3< 2, 3, 4, double >  t3_type_a;
		tensor3< 2, 3, 4, unsigned char >  t3_type_b;
		tensor3< 2, 3, 4, unsigned char >  t3_type_b_check;
		
		t3_type_a.fill(2.4); 
		t3_type_a.at(0,0,0) = 2.67; 
		t3_type_a.at(0,2,0) = 2.67; 
		t3_type_b_check.fill(2);
		t3_type_b_check.at(0,0,0) = 3; 
		t3_type_b_check.at(0,2,0) = 3; 
		t3_type_b.float_t_to_uint_t( t3_type_a );
		
		if (t3_type_b_check == t3_type_b )
		{	
			log( "from float_t to uint_t ", true  );
		} else
		{
			std::stringstream error;
			error << "from float_t to uint_t - tensor3 type double: " << std::endl << t3_type_a << std::endl
			<< " tensor3 type int should be: " << std::endl << t3_type_b_check << std::endl
			<< " is: " << t3_type_b << std::endl;
			log_error( error.str() );
		}
	}
	
	//export
	tensor3< 2, 3, 4, double >  t3_export_import;
	t3_export_import.fill_increasing_values();
	std::vector< double > export_data;
	t3_export_import.export_to( export_data );
	
	double export_data_check[] = { 0, 3, 1, 4, 2, 5, 6, 9, 7, 10, 8, 11, 12, 15, 13, 16, 14, 17, 18, 21, 19, 22, 20, 23 };
	precision = 1.0e-1;
	ok = true;
	for (int i = 0; i < 24 && ok; ++i )
	{
		ok = (abs(export_data_check[i]) - abs(export_data[i])) < precision;
	}
	
	log( "export tensor3", ok  );
	
	//import tucker3 from vector
	std::vector< double > in_data = export_data;
	tensor3< 2, 3, 4, double >  t3_import_check;
	t3_import_check.fill_increasing_values();
	t3_export_import.zero();
	t3_export_import.import_from( in_data );
		
	if ( t3_export_import.equals( t3_import_check, precision ) )	{	
		log( "import tensor3" , true  );
	} else
	{
		std::stringstream error;
		error 
		<< "import tensor3: " << std::endl
		<< "import tensor3 should be: " << std::endl << t3_import_check << std::endl
		<< "is: " << std::endl << t3_export_import << std::endl;		
		
		log_error( error.str() );
	}		
	
	//get_min() + get_max()
	
	tensor3< 3, 2, 5, int >  t3_get_min_max;
	t3_get_min_max.fill_increasing_values();
	
	int t3_min = t3_get_min_max.get_min();
	int t3_max = t3_get_min_max.get_max();
	
	if ( (t3_min == 0) && (t3_max == 29) )	{	
		log( "get min/max" , true  );
	} else
	{
		std::stringstream error;
		error 
		<< "get min/max: " << std::endl
		<< "min should be: " << std::endl << 0 << " is: " << t3_min << std::endl 
		<< "max should be: " << std::endl << 29 << " is: " << t3_max << std::endl;
		
		log_error( error.str() );
	}		
	
	//quantize
	{
		tensor3< 2, 4, 3, float >  t3_raw;
		t3_raw.fill(0.45692);
		t3_raw.at( 0,2,2) = 0.67777; t3_raw.at(1,0,1) = 0.111111; t3_raw.at(1,2,0) = -0.23; t3_raw.at(1,3,0) = -0.99;
		t3_raw.at(1,0,0) = -0.8; t3_raw.at(1,2,1) = 0.0; t3_raw.at(0,3,2) = 0.99; t3_raw.at(0,1,0) = 0.23;
		tensor3< 2, 4, 3, unsigned char >  t3_quant; t3_quant.zero();
		tensor3< 2, 4, 3, unsigned char >  t3_quant_check; 
		int data_unsigned[] = {186, 157, 186, 186, 24, 186, 98, 0, 186, 186, 186, 186, 142, 186, 128, 186, 186, 186, 215, 255, 186, 186, 186, 186};
		t3_quant_check.set(data_unsigned, data_unsigned+24);
		
		float min_value = 50;
		float max_value = -50;
		
		t3_raw.quantize( t3_quant, min_value, max_value );
		
		tensor3< 2, 4, 3, char >  t3_quant_sign; t3_quant_sign.zero();
		tensor3< 2, 4, 3, char >  t3_quant_sign_check; 
		int data_signed[] = { 59, 30, 59, 59, -102, 59, -29, -127, 59, 59, 59, 59, 14, 59, 0, 59, 59, 59, 87, 127, 59, 59, 59, 59 };
		t3_quant_sign_check.set(data_signed, data_signed +24 );
		
		t3_raw.quantize( t3_quant_sign, min_value, max_value );
		
		
		//dequantize
		tensor3< 2, 4, 3, float >  t3_dequant;
		tensor3< 2, 4, 3, float >  t3_dequant_sign;
		
		t3_quant.dequantize( t3_dequant, min_value, max_value );
		t3_quant_sign.dequantize( t3_dequant_sign, min_value, max_value );
		
		ok = ( t3_quant_check == t3_quant ) && ( t3_quant_sign_check == t3_quant_sign ) && t3_dequant.equals(t3_dequant_sign, 0.01);
#if 0
		std::cout << " quantization: is " << ok << std::endl 
		<< "original: " << t3_raw << std::endl
		<< "linear: " << t3_quant << std::endl
		<< "linear signed: " << t3_quant_sign << std::endl
		<< "deq. from linear: " << std::endl << t3_dequant << std::endl;
#endif	
		
		//logarithmic quantization
		float lmin_value = 50;
		float lmax_value = -50;
		tensor3< 2, 4, 3, unsigned char >  t3_quant_log; t3_quant.zero();
		tensor3< 2, 4, 3, float >  t3_raw2(t3_raw);
		t3_raw2.at(0,2,1) = 200;
		tensor3< 2, 4, 3, unsigned char >  t3_quant2; t3_quant2.zero();
		t3_raw2.quantize( t3_quant2, lmin_value, lmax_value );
		tensor3< 2, 4, 3, float >  t3_dequant2;
		t3_quant2.dequantize( t3_dequant2, lmin_value, lmax_value );
		
		unsigned char tt_range = 127;
		tensor3< 2, 4, 3, char >  signs;
		t3_raw2.quantize_log( t3_quant_log, signs, lmin_value, lmax_value, tt_range );
		tensor3< 2, 4, 3, float >  t3_dequant_log;
		t3_quant_log.dequantize_log( t3_dequant_log, signs, lmin_value, lmax_value );
		
		float deq_log_check[] = {
			0.456192, 0.232188, 0.456192, 0.456192,
			-0.794302, 0.456192, -0.232188, -0.950592,
			 0.456192, 0.456192, 200, 0.456192, 
			 0.13346, 0.456192, 0, 0.456192,
			 0.456192, 0.456192, 0.650535, 0.950592,
			 0.456192, 0.456192, 0.456192, 0.456192
		};
		tensor3< 2, 4, 3, float >  t3_dequant_log_check;
		t3_dequant_log_check.set(deq_log_check, deq_log_check +24 );
		
		ok = ok && t3_dequant_log_check.equals( t3_dequant_log, 0.001 );
		
#if 0
		std::cout << " quantization: is " << ok << std::endl 
		<< "original: " << t3_raw2 << std::endl
		<< "linear: " << t3_quant2 << std::endl
		<< "log-scale: " << std::endl << t3_quant_log << std::endl
		<< "signs (0 = neg, 1 = pos): " << std::endl << signs << std::endl
		<< "deq. from log-scale: " << std::endl << t3_dequant_log << std::endl
		<< "deq. from linear: " << std::endl << t3_dequant2 << std::endl;
#endif	

		//linear quantization with separate sign encoding
		tensor3< 2, 4, 3, float >  t3_raw3;
		signs.zero(); t3_quant.zero();
		t3_raw.quantize_to( t3_quant, signs, min_value, max_value, tt_range );
		
		t3_quant.dequantize( t3_raw3, signs, min_value, max_value );

		float t3_raw3_data[] = {
			0.459921, 0.233858, 0.459921, 0.459921,
			-0.802913, 0.459921, -0.233858, -0.99,
			0.459921, 0.459921, 0.459921, 0.459921,
			0.109134, 0.459921, 0, 0.459921,
			0.459921, 0.459921, 0.678189, 0.99,
			0.459921, 0.459921, 0.459921, 0.459921
		};
		tensor3< 2, 4, 3, float >  t3_raw3_check;
		t3_raw3_check.set( t3_raw3_data, t3_raw3_data +24 );
		
		ok = ok && t3_raw3.equals( t3_raw3_check, 0.001 );
		
		if ( ok )	{	
			log( "quantize/dequantize" , ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "quantize/dequantize " << std::endl
			<< "raw is: " << std::endl << t3_raw << std::endl 
			<< "unsigned quantized is: " << std::endl << t3_quant << std::endl
			<< "signed quantized is: " << std::endl << t3_quant_sign << std::endl
			<< "dequantized unsigned is: " << std::endl << t3_dequant << std::endl 
			<< "dequantized signed is: " << std::endl << t3_dequant_sign << std::endl 
			<< "min_value: " << min_value << " max_value: " << max_value << std::endl;
			
			log_error( error.str() );
		}	
	}
	{
		//number of nonzeros
		
		tensor3< 4, 5, 6, int > t3_nnz;
		t3_nnz.fill_increasing_values();
		t3_nnz.at( 3,3,3) = 0; t3_nnz.at( 2,3,3) = -4;
		size_t number_nonzeros = t3_nnz.nnz();
				
		tensor3< 4, 4, 4, float > t3_nnz2;
		t3_nnz2.fill( 0.9878);
		t3_nnz2.at( 3,3,3) = 0; t3_nnz2.at( 2,3,3) = -1; t3_nnz2.at( 2,2,3) = 0.045; t3_nnz2.at( 1,2,3) = -0.085;
		t3_nnz2.at( 0,2,3) = 0.00000035; t3_nnz2.at( 0,1,3) = -0.00000035;
		size_t number_nonzeros2 = t3_nnz2.nnz( 0.00001 );
		
		ok = ( number_nonzeros == 118 ) && (number_nonzeros2 == 61);
		log( "get number of nonzeros" , ok  );
		
	}
	
	{
		//thresholding
		
		tensor3< 4, 5, 6, float > t3_thresh;
		t3_thresh.fill(0.5673);
		t3_thresh.at( 3,3,3) = 0; t3_thresh.at( 2,3,3) = -1; t3_thresh.at( 2,2,3) = 0.045; t3_thresh.at( 1,2,3) = -0.085;
		t3_thresh.at( 0,2,3) = 0.00000035; t3_thresh.at( 0,1,3) = -0.00000035;
		t3_thresh.at( 3,3,3) = 0; t3_thresh.at( 2,3,3) = 0.001; t3_thresh.at( 0,3,3) = 0.00001;
				
		tensor3< 4, 5, 6, float > t3_thresh2(t3_thresh);
		t3_thresh.threshold( 0.001f );
		size_t number_nonzeros = t3_thresh.nnz();
		
		tensor3< 4, 5, 6, unsigned char > t3_thresh_char;
		t3_thresh_char.fill( 6);
		t3_thresh_char.at(0,0,0) = 3;
		t3_thresh_char.threshold( 4 );
		size_t number_nonzeros_char = t3_thresh_char.nnz();
		
		ok = (number_nonzeros == 115 ) && (number_nonzeros_char == 119);
		log( "thresholding" , ok  );
	}
	
	{
		//get_sub_tensor3
		
		tensor3< 4, 5, 6, int > t3_bigger;
		t3_bigger.fill_increasing_values();
		
		tensor3< 2, 3, 2, int > t3_sub_check;
		int data_sub_check[] = {32, 33, 34, 37, 38, 39, 52, 53, 54, 57, 58, 59};
		t3_sub_check.set( data_sub_check, data_sub_check + 12);
		
		tensor3< 2, 3, 2, int > t3_sub;
		t3_bigger.get_sub_tensor3( t3_sub, 2, 2, 1 );
		ok = t3_sub == t3_sub_check;
		
		t3_sub.fill_increasing_values();
		t3_bigger.set_sub_tensor3( t3_sub, 0, 1, 1 );
		
		tensor3< 4, 5, 6, int > t3_bigger_check;
		t3_bigger_check = t3_bigger;
		t3_bigger_check.at(0, 1, 1 ) = 0; t3_bigger_check.at(0, 2, 1 ) = 1; t3_bigger_check.at(0, 3, 1 ) = 2;
		t3_bigger_check.at(1, 1, 1 ) = 3; t3_bigger_check.at(1, 2, 1 ) = 4; t3_bigger_check.at(1, 3, 1 ) = 5;
		t3_bigger_check.at(0, 1, 2 ) = 6; t3_bigger_check.at(0, 2, 2 ) = 7; t3_bigger_check.at(0, 3, 2 ) = 8;
		t3_bigger_check.at(1, 1, 2 ) = 9; t3_bigger_check.at(1, 2, 2 ) = 10; t3_bigger_check.at(1, 3, 2 ) = 11;
		
		ok = ok && ( t3_bigger == t3_bigger_check);
		if ( ok )	{	
			log( "get/set sub tensor3" , ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "get/set sub tensor3 " << std::endl
			<< "sub tensor3 is: " << std::endl << t3_sub << std::endl 
			<< "should be: " << std::endl << t3_sub_check << std::endl;			
			log_error( error.str() );
		}	
	}


    {   // operator(i,j,k)
    	tensor3< 4, 3, 4, int >  t3x;

        t3x.zero();
        t3x.get_frontal_slice_fwd( 0 )( 1, 2 ) = 12;
        t3x.get_frontal_slice_fwd( 0 )( 2, 1 ) = 21;

        if ( t3x( 1, 2, 0 ) == 12 && t3x( 2, 1, 0 ) == 21 )
        {
            log( "operator( i1, i2, i3 )", true  );
        }
        else
        {
            std::stringstream error;
            error 
                << "operator(i1,i2,i3) failed. Tensor: " << std::endl << t3x
                << std::endl;
            log_error( error.str() );
        }
    }

    {   // get_array_ptr
    	tensor3< 4, 3, 4, int >  t3x;
        
        t3x.zero();
        t3x.get_frontal_slice_fwd( 0 )( 0, 0 ) = 23;
        t3x.get_frontal_slice_fwd( 0 )( 1, 0 ) = 12;
        
        int* array = t3x.get_array_ptr();
        
        if ( array[ 0 ] == 23 && array[ 1 ] == 12 )
        {
            log( "get_array_ptr", true  );
        }
        else
        {
            std::stringstream error;
            error 
            << "get_array_ptr() failed. Tensor: \n" << t3x
            << std::endl;
            log_error( error.str() );
        }
    }


    {   // get_array_ptr
    
    	tensor3< 4, 3, 4, int >  t3x;
        
        t3x.zero();
        t3x.get_frontal_slice_fwd( 0 )( 0, 0 ) = 23;
        t3x.get_frontal_slice_fwd( 0 )( 1, 0 ) = 12;
        t3x.get_frontal_slice_fwd( 0 )( 0, 1 ) = 11;
        t3x.get_frontal_slice_fwd( 2 )( 1, 1 ) = 13;
        
            
        tensor3< 4, 3, 4, int >::iterator 
            it      = t3x.begin(),
            it_end  = t3x.end();
        
        bool ok = *it == 23;
        if ( ! ok ) 
        {
            std::cout << "*it should be " << 23 << " but is " << *it << std::endl;
        }
        
        for( size_t index = 0; index < 1; ++it, ++index ) {}
        
        if ( ok ) 
            ok = *it == 12;

        if ( ! ok ) 
        {
            std::cout << "*it should be " << 12 << " but is " << *it << std::endl;
        }

        for( size_t index = 0; index < 3; ++it, ++index ) {}
        
        if ( ok ) 
            ok = *it == 11;
        
        if ( ! ok ) 
        {
            std::cout << "*it should be " << 11 << " but is " << *it << std::endl;
        }
        
        for( size_t index = 0; index < 8 + 12 + 5; ++it, ++index ) {}
        
        if ( ok ) 
            ok = *it == 13;

        if ( ! ok ) 
        {
            std::cout << "*it should be " << 13 << " but is " << *it << std::endl;
        }
        
        if ( ok )
        {
            log( "iterator", true  );
        }
        else
        {
            std::stringstream error;
            error 
            << "iterator test failed. Tensor: \n" << t3x
            << std::endl;
            log_error( error.str() );
        }
    }
	
	
	{
		
		
		unsigned int data_uct[] = { 
			0, 27, 51, 8, 26, 59, 64, 53,
			20, 66, 15, 54, 96, 219, 28, 4,
			39, 17, 17, 187, 226, 199, 33, 14,
			70, 39, 136, 58, 159, 6, 209, 23,
			5, 165, 86, 83, 15, 54, 230, 64,
			33, 70, 126, 155, 13, 148, 50, 31,
			7, 59, 181, 72, 236, 13, 13, 63,
			27, 2, 5, 4, 1, 10, 66, 40,
			40, 38, 23, 20, 42, 41, 36, 74,
			 7, 74, 0, 27, 186, 146, 42, 61,
			 36, 60, 32, 220, 39, 106, 254, 63,
			 50, 182, 129, 114, 82, 67, 187, 9,
			 20, 134, 241, 206, 76, 70, 70, 62,
			 39, 102, 9, 27, 189, 251, 179, 0,
			 75, 34, 226, 185, 75, 163, 55, 2,
			 56, 39, 11, 59, 67, 73, 25, 49
		};
		
		tensor3< 8,8, 2, unsigned char > uct_t3_check;
		uct_t3_check.set( data_uct, data_uct + 128);
		
		tensor3< 8,8, 2, unsigned char > uct_t3;
		uct_t3.fill_random(8);
		
		//std::cout << "t3 is: " << std::endl << uct_t3 << std::endl << "remove uct cylinder, now t3 is" << std::endl;
		
		std::string dir = ".";
		std::string in_filename = "in.raw";
		std::string out_filename = "out.raw";
		uct_t3.write_to_raw( dir, in_filename );
		
		double sigma = uct_t3.stdev();
		
		t3_converter<8,8,8, unsigned char>::remove_uct_cylinder( dir, in_filename, out_filename, sigma, 0, 0, 2 ); //remove_uct_cylinder
		
		uct_t3.read_from_raw( dir, out_filename );
		
		ok = uct_t3 == uct_t3_check;
		
		remove("in.raw");
		remove("out.raw");
		
		if ( ok )	{	
			log( "t3_converter: remove uct ring" , ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "remove uct ring " << std::endl
			<< "t3 is: " << std::endl << uct_t3 << std::endl 
			<< "t3 should be: " << std::endl << uct_t3_check 
			<< std::endl;
			
			log_error( error.str() );
		}	
		
	}
	
	{
		//average 8voxels to 1voxel
		tensor3< 5, 5, 5, unsigned int > t3;
		t3.fill_increasing_values();
		tensor3< 2, 2, 2,  unsigned int > t3_sub_check;
		unsigned int t3_sub_data[] = { 16, 18, 26, 28, 66, 68, 76, 78 }; 
		t3_sub_check.set( t3_sub_data, t3_sub_data + 8 );
		
		tensor3< 2, 2, 2, unsigned int > t3_sub;
		t3.average_8to1( t3_sub );
						
		ok = ( t3_sub == t3_sub_check ) ;
		
		log( "subsample tensor3 (average 8 voxels to 1 voxel)" , ok  );
		
	}
	
	
	{
		//fill slice symmetric tensor
		tensor3< 3, 3, 2, unsigned char > t3;
		t3.fill_rand_sym_slices( 3 );
				
		tensor3< 3, 3, 3, unsigned char > t3_2;
		t3_2.fill_rand_sym( 6 );
		
		
		ok = t3.at( 0, 1, 0) == t3.at( 1, 0, 0 );
		ok = ok && (t3.at( 1, 2, 1 ) == t3.at( 2, 1, 1));
		
		unsigned char val_102 =  t3_2.at( 1, 0, 2 );
		unsigned char val_012 =  t3_2.at( 0, 1, 2 );
		unsigned char val_021 =  t3_2.at( 0, 2, 1 );
		unsigned char val_201 =  t3_2.at( 2, 0, 1 );
		unsigned char val_120 =  t3_2.at( 1, 2, 0 );
		unsigned char val_210 =  t3_2.at( 2, 1, 0 );
		
		ok = ok && val_102 == val_012 && val_102 == val_021 && val_102 == val_201; 
		ok = ok && val_102 == val_120 && val_102 == val_210;
		
		log( "fill tensor3 with symmetric random values" , ok  );
		
	}
	
	{
		//mean, variance and stdev
		tensor3< 3, 3, 2, unsigned char > t3;
		t3.fill_random( 3 );
		
		double mean_val = t3.mean();
		double var = t3.variance();
		double sigma = t3.stdev();
						
		ok = (mean_val - 76.0556 <= 0.01 );
		ok = ok && ((4469.232026 - var) <= 0.01);
		ok = ok && ((66.85231504 - sigma) <= 0.01);
		
		log( "mean, variance and standard deviation" , ok  );
		
	}
	
	{
		//load mmap for tensor3
		
		//create test data
		std::string dir = ".";
		std::string filename = "mmap_testdata.raw";
		tensor3< 4,4,4, unsigned char > t3_testdata;
		t3_testdata.fill_random( 3 );
		t3_testdata.write_to_raw( dir, filename );
		
		tensor3< 4,4,4, unsigned char > t3( dir, filename, false );
		
		unsigned int data_mmp_t3[] = { 
			0, 152, 9, 125,
			100, 167, 205, 26,
			68, 35, 38, 40,
			95, 9, 142, 150, 
			3, 64, 137, 63,
			5, 15, 148, 26,
			38, 195, 70, 186,
			51, 201, 245, 73,
			200, 228, 188, 243,
			36, 68, 241, 55,
			53, 248, 42, 228,
			251, 24, 66, 166,
			208, 181, 117, 182,
			78, 210, 176, 130,
			76, 19, 185, 139,
			110, 127, 46, 244
		};
		tensor3< 4,4,4, unsigned char > t3_check;
		t3_check.set( data_mmp_t3, data_mmp_t3 + 64);
		
		//double fnorm = t3.frobenius_norm();
		//std::cout << "frobnorm of mmapped file " << fnorm << std::endl;
		
		//tensor3< 4,4,4, unsigned char > t33(t3);
		//std::cout << "t33\n" << t33 << std::endl;

		//tensor3< 4,4,4, unsigned short > t333(t3);
		//std::cout << "t333\n" << t333 << std::endl;

		
		ok = t3_check == t3;
		log( "load tensor3 from memory mapped file" , ok  );
		
		remove("mmap_testdata.raw");
	}
	
	{
		
		//create test data
		std::string dir = ".";
		std::string in_filename = "in.raw";
		std::string out_filename = "out.raw";
		tensor3< 4,4,4, unsigned char > t3_in;
		t3_in.fill_random( 3 );
		t3_in.write_to_raw( dir, in_filename );
		
		t3_converter<4,4,4, unsigned char>::convert_raw<float>( dir, in_filename, out_filename );
		
		tensor3< 4,4,4, float > t3_out;
		t3_out.read_from_raw( dir, out_filename );

		ok = t3_out.equals( t3_in, 0.001);
		log( "t3 converter: convert raw" , ok  );
		
		
		remove("in.raw");
		remove("out.raw");
		
	}
		
	ok = true;
    return ok;
}


	
} // namespace vmml

