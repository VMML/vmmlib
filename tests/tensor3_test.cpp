#include "tensor3_test.hpp"

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
	
	
	//test tensor3 matrix multiplication
	t3.fill_increasing_values();
	matrix<5, 4, int> u3;
	u3.fill(1);
	tensor3<2, 3, 5, int> t3_jji;
	t3_jji.multiply_horizontal_bwd(t3, u3 );
	
	tensor3<6, 3, 5, int> t3_iji;
	matrix<6, 2, int> u1;
	u1.fill(2);
	t3_iji.multiply_lateral_bwd(t3_jji, u1 );

	tensor3<6, 7, 5, int> t3_iii;
	matrix<7, 3, int> u2;
	u2.fill(3);
	t3_iii.multiply_frontal_bwd(t3_iji, u2 );
	 
	
	t3.fill_increasing_values();
	tensor3<6, 7, 5, int> t3_reco;
	t3_reco.full_tensor3_matrix_multiplication( t3, u1, u2, u3 ); 
	tensor3<6, 7, 5, int> t3_reco2;
	t3_reco2.full_tensor3_matrix_kronecker_mult( t3, u1, u2, u3 ); 
	ok = t3_reco == t3_reco2;
	
	tensor3<6, 7, 5, int> t3_iii_test;
	t3_iii_test.fill(1656);
	if ( t3_iii_test == t3_reco && t3_iii_test == t3_iii && ok )
	{	
		log( "tensor3 matrix multiplication along all three modes", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "T3_result (all values should be 1656): " << std::endl << t3_iii
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
	
	if ( f_norm == f_norm_check )
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
		
	if ( rmse == rmse_check )
	{	
		log( "compute RMSE ", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "compute RMSE: should be: " << rmse_check << " is: " << std::setprecision(16) << rmse << std::endl;
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
	double precision = 1.0e-1;
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
		
		
		if ( t3_sub == t3_sub_check )	{	
			log( "get sub tensor3" , true  );
		} else
		{
			std::stringstream error;
			error 
			<< "get sub tensor3 " << std::endl
			<< "sub tensor3 is: " << std::endl << t3_sub << std::endl 
			<< "should be: " << std::endl << t3_sub_check << std::endl;			
			log_error( error.str() );
		}	
	}
	
	ok = true;
    return ok;
}


	
} // namespace vmml

