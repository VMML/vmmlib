#include "tensor3_test.hpp"

#include <stdint.h>
#include <vmmlib/tensor3.hpp>

#include <sstream>

namespace vmml
{

bool
tensor3_test::run()
{
    bool ok = false;
	        
	
	tensor3< 2, 3, 4, uint16_t >  t3;
    tensor3< 2, 3, 4, uint16_t >  t3_tmp;

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
	uint16_t data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };
	t3.set(data, data + 24);
	ok = false;
	uint16_t check_value = 1;
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
    tensor3< 2, 3, 4, uint16_t >  t3_2;
    tensor3< 2, 3, 4, uint16_t >  t3_3;
    tensor3< 2, 3, 4, uint16_t >  t3_result;
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
	uint16_t i1 = 1;
	uint16_t i2 = 2;
	uint16_t i3 = 3;
	t3.fill_increasing_values();
	t3 += 1;

	vmml::vector< 3, uint16_t > test_I2_data ; 
	vmml::vector< 3, uint16_t > I2_data ; 
	t3.get_row( i1, i3, I2_data );
	test_I2_data.set(22, 23, 24);
	
	vmml::vector< 2, uint16_t > test_I1_data ; 
	vmml::vector< 2, uint16_t > I1_data ; 
	t3.get_column( i2, i3, I1_data );
	test_I1_data.set(21, 24);
	
	vmml::vector< 4, uint16_t > test_I3_data ; 
	vmml::vector< 4, uint16_t > I3_data ; 
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
	vmml::vector< 3, uint16_t > I2_data_2 ; 
	t3.get_row( i1, i3, I2_data_2 );
	
	I1_data.set( 2 );
	test_I1_data.set( 2, 2 );
	t3.set_column( i2, i3, I1_data );
	vmml::vector< 2, uint16_t > I1_data_2 ; 
	t3.get_column( i2, i3, I1_data_2 );
	
	I3_data.set( 3 );
	test_I3_data.set( 3, 3, 3, 3 );
	t3.set_tube( i1, i2, I3_data );
	vmml::vector< 4, uint16_t > I3_data_2 ; 
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
	matrix< 2, 3, uint16_t > mat_frontal;
	matrix< 2, 3, uint16_t > test_mat_frontal;
	uint16_t data2[] = { 13, 14, 15, 16, 17, 3 };
	test_mat_frontal.set(data2, data2 + 6);
	t3.get_frontal_slice( 2, mat_frontal );

	matrix< 2, 3, uint16_t > mat_frontal_2;
	matrix< 2, 3, uint16_t > test_mat_frontal_2;
	test_mat_frontal_2.fill(7);
	t3.set_frontal_slice( 2, test_mat_frontal_2);
	t3.get_frontal_slice( 2, mat_frontal_2);

	if (mat_frontal == test_mat_frontal && mat_frontal_2 == test_mat_frontal_2)
	{	
		log( "get/set_frontal_slice()", true  );
	} else
	{
		std::stringstream error;
		error 
		     << "after get_frontal_slice at i3 = 2: " << mat_frontal << std::endl
		     << "after set_frontal_slice after i3 = 2: " << mat_frontal_2 << std::endl;
		log_error( error.str() );
	}
	
	
	//lateral slice
	matrix< 2, 4, uint16_t > mat_lateral;
	matrix< 2, 4, uint16_t > test_mat_lateral;
	uint16_t data3[] = { 1, 7, 7, 19, 4, 10, 7, 1 };
	test_mat_lateral.set(data3, data3 + 8);
	t3.get_lateral_slice( 0, mat_lateral );
	
	matrix< 2, 4, uint16_t > mat_lateral_2;
	matrix< 2, 4, uint16_t > test_mat_lateral_2;
	test_mat_lateral_2.fill(6);
	t3.set_lateral_slice( 0, test_mat_lateral_2 );
	t3.get_lateral_slice( 0, mat_lateral_2 );
	
	if (mat_lateral == test_mat_lateral && mat_lateral_2 == test_mat_lateral_2)
	{	
		log( "get/set_lateral_slice()", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_lateral_slice i2 = 0: " << mat_frontal << std::endl
		<< "after set_lateral_slice i2 = 0: " << mat_frontal_2 << std::endl;
		log_error( error.str() );
	}
	
	//horizontal slice
	matrix< 3, 4, uint16_t > mat_horizontal;
	matrix< 3, 4, uint16_t > test_mat_horizontal;
	//uint16_t data4[] = { 6, 5, 3, 6, 11, 3, 6, 7, 7, 6, 1, 3 };
	uint16_t data4[] = { 6, 6, 6, 6, 5, 11, 7, 1, 3, 3, 7, 3 };
	test_mat_horizontal.set(data4, data4 + 12);
	t3.get_horizontal_slice( 1, mat_horizontal );
	
	matrix< 3, 4, uint16_t > mat_horizontal_2;
	matrix< 3, 4, uint16_t > test_mat_horizontal_2;
	test_mat_lateral_2.fill(5);
	t3.set_horizontal_slice( 1, test_mat_horizontal_2 );
	t3.get_horizontal_slice( 1, mat_horizontal_2 );
	
	if (mat_horizontal == test_mat_horizontal &&  mat_horizontal_2 == test_mat_horizontal_2 )
	{	
		log( "get/set_horizontal_slice()", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "after get_horizontal_slice at i1 = 0: " << mat_horizontal << std::endl
		<< "after set_horizontal_slice at i1 = 0: " << mat_horizontal_2 << std::endl;
		log_error( error.str() );
	}
	
	
	//test tensor3 matrix multiplication
	t3.fill_increasing_values();
	matrix<5, 4, uint16_t> u3;
	u3.fill(1);
	tensor3<2, 3, 5, uint16_t> t3_jji;
	t3_jji.multiply_horizontal(t3, u3 );
	
	tensor3<6, 3, 5, uint16_t> t3_iji;
	matrix<6, 2, uint16_t> u1;
	u1.fill(2);
	t3_iji.multiply_lateral(t3_jji, u1 );

	tensor3<6, 7, 5, uint16_t> t3_iii;
	matrix<7, 3, uint16_t> u2;
	u2.fill(3);
	t3_iii.multiply_frontal(t3_iji, u2 );
	 
	
	t3.fill_increasing_values();
	tensor3<6, 7, 5, uint16_t> t3_reco;
	t3_reco.full_tensor3_matrix_multiplication( t3, u1, u2, u3 ); 
	
	tensor3<6, 7, 5, uint16_t> t3_iii_test;
	t3_iii_test.fill(1656);
	if ( t3_iii_test == t3_reco && t3_iii_test == t3_iii )
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
	
	//matricization
	//matrix< I3, I1*I2, T> m_horizontal;
	//matrix< I1, I2*I3, T> m_lateral;
	//matrix< I2, I1*I3, T> m_frontal;
	matrix< 4, 6, uint16_t> m_horizontal;
	matrix< 2, 12, uint16_t> m_lateral;
	matrix< 3, 8, uint16_t> m_frontal;
	
	t3.horizontal_matricization( m_horizontal);
	t3.lateral_matricization( m_lateral);
	t3.frontal_matricization(  m_frontal);
	
	matrix< 4, 6, uint16_t> m_horizontal_test;
	uint16_t data5[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
	m_horizontal_test.set(data5, data5 + 24);
	matrix< 2, 12, uint16_t> m_lateral_test;
	uint16_t data6[] = { 0, 6, 12, 18, 1, 7, 13, 19, 2, 8, 14, 20, 3, 9, 15, 21, 4, 10, 16, 22, 5, 11, 17, 23 };
	m_lateral_test.set(data6, data6 + 24);
	matrix< 3, 8, uint16_t> m_frontal_test;
	uint16_t data7[] = { 0, 3, 6, 9, 12, 15, 18, 21, 1, 4, 7, 10, 13, 16, 19, 22, 2, 5, 8, 11, 14, 17, 20, 23 };
	m_frontal_test.set(data7, data7 + 24);
	
	if ( m_horizontal_test == m_horizontal && m_lateral_test == m_lateral && m_frontal_test == m_frontal )
	{	
		log( "matricization along all modes", true  );
	} else
	{
		std::stringstream error;
		error 
		<< "matricization_horizontal should be: " << std::endl << m_horizontal_test << std::endl
		<< "matricization_horizontal is: " << std::endl << m_horizontal << std::endl
		<< "matricization_lateral should be: " << std::endl << m_lateral_test << std::endl
		<< "matricization_lateral is: " << std::endl << m_lateral << std::endl
		<< "matricization_frontal should be: " << std::endl << m_frontal_test << std::endl
		<< "matricization_frontal is: " << std::endl << m_frontal << std::endl
		<< std::endl;
		log_error( error.str() );
	}
	
	
	//compute frobenius norm of a tensor3
	float_t f_norm_check = 65.75712889109438;
	
	t3.fill_increasing_values();
	float_t f_norm = t3.frobenius_norm();
	
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
	
	//set diagonal values in a cubic tensor3, i.e., R=I1, R=I2, R=I3
	tensor3< 3, 3, 3, uint16_t >  t3_diag;
	tensor3< 3, 3, 3, uint16_t >  t3_diag_check;
	t3_diag_check.zero(); t3_diag_check.at(1, 1, 1) = 1; t3_diag_check.at(2, 2, 2) = 2;
	vector< 3, uint16_t > diag_values;
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
	
	//tensor3 type conversion
	tensor3< 2, 3, 4, float_t >  t3_type_a;
    tensor3< 2, 3, 4, uint16_t >  t3_type_b;
    tensor3< 2, 3, 4, uint16_t >  t3_type_b_check;
	
	t3_type_a.fill(2.4); 
	t3_type_b_check.fill(2);
	t3_type_b.convert_from_type( t3_type_a );
	
	if (t3_type_b_check == t3_type_b )
	{	
		log( "type conversion ", true  );
	} else
	{
		std::stringstream error;
		error << "type conversion - tensor3 type float_t: " << std::endl << t3_type_a << std::endl
		<< " tensor3 type uint16_t should be: " << std::endl << t3_type_b_check << std::endl
		<< " is: " << t3_type_b << std::endl;
		log_error( error.str() );
	}
	
	
	//export
	tensor3< 2, 3, 4, float_t >  t3_export_import;
	t3_export_import.fill_increasing_values();
	std::vector< float_t > export_data;
	t3_export_import.export_to( export_data );
	
	float_t export_data_check[] = { 0, 3, 1, 4, 2, 5, 6, 9, 7, 10, 8, 11, 12, 15, 13, 16, 14, 17, 18, 21, 19, 22, 20, 23 };
	float_t precision = 1.0e-1;
	ok = true;
	for (int i = 0; i < 24 && ok; ++i )
	{
		ok = (abs(export_data_check[i]) - abs(export_data[i])) < precision;
	}
	
	log( "export tensor3", ok  );
	
	//import tucker3 from vector
	std::vector< float_t > in_data = export_data;
	tensor3< 2, 3, 4, float_t >  t3_import_check;
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
	
	
	
	
	
	ok = true;
    return ok;
}


	
} // namespace vmml

