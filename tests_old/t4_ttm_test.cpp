#include "t4_ttm_test.hpp"

#include <vmmlib/t4_ttm.hpp>

namespace vmml
{

	bool
	t4_ttm_test::run()
	{
        bool global_ok = true;
		bool ok = false;


		tensor4< 2, 2, 2, 2, int> t4_input;
		int data_input[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
		t4_input.set( data_input, data_input + 16 );
        matrix< 2, 2, int> matrix_factor;
        int data_matrix_factor[] = { 6, 7, 8, 9 };
        matrix_factor.set(data_matrix_factor, data_matrix_factor+4);

        tensor4< 2, 2, 2, 2, int> t4_result;
        tensor4< 2, 2, 2, 2, int> t4_check;


        int data_check1[] = { 14, 27, 18, 35, 66, 79, 86, 103, 118, 131, 154, 171, 170, 183, 222, 239 };
        t4_check.set( data_check1, data_check1 + 16 );
        t4_ttm::mode1_multiply_fwd(t4_input, matrix_factor, t4_result);
        TEST( t4_result == t4_check );
        if ( ok )
		{
			log( "tensor4 matrix multiplication along mode 1", true  );
		} else
		{
			std::stringstream error;
			error
			<< "Result is: " << std::endl
            << t4_result << std::endl
            << "Result should be: " << std::endl
            << t4_check	<< std::endl;
			log_error( error.str() );
		}

        int data_check2[] = { 7, 9, 33, 43, 59, 77, 85, 111, 111, 145, 137, 179, 163, 213, 189, 247 };
        t4_check.set( data_check2, data_check2 + 16 );
        t4_ttm::mode2_multiply_fwd(t4_input, matrix_factor, t4_result);
        TEST( t4_result == t4_check );
        if ( ok )
		{
			log( "tensor4 matrix multiplication along mode 2", true  );
		} else
		{
			std::stringstream error;
			error
			<< "Result is: " << std::endl
            << t4_result << std::endl
            << "Result should be: " << std::endl
            << t4_check	<< std::endl;
			log_error( error.str() );
		}

        int data_check3[] = { 28, 41, 54, 67, 36, 53, 70, 87, 132, 145, 158, 171, 172, 189, 206, 223 };
        t4_check.set( data_check3, data_check3 + 16 );
        t4_ttm::mode3_multiply_fwd(t4_input, matrix_factor, t4_result);
        TEST( t4_result == t4_check );
        if ( ok )
		{
			log( "tensor4 matrix multiplication along mode 3", true  );
		} else
		{
			std::stringstream error;
			error
			<< "Result is: " << std::endl
            << t4_result << std::endl
            << "Result should be: " << std::endl
            << t4_check	<< std::endl;
			log_error( error.str() );
		}

        int data_check4[] = { 56, 69, 82, 95, 108, 121, 134, 147, 72, 89, 106, 123, 140, 157, 174, 191 };
        t4_check.set( data_check4, data_check4 + 16 );
        t4_ttm::mode4_multiply_fwd(t4_input, matrix_factor, t4_result);
        TEST( t4_result == t4_check );
        if ( ok )
		{
			log( "tensor4 matrix multiplication along mode 4", true  );
		} else
		{
			std::stringstream error;
			error
			<< "Result is: " << std::endl
            << t4_result << std::endl
            << "Result should be: " << std::endl
            << t4_check	<< std::endl;
			log_error( error.str() );
		}

		return global_ok;
	}

} // namespace vmml
