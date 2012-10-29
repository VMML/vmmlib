#include "t3_ihopm_test.hpp"
#include "vmmlib/t3_ihopm.hpp"
#include <sstream>

namespace vmml
{
	
	bool
	t3_ihopm_test::run()
	{
		bool ok = false;
		
		double precision = 0.001;
		
		tensor3< 4, 4, 4, double > t3_cp_input;
		double data_in_cp[] = { 
			0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858, 
			0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
			0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
			0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071 };
		t3_cp_input.set(data_in_cp, data_in_cp + 64);
#define W 2
#define L 3
		typedef matrix< 4, W*L, double > cp3_u_type;
		typedef vector< W*L, double> cp3_lambda_type;

		cp3_lambda_type lambda;
		cp3_u_type u1;
		cp3_u_type u2;
		cp3_u_type u3;

		t3_ihopm< W*L, L, 4, 4, 4, double, double >::incremental_als( t3_cp_input, u1, u2, u3, lambda, 20 );
		
		//std::cout << "u1:\n" << u1 << std::endl << "u2:\n" << u2 << std::endl << "u3:\n" << u3 << std::endl << "lambda\n" << lambda << std::endl;


		//check test data
		cp3_lambda_type lambda_check;
		cp3_u_type u1_check;
		cp3_u_type u2_check;
		cp3_u_type u3_check;
                
		lambda_check.at(0) = 2.8700; lambda_check.at(1) = 1.5134;
		lambda_check.at(2) = 0.7419; lambda_check.at(3) = 0.5627;
		lambda_check.at(4) = 0.3273; lambda_check.at(5) = 0.3248;

		double data_u1_check[] = {
                0.5130,    0.0172,   -0.2130,    -0.3809,   0.0379,    0.2447,
                0.4671,   -0.0112,   -0.2093,    -0.7284,    -0.4237,   -0.5201,
                0.3632,    0.0641,   -0.1125,    -0.2532,    -0.8201,   -0.5425,
                0.6218,    0.9977,    0.9477,   0.5102,    -0.3827,   -0.6127,
		};
		u1_check.set( data_u1_check, data_u1_check + 24 );
		
		double data_u2_check[] = {
                -0.4996,    -0.0615,    0.8666,   -0.6119,    0.4305,    0.3545,
                -0.4958,    -0.5782,    0.2842,   -0.1137,   -0.1430,   -0.2926,
                -0.5648,    -0.8075,   -0.1975,    0.2110,   -0.3390,   -0.1387,
                -0.4309,    -0.0994,   -0.3596,    0.7538,    0.8242,    0.8772,
		};
		u2_check.set( data_u2_check, data_u2_check + 24 );
		
		double data_u3_check[] = { 
                -0.3748,    -0.9273,   -0.0100,   0.2564,    -0.1903,    0.5207,
                -0.4621,    -0.0808,    0.9288,   0.3269,    -0.6436,    0.4204,
                -0.6483,   0.3165,   -0.1998,    -0.8255,   0.0154,    0.1232,
                -0.4750,   0.1828,   -0.3119,    -0.3821,    -0.7412,    0.7327,
		};
		
		u3_check.set( data_u3_check, data_u3_check + 24 );
		
		ok = u1.equals( u1_check, precision );
		ok = u2.equals( u2_check, precision ) && ok;
		ok = u3.equals( u3_check, precision ) && ok;
		ok = lambda.equals( lambda_check, precision ) && ok;
		
		if( ok)
		{	
			log( "incremental rank-R approximation", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "incremental rank-R approximation" << std::setprecision(16) << std::endl
			<< "lambda should be: " << lambda_check << std::endl
			<< "lambda is: " << lambda	<< std::endl
			<< "u1 should be: " << std::endl << u1_check << std::endl
			<< "u1 is: " << std::endl << u1 << std::endl
			<< "u2 should be: " << std::endl << u2_check << std::endl
			<< "u2 is: " << std::endl << u2 << std::endl
			<< "u3 should be: " << std::endl << u3_check << std::endl
			<< "u3 is: " << std::endl << u3 << std::endl;
			
			log_error( error.str() );
		}
		
		return ok;
	}
	
} //end vmml namespace
