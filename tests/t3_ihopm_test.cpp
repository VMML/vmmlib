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
		
		typedef tensor3< 4, 4, 4, double > cp_t3_type;
		
		cp_t3_type t3_cp_input;
		double data_in_cp[] = { 
			0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858, 
			0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
			0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
			0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071 };
		t3_cp_input.set(data_in_cp, data_in_cp + 64);

		typedef matrix< 4, 3, double > cp3_u_type;
		typedef vector< 3, double> cp3_lambda_type;

		cp3_lambda_type lambda;
		cp3_u_type u1;
		cp3_u_type u2;
		cp3_u_type u3;

		t3_ihopm< 1, 3, 4, 4, 4, double, double >::incremental_als( t3_cp_input, u1, u2, u3, lambda );
		
		//std::cout << "u1:\n" << u1 << std::endl << "u2:\n" << u2 << std::endl << "u3:\n" << u3 << std::endl << "lambda\n" << lambda << std::endl;


		//check test data
		cp3_lambda_type lambda_check;
		cp3_u_type u1_check;
		cp3_u_type u2_check;
		cp3_u_type u3_check;
		lambda_check.at(0) = 2.7868; lambda_check.at(1) = 1.592; 
		//lambda_check.at(2) = 2.3321; lambda_check.at(3) = 2.1248;       
		//lambda_check.at(4) = 0.39103; lambda_check.at(5) = 0.37917;       
		
		double data_u1_check[] = {
			0.527015410617948, 0.056351846996132, -0.270070568153916, -0.313500982905647, 0.643680959587782,  0.718173845512749,
			0.480746418639447,  0.030284666943727, -0.415759124642206, -0.549335893818258,  0.055046229660312, -0.112051357316600,
			0.371207457664265,  0.088601553082487, -0.152884114581588, -0.185456056239643,  0.739465596886156,  0.559308690697938,
			0.594426329590955,  0.994010600084126,  0.854887528276306,  0.752032752401201, -0.189302313472758, -0.398553145081928 };
		u1_check.set( data_u1_check, data_u1_check + 24 );
		
		double data_u2_check[] = {
			0.504403339491653,  0.114653823530311,  0.828382311442367,  0.781391033358520,  0.286987251907120,  0.429484828957028,
			0.487454156012892,  0.587181630899377,  0.196061458997553,  0.152685111239459, -0.276096237356745, -0.156053868043593,
			0.554624618255953,  0.795085965591298, -0.254655993603039, -0.262100760689507, -0.218430969484626, -0.162045823851050,
			0.447612834623702,  0.099551697142648, -0.458795134354323, -0.545360890639107,  0.890896793422366,  0.874603408944202 };
		u2_check.set( data_u2_check, data_u2_check + 24 );
		
		double data_u3_check[] = { 
			0.338579272657363,  0.943081819028568,  0.050178095908785, -0.136825354681717, -0.813169091174282,  0.893225962438723,
			0.451425279509148,  0.136957766308855,  0.841276694398122, -0.695287017909927,  0.577211828569151, -0.375455071169286,
			0.664949661681733, -0.261520371343423, -0.431693562496846,  0.614846308347820, -0.073109135396691,  0.071281938453576,
			0.489306898149053, -0.153122004420833, -0.321521928104942,  0.346148526146221, -0.015413904119256,  0.236853868047333 };
		
		u3_check.set( data_u3_check, data_u3_check + 24 );
		
		ok = u1.equals( u1_check, precision );
		ok = u2.equals( u2_check, precision ) && ok;
		ok = u3.equals( u3_check, precision ) && ok;
		ok = lambda.equals( lambda_check, precision ) && ok;
		
		//if( ok)
		{	
			log( "todo: incremental rank-R approximation (CP-ALS/HOPM)", ok  );
		} 
		
		//FIX: HOPM problem different number of iterations (convergence criteria?)
#if 0
		else
		{
			std::stringstream error;
			error 
			<< "incremental rank-R approximation (CP-ALS/HOPM)" << std::setprecision(16) << std::endl
			<< " lambda should be: " << lambda_check << "lambda is: " << lambda	<< std::endl
			<< " u1 should be: " << std::endl << u1_check << std::endl
			<< " u1 is: " << std::endl << u1 << std::endl
			<< " u2 should be: " << std::endl << u2_check << std::endl
			<< " u2 is: " << std::endl << u2 << std::endl
			<< " u3 should be: " << std::endl << u3_check << std::endl
			<< " u3 is: " << std::endl << u3 << std::endl;
			
			log_error( error.str() );
		}
#endif
		
		return ok;
	}
	
} //end vmml namespace