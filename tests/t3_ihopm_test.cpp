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
#define W 2
#define L 3
		typedef matrix< 4, W*L, double > cp3_u_type;
		typedef vector< W*L, double> cp3_lambda_type;

		cp3_lambda_type lambda;
		cp3_u_type u1;
		cp3_u_type u2;
		cp3_u_type u3;

		t3_ihopm< W, L, 4, 4, 4, double, double >::incremental_als( t3_cp_input, u1, u2, u3, lambda );
		
		//std::cout << "u1:\n" << u1 << std::endl << "u2:\n" << u2 << std::endl << "u3:\n" << u3 << std::endl << "lambda\n" << lambda << std::endl;


		//check test data
		cp3_lambda_type lambda_check;
		cp3_u_type u1_check;
		cp3_u_type u2_check;
		cp3_u_type u3_check;
		lambda_check.at(0) = 2.9053; lambda_check.at(1) = 1.4844;
		lambda_check.at(2) = 0.70494; lambda_check.at(3) = 0.53852;
		lambda_check.at(4) = 0.34382; lambda_check.at(5) = 0.34367;

		double data_u1_check[] = {
			0.506085561468077, 0.001956847654422, -0.209569603142863, -0.394921956232513, -0.022486848229877, 0.184594118898338,
			0.460495692379029, -0.028557624219977, -0.180356605196994, -0.734631384709598, -0.486345040443633, -0.556857121401639,
			0.359317314562872, 0.054450501667292, -0.110310141873280, -0.269443896857940, -0.745176978684274, -0.455743051491161,
			0.634596083528118, 0.998106094417997, 0.954664207457395, 0.481407689522252, -0.455712753531810, -0.669427687379266,
		};
		u1_check.set( data_u1_check, data_u1_check + 24 );
		
		double data_u2_check[] = {
			-0.496342174057493, -0.038914583436338, 0.866036783907366, -0.562676911251632, 0.417141865256690, 0.349220775883442,
			-0.499456882130393, -0.573265837144616, 0.302670396624642, -0.123572086450233, -0.189594213931581, -0.304319214625875,
			-0.569895344707819, -0.811943783843057, -0.177969059686749, 0.182161561359157, -0.405487465007436, -0.208074746650915,
			-0.423564121739347, -0.102952547503582, -0.355946532109365, 0.796832352855890, 0.790965621260291, 0.861475226055291,
		};
		u2_check.set( data_u2_check, data_u2_check + 24 );
		
		double data_u3_check[] = { 
			-0.388210367985676, -0.916926601731041, 0.004409907552016, 0.259294452470892, -0.156403724994960, 0.498783074179373,
			-0.466062051276868, -0.053679545922064, 0.922185157079488, 0.262424594049043, -0.686783735606445, 0.485694863120284,
			-0.641760774047993, 0.342843706058446, -0.200373280225951, -0.840021412092139, 0.013805150959998, 0.142358862753370,
			-0.469278151463868, 0.197033770211290, -0.330765229957981, -0.397823763220280, -0.709700918078263, 0.703597824787933,
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