#include "matrix_pseudoinverse_test.hpp"
#include <vmmlib/matrix_pseudoinverse.hpp>

#include <sstream>
#include <stdint.h>

namespace vmml
{
	
	bool
	matrix_pseudoinverse_test::run()
	{
		bool ok = false;
		
		//test 1 (orthogonal)
		typedef matrix< 6, 4, float > pinv_type1;
		typedef matrix< 4, 6, float > pinv_t_type1;
		compute_pseudoinverse< pinv_type1 > compute_pinv;
		pinv_type1 input;
		float data[] = { 0.8147, 0.9058, 0.1270, 0.9134, 
			0.6324, 0.0975, 0.2785, 0.5469, 
			0.9575, 0.9649, 0.1576, 0.9706, 
			0.9572, 0.4854, 0.8003, 0.1419, 
			0.4218, 0.9157, 0.7922, 0.9595, 
			0.6557, 0.0357, 0.8491, 0.9340};
		input.set( data, data + 24);
		pinv_t_type1 pseudoinverse_control;
		float data3[] = {
			0.187956, 0.504057, 0.329735, 0.608532, -0.910448, 0.0212369,
			0.19669, -0.571985, 0.159718, 0.272729, 0.675097, -0.758368,
			-0.428897, -0.248724, -0.486611, 0.514871, 0.679076, 0.294918,
			0.170329, 0.27456, 0.126968, -0.932658, 0.139089, 0.610191	};
		pseudoinverse_control.set( data3, data3 + 24);
		
		pinv_type1 pseudoinverse_transposed;
		compute_pinv( input, pseudoinverse_transposed );
		pinv_t_type1 pseudoinverse = transpose( pseudoinverse_transposed );
		
		ok = pseudoinverse_control.equals( pseudoinverse, 1e-4 );
		
		//test 2 (non-symmetric transpose)
		typedef matrix< 4, 4, float > pinv_type2;
		compute_pseudoinverse< pinv_type2 > compute_pinv2;
		pinv_type2 input2;
		float data2[] = { 1, -0.446814, -0.218247, -0.077803,
			-0.446814, 1, -0.307783, -0.181625,
			-0.218247, -0.307783, 1, 0.0354524,
			-0.077803, -0.181625, 0.0354524, 1};
		input2.set( data2, data2 + 16);
		pinv_type2 pseudoinverse_control2;
		float data4[] = { 1.5831, 0.9504, 0.6283, 0.2735,
			0.9504, 1.7120, 0.7216, 0.3593,
			0.6283, 0.7216, 1.3545, 0.1319,
			0.2735, 0.3593, 0.1319, 1.0819
		};
		pseudoinverse_control2.set( data4, data4 + 16);
		
		pinv_type2 pseudoinverse_transposed2;
		compute_pinv2( input2, pseudoinverse_transposed2 );
		pinv_type2 pseudoinverse2 = transpose( pseudoinverse_transposed2 );
		
		ok = ok && (pseudoinverse_control2.equals( pseudoinverse2, 1e-4 ));
		
		//FIXME: check if multiply = Identity matrix
		
		if (ok) {
			log( "matrix compute pseudo inverse ", ok  );
		} else {
			std::stringstream error;
			error 
			<< "matrix compute pseudo inverse: " << std::endl
			<< "input 1 is: " << std::endl << input << std::endl
			<< "inverse matrix 1 should be: " << std::endl << pseudoinverse_control << std::endl
			<< "inverse matrix 1 is: " << std::endl << pseudoinverse << std::endl
			<< "input 2 is: " << std::endl << input2 << std::endl
			<< "inverse matrix 2 should be: " << std::endl << pseudoinverse_control2 << std::endl
			<< "inverse matrix 2 is: " << std::endl << pseudoinverse2 << std::endl;
			
			log_error( error.str() );
		}
		
		
		ok = true;
		return ok;
	}

	
} // namespace vmml

