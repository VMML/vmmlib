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
		
		
		compute_pseudoinverse< matrix< 6, 4, float > > compute_pinv;
		matrix< 6, 4, float > input;
		float data[] = { 0.8147, 0.9058, 0.1270, 0.9134, 
			0.6324, 0.0975, 0.2785, 0.5469, 
			0.9575, 0.9649, 0.1576, 0.9706, 
			0.9572, 0.4854, 0.8003, 0.1419, 
			0.4218, 0.9157, 0.7922, 0.9595, 
			0.6557, 0.0357, 0.8491, 0.9340};
		input.set( data, data + 24);
		matrix< 4, 6, float > pseudoinverse_control;
		float data3[] = {
			0.8147, 0.6324, 0.9575, 0.9572, 0.4218, 0.6557,
			0.9058, 0.0975, 0.9649, 0.4854, 0.9157, 0.0357,
			0.1270, 0.2785, 0.1576, 0.8003, 0.7922, 0.8491,
			0.9134, 0.5469, 0.9706, 0.1419, 0.9595,  0.9340
		};
		
		pseudoinverse_control.set( data3, data3 + 24);
		
		matrix< 6, 4, float > pseudoinverse_transposed;
		compute_pinv( input, pseudoinverse_transposed );

		matrix< 4, 6, float > pseudoinverse = transpose( pseudoinverse_transposed );
		
		if ( pseudoinverse_control.equals( pseudoinverse, 1e-4 ))
		{	
			log( "matrix compute pseudo inverse ", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "matrix compute pseudo inverse: " << std::endl
			<< "input is: " << std::endl << input << std::endl
			<< "inverse matrix should be: " << std::endl << pseudoinverse_control << std::endl
			<< "inverse matrix is: " << std::endl << pseudoinverse << std::endl;
			
			log_error( error.str() );
		}
		
		
		ok = true;
		return ok;
	}

	
} // namespace vmml