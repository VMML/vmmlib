#include "matrix_pseudoinverse_test.hpp"
#include <vmmlib/matrix_pseudoinverse.hpp>

#include <sstream>
#include <stdint.h>

namespace vmml
{

	bool
	matrix_pseudoinverse_test::run()
	{
        bool global_ok = true;
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

		TEST(pseudoinverse_control.equals( pseudoinverse, 1e-4 ));

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

		pinv_type2 pseudoinverse_transposed2( pinv_type2::IDENTITY );
		compute_pinv2( input2, pseudoinverse_transposed2 );
		pinv_type2 pseudoinverse2 = transpose( pseudoinverse_transposed2 );

		if (ok) TEST(pseudoinverse_control2.equals( pseudoinverse2, 1e-4 ));

		//FIXME: check if multiply = Identity matrix

		if (ok) {
			log( "matrix compute pseudo inverse (ROWS > COLS)", ok  );
		} else {
			std::stringstream error;
			error
			<< "matrix compute pseudo inverse (ROWS > COLS): " << std::endl
			<< "input 1 is: " << std::endl << input << std::endl
			<< "inverse matrix 1 should be: " << std::endl << pseudoinverse_control << std::endl
			<< "inverse matrix 1 is: " << std::endl << pseudoinverse << std::endl
			<< "input 2 is: " << std::endl << input2 << std::endl
			<< "inverse matrix 2 should be: " << std::endl << pseudoinverse_control2 << std::endl
			<< "inverse matrix 2 is: " << std::endl << pseudoinverse2 << std::endl;
			log_error( error.str() );
		}

		//cols > rows
		typedef matrix< 4, 6, float > pinv_type3;

		pinv_type3 input3;
		input3.set( data, data + 24);

		compute_pseudoinverse< pinv_type3 > compute_pinv3;

		pinv_type3 pseudoinverse_transposed3;
		compute_pinv3( input3, pseudoinverse_transposed3 );

		pinv_type3 pseudoinverse3_control;
		float data6[] = {
			0.5110213756561279, 0.1466298252344131, -0.2957937121391296, 0.5260885953903198, 0.08560904860496521, -0.4743495881557465,
			-0.6723353862762451, 0.1120621934533119, 0.3330267071723938, 0.5148389935493469, -0.1580927520990372, 0.3453871607780457,
			1.576845765113831, -0.9367456436157227, 0.3191587626934052, -0.09735386073589325, -0.6063382625579834, -0.04424500837922096,
			-0.9502310156822205, 0.8770782351493835, -0.1227491199970245, -0.5944689512252808, 0.7954459190368652, 0.3613623976707458
		};
		pseudoinverse3_control.set( data6, data6 + 24);

		TEST(pseudoinverse3_control.equals( pseudoinverse_transposed3, 1e-6 ));

		if (ok) {
			log( "matrix compute pseudo inverse (COLS > ROWS) ", ok  );
		} else {
			std::stringstream error;
			error
			<< "matrix compute pseudo inverse (COLS > ROWS): " << std::endl
			<< "input 3 is: " << std::endl << input3 << std::endl
			<< "pinv matrix 3 transposed should be: " << std::endl << pseudoinverse3_control << std::endl
			<< "pinv matrix 3 transposed is: " << std::endl << pseudoinverse_transposed3 << std::endl;
			log_error( error.str() );
		}



		return global_ok;
	}


} // namespace vmml
