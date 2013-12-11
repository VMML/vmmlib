#include <vmmlib/tensor4.hpp>
#include "t4_hooi_test.hpp"
#include "vmmlib/t4_hooi.hpp"
#include "vmmlib/tucker4_tensor.hpp"

#include <sstream>

namespace vmml
{

	bool
	t4_hooi_test::run()
	{
        bool global_ok = true;
		bool ok = false;

		double precision = 0.001;

		tensor4< 2, 2, 2, 2, double> t4_data;
		double data[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
		t4_data.set( data, data + 16 );

		matrix< 2, 1, double > u1;
		matrix< 2, 1, double > u2;
		matrix< 2, 1, double > u3;
        matrix< 2, 1, double > u4;
		matrix< 2, 1, double > u1_check;
		matrix< 2, 1, double > u2_check;
		matrix< 2, 1, double > u3_check;
        matrix< 2, 1, double > u4_check;

		double data_u1[] = { 0.6327, 0.7744 };
		u1_check.set( data_u1, data_u1 + 2);

		double data_u2[] = { 0.6709, 0.7415 };
		u2_check.set( data_u2, data_u2 + 2);

		double data_u3[] = { 0.5482, 0.8363 };
		u3_check.set( data_u3, data_u3 + 2);

        double data_u4[] = { 0.3182, 0.9480 };
		u4_check.set( data_u4, data_u4 + 2);

		tensor4< 1, 1, 1, 1, double > core;
		tensor4< 1, 1, 1, 1, double > core_check;
		double data_core[] = { 34.9529 };
		core_check.set( data_core, data_core + 1);

		typedef t4_hooi< 1, 1, 1, 1, 2, 2, 2, 2, double > hooi_type;
		hooi_type::als( t4_data, u1, u2, u3, u4, core, hooi_type::init_hosvd() );

		TEST(u1.equals( u1_check, precision ) &&
                u2.equals( u2_check, precision ) &&
                u3.equals( u3_check, precision ) &&
                u3.equals( u3_check, precision ) &&
                core.equals( core_check, precision  ));
		if ( ok )
		{
			log( "HOOI rank-(1,1,1) approximation" , ok  );
		} else
		{
			std::stringstream error;
			error
			<< "HOOI rank-(1,1,1) approximation: " << std::setprecision(16) << std::endl
			<< "U1 should be: " << std::endl << u1_check << std::endl
			<< "U1 is: " << std::endl << u1 << std::endl
			<< "U2 should be: " << std::endl << u2_check << std::endl
			<< "U2 is: " << std::endl << u2 << std::endl
			<< "U3 should be: " << std::endl << u3_check << std::endl
			<< "U3 is: " << std::endl << u3 << std::endl
			<< "core should be: " << std::endl << core_check << std::endl
			<< "core is: " << std::endl << core << std::endl;


			log_error( error.str() );
		}

		return global_ok;
	}

} //end vmml namespace
