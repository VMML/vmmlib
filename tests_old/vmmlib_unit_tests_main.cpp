#include <iostream>

#include "unit_test_globals.hpp"

#include "matrix_test.hpp"
#include "lowpass_filter_test.hpp"
#include "frustum_test.hpp"
#include "intersection_test.hpp"
#include "quaternion_test.hpp"
#include "qr_decomposition_test.hpp"
#include "util_test.hpp"

void run_and_log( vmml::unit_test& test )
{
    test.run();
    std::cout << test << std::endl;
}

int main( int, const char** )
{
    vmml::matrix_test matrix_test_;
    run_and_log( matrix_test_ );

    vmml::lowpass_filter_test lowpass_filter_test_;
    run_and_log( lowpass_filter_test_ );

    vmml::frustum_test frustum_test_;
    run_and_log( frustum_test_ );

    vmml::intersection_test intersection_test_;
    run_and_log( intersection_test_ );

    vmml::quaternion_test quaternion_test_;
    run_and_log( quaternion_test_ );

	vmml::qr_decomposition_test qr_test_;
	run_and_log( qr_test_ );

    vmml::util_test util_test_;
    run_and_log( util_test_ );

    std::cout << vmml::unit_test_globals::get_instance() << std::endl;

}
