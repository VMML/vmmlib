
#include "matrix_test.hpp"
#include "qr_decomposition_test.hpp"


#include <iostream>

void
run_and_log( vmml::unit_test& test )
{
    test.run();
    std::cout << test << std::endl; 
}

int
main( int argc, const char* argv[] )
{
    vmml::matrix_test matrix_test_;
    run_and_log( matrix_test_ );

	vmml::qr_decomposition_test qr_test_;
	run_and_log( qr_test_ );

}
