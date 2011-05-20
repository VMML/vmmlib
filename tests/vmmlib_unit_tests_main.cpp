
#include "vector_test.hpp"
#include "matrix_test.hpp"
#include "quaternion_test.hpp"
#include "qr_decomposition_test.hpp"
#include "svd_test.hpp"
#include "tensor3_test.hpp"
#include "tensor3_iterator_test.hpp"
#include "tucker3_tensor_test.hpp"
#include "cp3_tensor_test.hpp"
#include "matrix_pseudoinverse_test.hpp"

#include <iostream>

//#define VMMLIB_USE_LAPACK

#ifdef VMMLIB_USE_LAPACK
#include "lapack_linear_least_squares_test.hpp"
#include "lapack_gaussian_elimination_test.hpp"
#include "lapack_svd_test.hpp"
#include "lapack_sym_eigs_test.hpp"
#endif

void
run_and_log( vmml::unit_test& test )
{
    test.run();
    std::cout << test << std::endl; 
}

int
main( int argc, const char* argv[] )
{

    vmml::vector_test vector_test_;
    run_and_log( vector_test_ );

    vmml::matrix_test matrix_test_;
    run_and_log( matrix_test_ );

    vmml::quaternion_test quaternion_test_;
    run_and_log( quaternion_test_ );
    
	vmml::qr_decomposition_test qr_test_;
	run_and_log( qr_test_ );

	vmml::svd_test svd_test_;
	run_and_log( svd_test_ );
    
#ifdef VMMLIB_USE_LAPACK
    vmml::lapack_svd_test lapack_svd_test_;
    run_and_log( lapack_svd_test_ );

    vmml::lapack_linear_least_squares_test lapack_llsq_test;
    run_and_log( lapack_llsq_test );

    vmml::lapack_gaussian_elimination_test lapack_ge_test;
    run_and_log( lapack_ge_test );
	
	vmml::lapack_sym_eigs_test lapack_sym_eigs_test_;
    run_and_log( lapack_sym_eigs_test_ );

#endif
	
	//problems on some unix machines
    vmml::matrix_pseudoinverse_test m_pinv;
    run_and_log( m_pinv );
	
    vmml::tensor3_test t3t;
    run_and_log( t3t );
	
    vmml::tensor3_iterator_test t3it;
    run_and_log( t3it );
	
    vmml::tucker3_tensor_test tt3t;
    run_and_log( tt3t );
	
	vmml::cp3_tensor_test cp3t;
    run_and_log( cp3t );
	
}
