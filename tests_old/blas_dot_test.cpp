#include "blas_dot_test.hpp"

#include <vmmlib/blas_dot.hpp>

namespace vmml
{

	bool
	blas_dot_test::run()
	{
        bool global_ok = true;
		bool ok = false;

		vector< 4, double > A;
		vector< 4, double > B;

		double AData[] = { 1, 2, 3, 4};
		A = AData;
		double BData[] = { 5, 6, 7, 8};
		B = BData;

		double dot_prod = 0;
		double dot_prod_check = 70;

		blas_dot< 4, double > blas_dot1;
		blas_dot1.compute( A, B, dot_prod );

		TEST(dot_prod == dot_prod_check);

		log( "dot product", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "dot product of \n" << A << "\n"
			<< "and \n" << B << "\n"
			<< "should be\n" << dot_prod_check << "\n"
            << "is\n" << dot_prod << "\n"
            << std::endl;
			log_error( ss.str() );
		}


		return global_ok;
	}



} // namespace vmml
