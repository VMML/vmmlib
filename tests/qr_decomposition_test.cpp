#include "qr_decomposition_test.hpp"

#include <vmmlib/qr_decomposition.hpp>
#include <vmmlib/matrix.hpp>

#include <sstream>

namespace vmml
{
bool
qr_decomposition_test::run()
{
    bool ok = true;
    
    // tests copyFrom1DimCArray function with row_by_row data
	ok = true;
	{
            matrix< 3, 3, double > A;
            matrix< 3, 3, double > Q;
            matrix< 3, 3, double > Q_correct;
            matrix< 3, 3, double > R;
            matrix< 3, 3, double > R_correct;

            double data[ 3 * 3 ] = { 12, -51, 4, 6, 167, -68, -4, 24, -41 };///, 2, 5, -23 };
            A = data;
                      
            qr_decompose( A, Q, R );

            double correct_solution_R[ 3 * 3 ] = { 14, 21, -14, 0, 175, -70, 0, 0, 35 };
            R_correct.copyFrom1DimCArray( correct_solution_R );

            double correct_solution_Q[ 3 * 3 ] = 
                { 0.857143, 0.428571, -0.285714, -0.394286, 0.902857, 0.171429, -0.331429, 0.034286, -0.942857 };
            Q_correct.copyFrom1DimCArray( correct_solution_Q, false );

            ok = ( R == R_correct && Q == Q_correct );
            log( "qr_decomposition of a 3x3 matrix using stabilized Gram-Schmidt (maximum precision)", ok, true  );

            if ( ! ok )
            {
                ok = R.isEqualTo( R_correct, 1e-9 ) && Q.isEqualTo( Q_correct, 1e-9 );
                
                log( "qr_decomposition a 3x3 matrix using stabilized Gram-Schmidt with precision tolerance 1e-9", ok );
            }

		if ( ! ok )
		{
			std::stringstream error;
			error << " A " << A << std::endl;
			error << " Q " << Q << std::endl;
			error << " Q correct" << Q_correct << std::endl;
			error << " R " << R << std::endl;
			error << " R correct " << R_correct << std::endl;
			log_error( error.str() );
		}
	
	}
    
    return ok;
}


} // namespace vmml

