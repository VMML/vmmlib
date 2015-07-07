#include "qr_decomposition_test.hpp"

#include <vmmlib/qr_decomposition.hpp>
#include <vmmlib/matrix.hpp>

#include <sstream>

namespace vmml
{
bool
qr_decomposition_test::run()
{
    bool global_ok = true;
    bool ok = true;

    // tests qr decomposition using modified gram-schmidt
    {
        Matrix< 3, 3, double > A, Q, R;
        double Adata[] = { 12, -51, 4, 6, 167, -68, -4, 24, -41 };
        A.set( Adata, Adata + 9 );
        qr_decompose_gram_schmidt( A, Q, R );

        double Qcorrect[] = {
            6./7, -69./175, 58./175,
            3./7, 158./175, -6./175,
            -2./7, 6./35, 33./35 };

        double Rcorrect[] = { 14, 21, -14, 0, 175, -70, 0, 0, -35 };
        Matrix< 3, 3, double > Qc, Rc;
        Qc.set( Qcorrect, Qcorrect + 9 );
        Rc.set( Rcorrect, Rcorrect + 9 );

        for( size_t index = 0; ok && index < 3; ++index )
        {
            Vector< 3, double > q = Q.get_column( index );
            Vector< 3, double > qc = Qc.get_column( index );
            TEST(q == qc);
            if ( ! ok )
            {
                q *= -1.0;
                TEST(q == qc);
            }
        }
        for( size_t index = 0; ok && index < 3; ++index )
        {
            Vector< 3, double > r = R.get_row( index );
            Vector< 3, double > rc = Rc.get_row( index );
            TEST(r == rc);
            if ( ! ok )
            {
                r *= -1.0;
                TEST(r == rc);
            }
        }

        log( "QR decomposition using modified gram-schmidt, maximal precision", ok, true );

        if ( ! ok )
        {
            ok = true;
            double tolerance = 1e-9;
            for( size_t index = 0; ok && index < 3; ++index )
            {
                Vector< 3, double > q = Q.get_column( index );
                Vector< 3, double > qc = Qc.get_column( index );
                TEST(q.equals( qc, tolerance ));
                if ( ! ok )
                {
                    q *= -1.0;
                    TEST(q.equals( qc, tolerance ));
                }
            }
            for( size_t index = 0; ok && index < 3; ++index )
            {
                Vector< 3, double > r = R.get_row( index );
                Vector< 3, double > rc = Rc.get_row( index );
                TEST(r.equals( rc, tolerance ));
                if ( ! ok )
                {
                    r *= -1.0;
                    TEST(r.equals( rc ));
                }
            }
            log( "QR decomposition using modified gram-schmidt, tolerance 1e-9", ok );
            if ( ! ok )
            {
                std::stringstream error;
                error << " A " << A << std::endl;
                error << " Q " << Q << std::endl;
                error << " Qc " << Qc << std::endl;
                error << " diff " << Q - Qc << std::endl;
                error << " R " << R << std::endl;
                error << " Rc " << Rc << std::endl;
                error << " diff " << R - Rc << std::endl;
                std::cout << error.str() << std::endl;
            }
        }
    }
    return global_ok;
}


} // namespace vmml
