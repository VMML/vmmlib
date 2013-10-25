#include "lapack_gaussian_elimination_test.hpp"


#include <vmmlib/lapack_gaussian_elimination.hpp>

#define TEST( x ) \
{ \
    ok = x; \
    global_ok &= ok; \
}

namespace vmml
{

bool
lapack_gaussian_elimination_test::run()
{
    bool global_ok = true;

    {
        bool ok = true;

        matrix< 3, 3, float >   A;
        vector< 3, float >      B;
        vector< 3, float >      X;

        A( 0, 0 )   = 2;
        A( 0, 1 )   = 1;
        A( 0, 2 )   = -1;
        B( 0 )      = 8;

        A( 1, 0 )   = -3;
        A( 1, 1 )   = -1;
        A( 1, 2 )   = 2;
        B( 1 )      = -11;

        A( 2, 0 )   = -2;
        A( 2, 1 )   = 1;
        A( 2, 2 )   = 2;
        B( 2 )      = -3;


        vmml::lapack::gaussian_elimination< 1, 3, float > ge;
        try
        {
            X   = B;
            ge.compute( A, X );
            float epsilon = 1e-6;

            //std::cout << " X " << X << std::endl;

            TEST( fabs( X( 0 ) - 2.0f ) < epsilon
                && fabs( X( 1 ) - 3.0f ) < epsilon
                &&  fabs( X( 2 ) - -1.0f ) < epsilon
                );
        }
        catch(...)
        {
            TEST(false);
            std::cout << ge.get_params() << std::endl;
        }

        log( "gaussian elimination using lapack xGESV", ok );
        if ( ! ok )
        {
            std::cout << A << std::endl;
            std::cout << B << std::endl;
            std::cout << X << std::endl;
        }

    }

    return global_ok;
}

} // namespace vmml

