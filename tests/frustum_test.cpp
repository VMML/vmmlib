#include "frustum_test.hpp"

#include <vmmlib/frustum.hpp>
#include <vmmlib/frustum_culler.hpp>

#define TEST( x, info )                                               \
    {                                                                 \
        if( !(x) )                                                    \
        {                                                             \
            std::ostringstream os;                                      \
            os << #x << ": " << info;                                   \
            log( os.str(), false );                                     \
            return false;                                               \
        }                                                               \
        else                                                            \
        {                                                               \
            log( #x, true );                                            \
        }                                                               \
    }

namespace vmml
{

bool frustum_test::run() 
{
    const frustum< float > f( -1.f, 1., -1.f, 1., 1.f, 100.f );
    const matrix< 4, 4, float > mvp = f.compute_matrix();

    frustum_culler< float > fc;
    fc.setup( mvp );

    const vector< 4, float > sphereIn( 0.f, 0.f, -10.f, 1.f );
    const vector< 4, float > sphereOut( 0.f, 0.f, 0.f, .5f );
    const vector< 4, float > sphereBorder( 0.f, 0.f, -1.f, 1.f );

    TEST( fc.test_sphere( sphereIn ) == VISIBILITY_FULL,
          fc.test_sphere( sphereIn ));
    TEST( fc.test_sphere( sphereOut ) == VISIBILITY_NONE,
          fc.test_sphere( sphereOut ));
    TEST( fc.test_sphere( sphereBorder ) == VISIBILITY_PARTIAL,
          fc.test_sphere( sphereBorder ));
    return true;
}


} // namespace vmml

