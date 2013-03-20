#include "intersection_test.hpp"

#include <vmmlib/intersection.hpp>
#include <vmmlib/vector.hpp>
#include <sstream>

namespace vmml
{
bool intersection_test::run()
{
    const vector< 4, float > sphere = vector< 4, float >(10.f, 0.f, 0.f, 0.1f);
    const vector< 3, float > rayOrigin = vector< 3, float >(0.f, 0.f, 0.f);
    const vector< 3, float > rayDir = vector< 3, float >(1.f, 0.f, 0.f);

    intersection< float > intersection( rayOrigin, rayDir);

    float t = 0, tmp = 9.9f;
    intersection.test_sphere(sphere, t);

    bool ok = t == tmp;

    log( "intersection, raySphereIntersect( sphere, t)", ok );
    if ( ! ok )
    {
        std::stringstream error;
        error << "Intersection:\n"
              << "result should be " << tmp << ",\n"
              << "result is: " << t << std::endl;
        log_error( error.str() );
    }

    return ok;
}

} // namespace vmml
