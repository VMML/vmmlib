#include "intersection_test.hpp"

#include <vmmlib/intersection.hpp>
#include <vmmlib/vector.hpp>
#include <sstream>

namespace vmml
{

bool intersection_test::run()
{
    bool ok = true;
    bool global_ok = true;

    const Vector< 4, float > sphere = Vector< 4, float >(10.f, 0.f, 0.f, 0.1f);
    const Vector< 3, float > rayOrigin = Vector< 3, float >(0.f, 0.f, 0.f);
    const Vector< 3, float > rayDir = Vector< 3, float >(1.f, 0.f, 0.f);

    Intersection< float > intersection( rayOrigin, rayDir);

    float t = 0, tmp = 9.9f;
    intersection.test_sphere(sphere, t);

    TEST(t == tmp);

    log( "intersection, raySphereIntersect( sphere, t)", ok );
    if ( ! ok )
    {
        std::stringstream error;
        error << "Intersection:\n"
              << "result should be " << tmp << ",\n"
              << "result is: " << t << std::endl;
        log_error( error.str() );
    }

    return global_ok;
}

} // namespace vmml
