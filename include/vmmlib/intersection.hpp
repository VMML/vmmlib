/*
 * VMMLib - Intersection classes
 *
 * @author Jafet Villafranca
 *
 * Implements ray-sphere intersection using the three-dimensional coordinates
 * of the ray origin and direction, the sphere center and its radius.
 *
 */

#ifndef __VMML__INTERSECTION__HPP__
#define __VMML__INTERSECTION__HPP__

#include <vmmlib/vector.hpp>

namespace vmml
{
template< typename T > class intersection
{
public:
    typedef vector< 2, T >    vec2;
    typedef vector< 3, T >    vec3;
    typedef vector< 4, T >    vec4;

    /**
      Constructors
     */
    intersection() {}
    intersection( const vec3& rayOrigin, const vec3& rayDir );
    ~intersection() {}

    /**
      Sets ray parameters

      @param[in]    origin      Ray origin
      @param[in]    direction   Ray direction
     */
    void setRay( const vec3& origin, const vec3& direction );

    /**
      Ray Sphere Intersection - Optimized solution
      "Real-time Rendering 3rd Edition"

      @param[in]    center      Sphere center
      @param[in]    radius      Sphere radius
      @param[out]   t           Intersection distance

      @return Whether the ray intersects the sphere
     */
    bool raySphereIntersect( const vec4 &sphere, T &t ) const;

private:
    vec3 _origin;
    vec3 _direction;

}; // class intersection


template< typename T >
intersection< T >::intersection( const vec3& rayOrigin, const vec3& rayDir )
{
    setRay( rayOrigin, rayDir );
}

template< typename T >
void intersection< T >::setRay( const vec3& origin, const vec3& direction )
{
    vec3 direction_unit = direction;
    direction_unit.normalize();

    _origin = origin;
    _direction = direction_unit;
}

template< typename T >
bool
intersection< T >::raySphereIntersect( const vec4 &sphere, T &t ) const
{
    const vec3 center = vec3(sphere.x(), sphere.y(), sphere.z());
    const T radius = sphere.w();

    const vec3 centerVec = center - _origin;
    const T vecProjection = centerVec.dot(_direction);

    const T sqDistance = centerVec.squared_length();
    const T sqRadius = radius * radius;

    /** Sphere behind the ray origin && ray origin outside the sphere */
    if(vecProjection < 0 && sqDistance > sqRadius)
        return false;

    /** Squared distance from sphere center to the projection */
    const T sqCenterToProj = sqDistance - vecProjection * vecProjection;

    if(sqCenterToProj > sqRadius)
        return false;

    /** Distance from the sphere center to the surface along the ray direction*/
    const T distSurface = sqrt(sqRadius - sqCenterToProj);

    if(sqDistance > sqRadius)
        t = vecProjection - distSurface;
    else
        t = vecProjection + distSurface;

    return true;
}


} // namespace vmml

#endif
