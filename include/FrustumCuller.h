/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Stefan Eilemann
*
* @license revised BSD license, check LICENSE
*/ 

#ifndef __VMML__FRUSTUM_CULLER__H__
#define __VMML__FRUSTUM_CULLER__H__

#include "Vector4.h" // member

// - declaration -

namespace vmml
{
    template< typename T > class Matrix4;

enum Visibility
{
    VISIBILITY_NONE     = 0,
    VISIBILITY_PARTIAL,
    VISIBILITY_FULL,
};

/** Helper class for OpenGL view frustum culling. */
template< class T > 
class FrustumCuller
{
public:
    // contructors
    FrustumCuller() {}// warning: components NOT initialised ( for performance )
    ~FrustumCuller(){}

    void setup( const Matrix4<T>& projection );
    Visibility testSphere( const Vector4<T>& sphere );

private:
    Vector4<T> _leftPlane;
    Vector4<T> _rightPlane;
    Vector4<T> _bottomPlane;
    Vector4<T> _topPlane;
    Vector4<T> _nearPlane;
    Vector4<T> _farPlane;
};

#ifndef VMMLIB_DISABLE_TYPEDEFS

    typedef FrustumCuller< float >  FrustumCullerf;
    typedef FrustumCuller< double > FrustumCullerd;

#endif
}

// - implementation - //
#include "Matrix4.h"

namespace vmml
{

/** 
 * Setup the culler by extracting the frustum planes from the projection
 * matrix. The projection matrix should contain the viewing transformation.
 */
template < class T > 
void FrustumCuller< T >::setup( const Matrix4<T>& projection )
{
    // See http://www2.ravensoft.com/users/ggribb/plane%20extraction.pdf pp.5
    _leftPlane   = projection.getRow(3) + projection.getRow(0);
    _rightPlane  = projection.getRow(3) - projection.getRow(0);
    _bottomPlane = projection.getRow(3) + projection.getRow(1);
    _topPlane    = projection.getRow(3) - projection.getRow(1);
    _nearPlane   = projection.getRow(3) + projection.getRow(2);
    _farPlane    = projection.getRow(3) - projection.getRow(2);

    _leftPlane.normalize();
    _rightPlane.normalize();
    _bottomPlane.normalize(); 
    _topPlane.normalize();
    _nearPlane.normalize();
    _farPlane.normalize();
}

template < class T > 
Visibility FrustumCuller< T >::testSphere( const Vector4<T>& sphere )
{
    Visibility visibility = VISIBILITY_FULL;

    // see http://www.flipcode.com/articles/article_frustumculling.shtml
    // distance = plane.normal . sphere.center + plane.distance
    // Test all planes:
    // - if sphere behind plane: not visible
    // - if sphere intersect one plane: partially visible
    // - else: fully visible

    T distance = _leftPlane.normal.x * sphere.center.x +
                 _leftPlane.normal.y * sphere.center.y +
                 _leftPlane.normal.z * sphere.center.z + _leftPlane.distance;
    if( distance <= -sphere.radius )
        return VISIBILITY_NONE;
    if( distance < sphere.radius )
        visibility = VISIBILITY_PARTIAL;

    distance = _rightPlane.normal.x * sphere.center.x +
               _rightPlane.normal.y * sphere.center.y +
               _rightPlane.normal.z * sphere.center.z + _rightPlane.distance;
    if( distance <= -sphere.radius )
        return VISIBILITY_NONE;
    if( distance < sphere.radius )
        visibility = VISIBILITY_PARTIAL;

    distance = _bottomPlane.normal.x * sphere.center.x +
               _bottomPlane.normal.y * sphere.center.y +
               _bottomPlane.normal.z * sphere.center.z + _bottomPlane.distance;
    if( distance <= -sphere.radius )
        return VISIBILITY_NONE;
    if( distance < sphere.radius )
        visibility = VISIBILITY_PARTIAL;

    distance = _topPlane.normal.x * sphere.center.x +
               _topPlane.normal.y * sphere.center.y +
               _topPlane.normal.z * sphere.center.z + _topPlane.distance;
    if( distance <= -sphere.radius )
        return VISIBILITY_NONE;
    if( distance < sphere.radius )
        visibility = VISIBILITY_PARTIAL;

    distance = _nearPlane.normal.x * sphere.center.x +
               _nearPlane.normal.y * sphere.center.y +
               _nearPlane.normal.z * sphere.center.z + _nearPlane.distance;
    if( distance <= -sphere.radius )
        return VISIBILITY_NONE;
    if( distance < sphere.radius )
        visibility = VISIBILITY_PARTIAL;

    distance = _farPlane.normal.x * sphere.center.x +
               _farPlane.normal.y * sphere.center.y +
               _farPlane.normal.z * sphere.center.z + _farPlane.distance;
    if( distance <= -sphere.radius )
        return VISIBILITY_NONE;
    if( distance < sphere.radius )
        visibility = VISIBILITY_PARTIAL;

    return visibility;
}	
}
#endif
