#ifndef _vmml_Sphere_H_
#define _vmml_Sphere_H_

#include "Vector3.h"

namespace vmml
{

template< typename T >
class Sphere
{

public:
    Sphere();
    Sphere( T cx, T cy, T cz, T radius );
    Sphere( const Vector3< T >& center, T radius );
    ~Sphere();
    
    inline T getRadius() const { return _radius; };
    void setRadius( T radius ) { _radius = radius; };

    inline const Vector3< T >& getCenter() const { return _center; };
    void setCenter( const Vector3< T >& center ) { _center = center; };
    
protected:
    Vector3< T > _center;
    T _radius;

};

template< typename T >
Sphere< T >::Sphere()
    : _center( .0, .0, .0 )
    , _radius( 1. )
{};

template< typename T >
Sphere< T >::Sphere( T cx, T cy, T cz, T radius )
    : _center( cx, cy, cz )
    , _radius( radius)
{};

template< typename T >
Sphere< T >::Sphere( const Vector3< T >& center, T radius )
    : _center( center )
    , _radius( radius )
{};

template< typename T >
Sphere< T >::~Sphere()
{};

};

#endif
