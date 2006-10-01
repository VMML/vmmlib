/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Jonas Boesch
* @author Stefan Eilemann
* @author Renato Pajarola
* @author David H. Eberly ( Wild Magic )
* @author Andrew Willmott ( VL )
*
* @license BSD license, check LICENSE
*
* parts of the source code of VMMLib were inspired by David Eberly's 
* Wild Magic and Andrew Willmott's VL.
* 
*/ 



#ifndef _Vector3_H_
#define _Vector3_H_

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

// - declaration -

namespace vmml
{

template< class Real > 
class Vector3
{
public:
    union
    {
        struct
        {
            Real  x, y, z;
        };
        struct
        {
            Real  r, g, b;
        };
        Real xyz[3];
        Real rgb[3];
    };

    // contructors
    Vector3(); // warning: components NOT initialised ( for performance )
    Vector3( const Real  a ); 
    Vector3( const Real  i, const Real  j, const Real  k ); 
    
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector3< double >
    //the pointer 'values must be a valid 3 component c array of the resp. type
    Vector3( const float* values );
    Vector3( const double* values );
    
    ~Vector3();

    void set( Real xx, Real yy, Real zz );
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector3< double >
    // the pointer 'values' must be a valid 3 component c array of the resp. type
    void set( const float* values );
    void set( const double* values );

    const Vector3& operator=( Real a ); 
    const Vector3& operator=( const Vector3& a ); 

    Real& operator[]( size_t position);
    const Real& operator[]( size_t position) const;

    // vector/scalar operations
    Vector3 operator+( const Real  a ) const;
    Vector3 operator-( const Real  a ) const; 
    Vector3 operator*( const Real  a ) const;
    Vector3 operator/( Real  a ) const;
     
    const Vector3& operator+=( Real  a );
    const Vector3& operator-=( Real  a );
    const Vector3& operator*=( Real  a );
    const Vector3& operator/=( Real  a ); 

    // vector/vector operations
    Vector3 operator+( const Vector3& a ) const; 
    Vector3 operator-( const Vector3& a ) const;
    Vector3 operator*( const Vector3& a ) const; 
    Vector3 operator/( const Vector3& a ) const; 
    Vector3 operator-() const;

    const Vector3& operator+=( const Vector3& a ); 
    const Vector3& operator-=( const Vector3& a ); 
    const Vector3& operator*=( const Vector3& a ); 
    const Vector3& operator/=( const Vector3& a ); 

    bool operator==( const Vector3 &a ) const;
    bool operator!=(const Vector3 &a ) const;

    Real length() const;
    Real lengthSquared() const;

    Real normalise();
    Real normalise( const float* source );
    void scale( Real scale_factor );

    // result = vec1.cross( vec2 ) => vec1 x vec2
    Vector3 cross( const Vector3& a ) const;
    // result.cross( vec1, vec2 ) => vec1 x vec2
    void cross( const Vector3 &a, const Vector3 &b);
    Real dot( const Vector3& a) const;
    static Real dot( const Vector3& a, const Vector3& b);

    void invert(); 

    // *this is the result
    void normal( const Vector3& aa, const Vector3& bb, const Vector3& cc );
    //returns the normal of *this and the two argument vectors
    Vector3 normal( const Vector3& aa, const Vector3& bb );

    Vector3 rotate( Real theta, Real  rx, Real  ry, Real  rz ); 

    Real min();
    Real max();

    friend std::ostream& operator<<( std::ostream& o, const Vector3& v )
    {
        o << "Vector3(" << v.x << "," << v.y << "," << v.z << ") ";
        return o;
    };        
};
    
// - implementation - //
       
template < class Real > 
Vector3< Real >::Vector3() 
{} 

template < class Real > 
Vector3< Real >::Vector3( const Real  a )
    : x(a)
    , y(a)
    , z(a) 
{} 

template < class Real > 
Vector3< Real >::Vector3( const Real  i, const Real  j, const Real  k )
    : x(i)
    , y(j)
    , z(k) 
{} 

template < class Real > 
Vector3< Real >::Vector3( const float* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
}

template < class Real > 
Vector3< Real >::Vector3( const double* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
}


template < class Real > 
Vector3< Real >::~Vector3()
{}

template < class Real > 
void Vector3< Real >::set( Real xx, Real yy, Real zz )
{ 
    x = xx; 
    y = yy; 
    z = zz; 
}

template < class Real > 
void Vector3< Real >::set( const float* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
}

template < class Real > 
void Vector3< Real >::set( const double* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
}

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator=( Real a )
{ 
    x = y = z = a; 
    return *this; 
} 

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator=( const Vector3& a ) 
{ 
    x = a.x; 
    y = a.y; 
    z = a.z; 
    return *this;
} 


template < class Real > 
Real& Vector3< Real >::operator[]( size_t position) 
{ 
    assert( position < 3 && "Vector3::operator[] Invalid component index!" ); 
    return * (&x+position); 
}
         
template < class Real > 
const Real& Vector3< Real >::operator[]( size_t position) const
{ 
    assert( position < 3 && "Vector3::operator[] Invalid component index!" ); 
    return * (&x+position); 
} 
	
template < class Real > 
Real  Vector3< Real >::length() const 
{ 
    Real l = lengthSquared();
    return ( l <= 0 ) ? 0 : sqrt( l ); 
} 

template < class Real > 
Real  Vector3< Real >::lengthSquared() const 
{ 
    return x * x + y * y + z * z; 
} 

template < class Real > 
Real  Vector3< Real >::normalise()
{ 
    Real l = length(); 
    if ( l == 0 ) 
        return 0; 
    l = 1.0f / l; 
    x *= l; 
    y *= l; 
    z *= l; 
    return l; 
} 

template < class Real > 
Real  Vector3< Real >::normalise( const float* source )
{
    Vector3< Real > a ( source );
    Real l = a.length();
    if ( l == 0 ) 
        return 0;
     
}


template < class Real >
void Vector3< Real >::scale( Real scale_factor )
{
    operator*=( scale_factor );
}

template < class Real > 
Vector3< Real > Vector3< Real >::operator+( const Real  a ) const 
{ 
    return Vector3( x + a, y + a, z + a ); 
} 

template < class Real > 
Vector3< Real > Vector3< Real >::operator-( const Real  a ) const 
{ 
    return Vector3( x - a, y - a, z - a ); 
}
 
template < class Real > 
Vector3< Real > Vector3< Real >::operator*( const Real  a ) const 
{ 
    return Vector3( x * a, y * a, z * a ); 
}

template < class Real > 
Vector3< Real > Vector3< Real >::operator/( Real  a ) const 
{ 
    assert( a != 0.0f ); 
    a = 1.0f / a; 
    return Vector3( x * a, y * a, z * a ); 
}

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator+=( Real  a ) 
{ 
    x += a; 
    y += a; 
    z += a; 
    return *this; 
} 

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator-=( Real  a ) 
{ 
    x -= a; 
    y -= a; 
    z -= a; 
    return *this; 
} 

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator*=( Real  a ) 
{ 
    x *= a; 
    y *= a; 
    z *= a; 
    return *this; 
}
 
template < class Real > 
const Vector3< Real >& Vector3< Real >::operator/=( Real  a ) 
{ 
    a = 1.0f / a; 
    x *= a; 
    y *= a; 
    z *= a; 
    return *this; 
} 

// vector/vector operations
template < class Real > 
Vector3< Real > Vector3< Real >::operator+( const Vector3 &a ) const 
{ 
    return Vector3( x + a.x, y + a.y, z + a.z ); 
}
 
template < class Real > 
Vector3< Real > Vector3< Real >::operator-( const Vector3 &a ) const 
{ 
    return Vector3( x - a.x, y - a.y, z - a.z ); 
}

template < class Real > 
Vector3< Real > Vector3< Real >::operator*( const Vector3 &a ) const 
{ 
    return Vector3( x * a.x, y * a.y, z * a.z ); 
} 

template < class Real > 
Vector3< Real > Vector3< Real >::operator/( const Vector3 &a ) const 
{ 
    return Vector3( x / a.x, y / a.y, z / a.z ); 
} 

template < class Real > 
Vector3< Real > Vector3< Real >::operator-() const 
{ 
    return Vector3( -x, -y, -z );
}

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator+=( const Vector3 &a ) 
{ 
    x += a.x; 
    y += a.y; 
    z += a.z; 
    return *this; 
} 

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator-=( const Vector3 &a ) 
{ 
    x -= a.x; 
    y -= a.y; 
    z -= a.z; 
    return *this; 
}

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator*=( const Vector3 &a ) 
{ 
    x *= a.x; 
    y *= a.y; 
    z *= a.z; 
    return *this; 
}

template < class Real > 
const Vector3< Real >& Vector3< Real >::operator/=( const Vector3 &a ) 
{ 
    x /= a.x; 
    y /= a.y; 
    z /= a.z; 
    return *this; 
}

// result = vec1.cross( vec2 ) => vec1 x vec2
template < class Real > 
Vector3< Real > Vector3< Real >::cross( const Vector3& a ) const
{ 
    return Vector3( y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x ); 
}

// result.cross( vec1, vec2 ) => vec1 x vec2
template < class Real > 
void Vector3< Real >::cross( const Vector3 &a, const Vector3 &b)
{ 
    x = a.y * b.z - a.z * b.y; 
    y = a.z * b.x - a.x * b.z; 
    z = a.x * b.y - a.y * b.x; 
}

template < class Real > 
Real Vector3< Real >::dot( const Vector3& a) const 
{ 
    return x * a.x + y * a.y + z * a.z; 
}

template < class Real > 
Real Vector3< Real >::dot( const Vector3& a, const Vector3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


template < class Real > 
bool Vector3< Real >::operator==( const Vector3 &a ) const 
{ 
    return ( x == a.x && y == a.y && z == a.z ); 
}

template < class Real > 
bool Vector3< Real >::operator!=(const Vector3 &a ) const 
{ 
    return ( x != a.x || y != a.y || z != a.z ); 
}

template < class Real > 
void Vector3< Real >::invert() 
{	
    x = -x; 
    y = -y; 
    z = -z; 
}

// *this is the result
template< class Real >
void Vector3< Real >::normal( const Vector3< Real >& aa, 
                              const Vector3< Real >& bb,
                              const Vector3< Real >& cc )
{
    Vector3< Real > u,v;

    // right hand system, CCW triangle
    u = bb - aa;
    v = cc - aa;
    cross( u, v );
    normalise();
}

//returns the normal of *this and the two argument vectors
template< class Real >
Vector3< Real > Vector3< Real >::normal( const Vector3< Real >& aa, 
                                         const Vector3< Real >& bb )
{
    Vector3< Real > tmp;
    tmp.normal( *this, aa, bb);
    return tmp;
}

template < class Real > 
Vector3< Real > Vector3< Real >::rotate( Real theta, Real rx, Real ry, 
                                         Real rz )      
{   
    Vector3 q( 0, 0, 0 ), r( rx, ry, rz );
    r.normalise();
    const Real costheta = ( Real ) cos( theta );
    const Real sintheta = ( Real ) sin( theta );
    q.x += ( costheta + ( 1.0f - costheta ) * r.x * r.x ) * x;
    q.x += ( ( 1 - costheta ) * r.x * r.y - r.z * sintheta ) * y;
    q.x += ( ( 1 - costheta ) * r.x * r.z + r.y * sintheta ) * z;
    q.y += ( ( 1 - costheta ) * r.x * r.y + r.z * sintheta ) * x;
    q.y += ( costheta + ( 1 - costheta ) * r.y * r.y ) * y;
    q.y += ( ( 1 - costheta ) * r.y * r.z - r.x * sintheta ) * z;
    q.z += ( ( 1 - costheta ) * r.x * r.z - r.y * sintheta ) * x;
    q.z += ( ( 1 - costheta ) * r.y * r.z + r.x * sintheta ) * y;
    q.z += ( costheta + ( 1 - costheta ) * r.z * r.z ) * z;
    return q; 
} 

#if 0
template<>
inline Vector3< float > Vector3< float >::rotate( float theta, float rx, 
                                                  float ry, float rz )      
{   
    Vector3 q( 0, 0, 0 ), r( rx, ry, rz );
    r.normalise();
    const float costheta = cosf( theta );
    const float sintheta = sinf( theta );
    q.x += ( costheta + ( 1.0f - costheta ) * r.x * r.x ) * x;
    q.x += ( ( 1 - costheta ) * r.x * r.y - r.z * sintheta ) * y;
    q.x += ( ( 1 - costheta ) * r.x * r.z + r.y * sintheta ) * z;
    q.y += ( ( 1 - costheta ) * r.x * r.y + r.z * sintheta ) * x;
    q.y += ( costheta + ( 1 - costheta ) * r.y * r.y ) * y;
    q.y += ( ( 1 - costheta ) * r.y * r.z - r.x * sintheta ) * z;
    q.z += ( ( 1 - costheta ) * r.x * r.z - r.y * sintheta ) * x;
    q.z += ( ( 1 - costheta ) * r.y * r.z + r.x * sintheta ) * y;
    q.z += ( costheta + ( 1 - costheta ) * r.z * r.z ) * z;
    return q; 
} 
#endif

template < class Real > 
Real Vector3< Real >::max() 
{ 
    Real m = std::max( x,y ); 
    m = std::max( m, z); 
    return m; 
}

template < class Real > 
Real Vector3< Real >::min() 
{ 
    Real m = std::min( x,y ); 
    m = std::min( m, z); 
    return m; 
} 
}	
#endif
