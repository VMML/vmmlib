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

#ifndef _Vector4_H_
#define _Vector4_H_

#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <assert.h>

// - declaration -
namespace vmml
{
template< typename T > class Vector3;
template< typename T > class Matrix4;

template< typename T > 
class Vector4
{
public:
    union
    {
        struct
        {
            T  x, y, z, w;
        };
        struct
        {
            T  r, g, b, a;
        };
        T xyzw[4];
        T rgba[4];
    };

    // contructors
    Vector4(); // warning: components NOT initialised ( for performance )
    Vector4( const T aa ); 
    Vector4( const T xx, const T yy, const T zz, const T ww ); 
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector< double >
    //the pointer 'values must be a valid 4 component c array of the resp. type
    Vector4( const float* aa );
    Vector4( const double* aa ); 
    Vector4( const Vector3< T >& xxyyzz, const T aa );   
    ~Vector4();

    void set( T xx, T yy, T zz, T ww );
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector< double >
    //the pointer 'values must be a valid 4 component c array of the resp. type
    void set( const float* aa );
    void set( const double* aa );

    const Vector4& operator=( T aa ); 
    const Vector4& operator=( const Vector4& aa ); 

    T& operator[]( size_t position);
    const T& operator[]( size_t position) const;

    T length() const; 
    T lengthSquared() const;

    T normalise();

    // vector/scalar operations
    Vector4 operator+( const T aa ) const;
    Vector4 operator-( const T aa ) const; 
    Vector4 operator*( const T aa ) const;
    Vector4 operator/( T aa ) const;
    Vector4 operator-() const;
     
    const Vector4& operator+=( T aa );
    const Vector4& operator-=( T aa );
    const Vector4& operator*=( T aa );
    const Vector4& operator/=( T aa ); 

    // vector/vector operations
    Vector4 operator+( const Vector4 &aa ) const; 
    Vector4 operator-( const Vector4 &aa ) const;
    Vector4 operator*( const Vector4 &aa ) const; 
    Vector4 operator/( const Vector4 &aa ) const; 

    const Vector4& operator+=( const Vector4 &aa ); 
    const Vector4& operator-=( const Vector4 &aa ); 
    const Vector4& operator*=( const Vector4 &aa ); 
    const Vector4& operator/=( const Vector4 &aa ); 

    // vector'/matrix operations
    Vector4 operator* (const Matrix4<T>& m ) const;

    T dot( const Vector4 &aa ) const;

    // other
    bool operator==( const Vector4 &aa ) const;
    bool operator!=(const Vector4 &aa ) const;
    bool isAkin( const Vector4& a, 
                 const T& delta = std::numeric_limits<T>::epsilon() );

    void invert(); 

    T min();
    T max();

    friend std::ostream& operator << ( std::ostream& os, const Vector4& v )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << "[" << std::setw(10) << v.x << " " << std::setw(10) << v.y 
           << " " << std::setw(10) << v.z << " " << std::setw(10) << v.w << "]";
        os.precision( prec );
        os.setf( flags );
        return os;
    };        
};

#ifndef VMMLIB_DISABLE_TYPEDEFS
    typedef Vector4<float>  Vector4f;
    typedef Vector4<double> Vector4d;
#endif
}
    
// - implementation - //
#include "Vector3.h"
#include "Matrix4.h"

namespace vmml
{
template < typename T > 
Vector4< T >::Vector4() 
{} 

template < typename T > 
Vector4< T >::Vector4( const T  aa )
    : x( aa )
    , y( aa )
    , z( aa ) 
    , w( aa )
{} 

template < typename T > 
Vector4< T >::Vector4( const T xx, const T yy, const T zz, 
                       const T ww )
    : x( xx )
    , y( yy )
    , z( zz ) 
    , w( ww )
{} 

template < typename T > 
Vector4< T >::Vector4( const float* values )
{
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
    w = static_cast< T > ( values[3] );
}

template < typename T > 
Vector4< T >::Vector4( const double* values )
{
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
    w = static_cast< T > ( values[3] );
}

template < typename T > 
Vector4< T >::Vector4( const Vector3< T >& v3, const T aa )
    : x ( v3.x )
    , y ( v3.y )
    , z ( v3.z )
    , w ( aa )
{} 

template < typename T > 
Vector4< T >::~Vector4()
{}

template < typename T > 
void Vector4< T >::set( T xx, T yy, T zz, T ww )
{ 
    x = xx; 
    y = yy; 
    z = zz; 
    w = ww;
}

template < typename T > 
void Vector4< T >::set( const float* values )
{ 
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
    w = static_cast< T > ( values[3] );
}

template < typename T > 
void Vector4< T >::set( const double* values )
{ 
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
    w = static_cast< T > ( values[3] );
}

template < typename T > 
const Vector4< T >& Vector4< T >::operator=( T aa )
{ 
    x = y = z = w = aa; 
    return *this; 
} 

template < typename T > 
const Vector4< T >& Vector4< T >::operator=( const Vector4& aa ) 
{ 
    x = aa.x; 
    y = aa.y; 
    z = aa.z; 
    w = aa.w;
    return *this;
} 


template < typename T > 
T& Vector4< T >::operator[]( size_t index ) 
{ 
    assert( index < 4 && "Vector4::operator[] Invalid component index!" ); 
    return xyzw[ index ]; 
}
         
template < typename T > 
const T& Vector4< T >::operator[]( size_t index ) const
{ 
    assert( index < 4 && "Vector4::operator[] Invalid component index!" ); 
    return xyzw[ index ]; 
} 
	
template < typename T > 
T  Vector4< T >::length() const 
{ 
    const T l = lengthSquared();
    return ( l <= 0 ) ? 0 : sqrt( l ); 
} 

template < typename T > 
T  Vector4< T >::lengthSquared() const 
{ 
    return x * x + y * y + z * z + w * w; 
} 

template < typename T > 
T Vector4< T >::normalise()
{ 
    T l = length(); 
    if ( l == 0 ) 
        return 0; 
    l = 1.0f / l; 
    x *= l; 
    y *= l; 
    z *= l; 
    w *= l;
    return l; 
} 

template < typename T > 
Vector4< T > Vector4< T >::operator+( const T  aa ) const 
{ 
    return Vector4( x + aa, y + aa, z + aa, w + aa ); 
} 

template < typename T > 
Vector4< T > Vector4< T >::operator-( const T  aa ) const 
{ 
    return Vector4( x - aa, y - aa, z - aa, w - aa ); 
}
 
template < typename T > 
Vector4< T > Vector4< T >::operator*( const T  aa ) const 
{ 
    return Vector4( x * aa, y * aa, z * aa, w * aa ); 
}

template < typename T > 
Vector4< T > Vector4< T >::operator/( T  aa ) const 
{ 
    assert( aa != 0.0f ); 
    aa = 1.0f / aa; 
    return Vector4( x * aa, y * aa, z * aa, w * aa ); 
}

template < typename T > 
Vector4< T > Vector4< T >::operator-() const 
{ 
    return Vector4( -x, -y, -z, -w );
}


template < typename T > 
const Vector4< T >& Vector4< T >::operator+=( T  aa ) 
{ 
    x += aa; 
    y += aa; 
    z += aa; 
    w += aa;
    return *this; 
} 

template < typename T > 
const Vector4< T >& Vector4< T >::operator-=( T  aa ) 
{ 
    x -= aa; 
    y -= aa; 
    z -= aa; 
    w -= aa;
    return *this; 
} 

template < typename T > 
const Vector4< T >& Vector4< T >::operator*=( T  aa ) 
{ 
    x *= aa; 
    y *= aa; 
    z *= aa;
    w *= aa;
    return *this; 
}
 
template < typename T > 
const Vector4< T >& Vector4< T >::operator/=( T  aa ) 
{ 
    aa = 1.0f / aa; 
    x *= aa; 
    y *= aa; 
    z *= aa; 
    w *= aa;
    return *this; 
} 

// vector/vector operations
template < typename T > 
Vector4< T > Vector4< T >::operator+( const Vector4 &aa ) const 
{ 
    return Vector4( x + aa.x, y + aa.y, z + aa.z, w + aa.w ); 
}
 
template < typename T > 
Vector4< T > Vector4< T >::operator-( const Vector4 &aa ) const 
{ 
    return Vector4( x - aa.x, y - aa.y, z - aa.z, w - aa.w ); 
}

template < typename T > 
Vector4< T > Vector4< T >::operator*( const Vector4 &aa ) const 
{ 
    return Vector4( x * aa.x, y * aa.y, z * aa.z, w * aa.w ); 
} 

template < typename T > 
Vector4< T > Vector4< T >::operator/( const Vector4 &aa ) const 
{ 
    return Vector4( x / aa.x, y / aa.y, z / aa.z, w / aa.w ); 
} 

template < typename T > 
const Vector4< T >& Vector4< T >::operator+=( const Vector4 &aa ) 
{ 
    x += aa.x; 
    y += aa.y; 
    z += aa.z;
    w += aa.w; 
    return *this; 
} 

template < typename T > 
const Vector4< T >& Vector4< T >::operator-=( const Vector4 &aa ) 
{ 
    x -= aa.x; 
    y -= aa.y; 
    z -= aa.z;
    w -= aa.w; 
    return *this; 
}

template < typename T > 
const Vector4< T >& Vector4< T >::operator*=( const Vector4 &aa ) 
{ 
    x *= aa.x; 
    y *= aa.y; 
    z *= aa.z; 
    w *= aa.w;
    return *this; 
}

template < typename T > 
const Vector4< T >& Vector4< T >::operator/=( const Vector4 &aa ) 
{ 
    x /= aa.x; 
    y /= aa.y; 
    z /= aa.z;
    w /= aa.w; 
    return *this; 
}

// vector'/matrix operations    
template< typename T > 
Vector4< T > Vector4< T >::operator * (const Matrix4<T>& m ) const
{
	return Vector4< T >( x * m.m00 + y * m.m01 + z * m.m02 + w * m.m03,
                         x * m.m10 + y * m.m11 + z * m.m12 + w * m.m13,
                         x * m.m20 + y * m.m21 + z * m.m22 + w * m.m23,
                         x * m.m30 + y * m.m31 + z * m.m32 + w * m.m33 );
}

template < typename T > 
T Vector4< T >::dot( const Vector4 &aa ) const 
{ 
    return x * aa.x + y * aa.y + z * aa.z + w * aa.w; 
}

template < typename T > 
bool Vector4< T >::operator==( const Vector4 &aa ) const 
{ 
    return ( x == aa.x && y == aa.y && z == aa.z && w = aa.w ); 
}

template < typename T > 
bool Vector4< T >::operator!=(const Vector4 &aa ) const 
{ 
    return ( x != aa.x || y != aa.y || z != aa.z || w != aa.w ); 
}
template < typename T > 
bool Vector4< T >::isAkin( const Vector4& a, const T& delta )
{
    if( fabs( x-a.x ) > delta || fabs( y-a.y ) > delta || 
        fabs( z-a.z ) > delta || fabs( w-a.w ) > delta )

        return false;
    return true;
}
template< > 
inline bool Vector4< float >::isAkin( const Vector4& a, const float& delta )
{
    if( fabsf( x-a.x ) > delta || fabsf( y-a.y ) > delta || 
        fabsf( z-a.z ) > delta || fabsf( w-a.w ) > delta )

        return false;
    return true;
}

template < typename T > 
void Vector4< T >::invert() 
{	
    x = -x; 
    y = -y; 
    z = -z; 
    w = -w;
}

template < typename T > 
T Vector4< T >::max() 
{ 
    T m = std::max( x, y ); 
    m = std::max( m, z); 
    m = std::max( m, w);
    return m; 
}

template < typename T > 
T Vector4< T >::min() 
{ 
    T m = std::min( x, y ); 
    m = std::min( m, z); 
    m = std::min( m, w);
    return m; 
} 
}	
#endif
