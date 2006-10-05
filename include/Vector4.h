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
#include <iostream>
#include <algorithm>
#include <assert.h>

#include "Vector3.h"

// - declaration -
namespace vmml
{

template< class Real > 
class Vector4
{
public:
    union
    {
        struct
        {
            Real  x, y, z, w;
        };
        struct
        {
            Real  r, g, b, a;
        };
        Real xyzw[4];
        Real rgba[4];
    };

    // contructors
    Vector4(); // warning: components NOT initialised ( for performance )
    Vector4( const Real aa ); 
    Vector4( const Real xx, const Real yy, const Real zz, const Real ww ); 
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector< double >
    //the pointer 'values must be a valid 4 component c array of the resp. type
    Vector4( const float* aa );
    Vector4( const double* aa ); 
    Vector4( const Vector3< Real >& xxyyzz, const Real aa );   
    ~Vector4();

    void set( Real xx, Real yy, Real zz, Real ww );
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector< double >
    //the pointer 'values must be a valid 4 component c array of the resp. type
    void set( const float* aa );
    void set( const double* aa );

    const Vector4& operator=( Real aa ); 
    const Vector4& operator=( const Vector4& aa ); 

    Real& operator[]( size_t position);
    const Real& operator[]( size_t position) const;

    Real length() const; 
    Real lengthSquared() const;

    Real normalise();

    // vector/scalar operations
    Vector4 operator+( const Real aa ) const;
    Vector4 operator-( const Real aa ) const; 
    Vector4 operator*( const Real aa ) const;
    Vector4 operator/( Real aa ) const;
    Vector4 operator-() const;
     
    const Vector4& operator+=( Real aa );
    const Vector4& operator-=( Real aa );
    const Vector4& operator*=( Real aa );
    const Vector4& operator/=( Real aa ); 

    // vector/vector operations
    Vector4 operator+( const Vector4 &aa ) const; 
    Vector4 operator-( const Vector4 &aa ) const;
    Vector4 operator*( const Vector4 &aa ) const; 
    Vector4 operator/( const Vector4 &aa ) const; 

    const Vector4& operator+=( const Vector4 &aa ); 
    const Vector4& operator-=( const Vector4 &aa ); 
    const Vector4& operator*=( const Vector4 &aa ); 
    const Vector4& operator/=( const Vector4 &aa ); 

    Real dot( const Vector4 &aa ) const;

    // other
    bool operator==( const Vector4 &aa ) const;
    bool operator!=(const Vector4 &aa ) const;
    void invert(); 

    Real min();
    Real max();

    friend std::ostream& operator << ( std::ostream& os, const Vector4& v )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << "[" << std::setw(7) << v.x << " " << std::setw(7) << v.y 
           << " " << std::setw(7) << v.z << " " << std::setw(7) << v.w << "]";
        os.precision( prec );
        os.setf( flags );
        return os;
    };        
};


    
// - implementation - //
       
template < class Real > 
Vector4< Real >::Vector4() 
{} 

template < class Real > 
Vector4< Real >::Vector4( const Real  aa )
    : x( aa )
    , y( aa )
    , z( aa ) 
    , w( aa )
{} 

template < class Real > 
Vector4< Real >::Vector4( const Real xx, const Real yy, const Real zz, 
                          const Real ww )
    : x( xx )
    , y( yy )
    , z( zz ) 
    , w( ww )
{} 

template < class Real > 
Vector4< Real >::Vector4( const float* values )
{
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
    w = static_cast< Real > ( values[3] );
}

template < class Real > 
Vector4< Real >::Vector4( const double* values )
{
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
    w = static_cast< Real > ( values[3] );
}

template < class Real > 
Vector4< Real >::Vector4( const Vector3< Real >& v3, const Real aa )
    : x ( v3.x )
    , y ( v3.y )
    , z ( v3.z )
    , w ( aa )
{} 

template < class Real > 
Vector4< Real >::~Vector4()
{}

template < class Real > 
void Vector4< Real >::set( Real xx, Real yy, Real zz, Real ww )
{ 
    x = xx; 
    y = yy; 
    z = zz; 
    w = ww;
}

template < class Real > 
void Vector4< Real >::set( const float* values )
{ 
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
    w = static_cast< Real > ( values[3] );
}

template < class Real > 
void Vector4< Real >::set( const double* values )
{ 
    assert( values && "Vector4: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
    z = static_cast< Real > ( values[2] );
    w = static_cast< Real > ( values[3] );
}

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator=( Real aa )
{ 
    x = y = z = w = aa; 
    return *this; 
} 

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator=( const Vector4& aa ) 
{ 
    x = aa.x; 
    y = aa.y; 
    z = aa.z; 
    w = aa.w;
    return *this;
} 


template < class Real > 
Real& Vector4< Real >::operator[]( size_t index ) 
{ 
    assert( index < 4 && "Vector4::operator[] Invalid component index!" ); 
    return xyzw[ index ]; 
}
         
template < class Real > 
const Real& Vector4< Real >::operator[]( size_t index ) const
{ 
    assert( index < 4 && "Vector4::operator[] Invalid component index!" ); 
    return xyzw[ index ]; 
} 
	
template < class Real > 
Real  Vector4< Real >::length() const 
{ 
    const Real l = lengthSquared();
    return ( l <= 0 ) ? 0 : sqrt( l ); 
} 

template < class Real > 
Real  Vector4< Real >::lengthSquared() const 
{ 
    return x * x + y * y + z * z + w * w; 
} 

template < class Real > 
Real Vector4< Real >::normalise()
{ 
    Real l = length(); 
    if ( l == 0 ) 
        return 0; 
    l = 1.0f / l; 
    x *= l; 
    y *= l; 
    z *= l; 
    w *= l;
    return l; 
} 

template < class Real > 
Vector4< Real > Vector4< Real >::operator+( const Real  aa ) const 
{ 
    return Vector4( x + aa, y + aa, z + aa, w + aa ); 
} 

template < class Real > 
Vector4< Real > Vector4< Real >::operator-( const Real  aa ) const 
{ 
    return Vector4( x - aa, y - aa, z - aa, w - aa ); 
}
 
template < class Real > 
Vector4< Real > Vector4< Real >::operator*( const Real  aa ) const 
{ 
    return Vector4( x * aa, y * aa, z * aa, w * aa ); 
}

template < class Real > 
Vector4< Real > Vector4< Real >::operator/( Real  aa ) const 
{ 
    assert( aa != 0.0f ); 
    aa = 1.0f / aa; 
    return Vector4( x * aa, y * aa, z * aa, w * aa ); 
}

template < class Real > 
Vector4< Real > Vector4< Real >::operator-() const 
{ 
    return Vector4( -x, -y, -z, -w );
}


template < class Real > 
const Vector4< Real >& Vector4< Real >::operator+=( Real  aa ) 
{ 
    x += aa; 
    y += aa; 
    z += aa; 
    w += aa;
    return *this; 
} 

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator-=( Real  aa ) 
{ 
    x -= aa; 
    y -= aa; 
    z -= aa; 
    w -= aa;
    return *this; 
} 

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator*=( Real  aa ) 
{ 
    x *= aa; 
    y *= aa; 
    z *= aa;
    w *= aa;
    return *this; 
}
 
template < class Real > 
const Vector4< Real >& Vector4< Real >::operator/=( Real  aa ) 
{ 
    aa = 1.0f / aa; 
    x *= aa; 
    y *= aa; 
    z *= aa; 
    w *= aa;
    return *this; 
} 

// vector/vector operations
template < class Real > 
Vector4< Real > Vector4< Real >::operator+( const Vector4 &aa ) const 
{ 
    return Vector4( x + aa.x, y + aa.y, z + aa.z, w + aa.w ); 
}
 
template < class Real > 
Vector4< Real > Vector4< Real >::operator-( const Vector4 &aa ) const 
{ 
    return Vector4( x - aa.x, y - aa.y, z - aa.z, w - aa.w ); 
}

template < class Real > 
Vector4< Real > Vector4< Real >::operator*( const Vector4 &aa ) const 
{ 
    return Vector4( x * aa.x, y * aa.y, z * aa.z, w * aa.w ); 
} 

template < class Real > 
Vector4< Real > Vector4< Real >::operator/( const Vector4 &aa ) const 
{ 
    return Vector4( x / aa.x, y / aa.y, z / aa.z, w / aa.w ); 
} 

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator+=( const Vector4 &aa ) 
{ 
    x += aa.x; 
    y += aa.y; 
    z += aa.z;
    w += aa.w; 
    return *this; 
} 

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator-=( const Vector4 &aa ) 
{ 
    x -= aa.x; 
    y -= aa.y; 
    z -= aa.z;
    w -= aa.w; 
    return *this; 
}

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator*=( const Vector4 &aa ) 
{ 
    x *= aa.x; 
    y *= aa.y; 
    z *= aa.z; 
    w *= aa.w;
    return *this; 
}

template < class Real > 
const Vector4< Real >& Vector4< Real >::operator/=( const Vector4 &aa ) 
{ 
    x /= aa.x; 
    y /= aa.y; 
    z /= aa.z;
    w /= aa.w; 
    return *this; 
}

template < class Real > 
Real Vector4< Real >::dot( const Vector4 &aa ) const 
{ 
    return x * aa.x + y * aa.y + z * aa.z + w * aa.w; 
}

template < class Real > 
bool Vector4< Real >::operator==( const Vector4 &aa ) const 
{ 
    return ( x == aa.x && y == aa.y && z == aa.z && w = aa.w ); 
}

template < class Real > 
bool Vector4< Real >::operator!=(const Vector4 &aa ) const 
{ 
    return ( x != aa.x || y != aa.y || z != aa.z || w != aa.w ); 
}

template < class Real > 
void Vector4< Real >::invert() 
{	
    x = -x; 
    y = -y; 
    z = -z; 
    w = -w;
}

template < class Real > 
Real Vector4< Real >::max() 
{ 
    Real m = std::max( x, y ); 
    m = std::max( m, z); 
    m = std::max( m, w);
    return m; 
}

template < class Real > 
Real Vector4< Real >::min() 
{ 
    Real m = std::min( x, y ); 
    m = std::min( m, z); 
    m = std::min( m, w);
    return m; 
} 
}	
#endif
