/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Stefan Eilemann
*
* @license BSD license, check LICENSE
*/ 

#ifndef _vmml_Vector2_H_
#define _vmml_Vector2_H_

#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <assert.h>

// - declaration -

namespace vmml
{

template< class Real > 
class Vector2
{
public:
    union
    {
        struct
        {
            Real  x, y;
        };
        Real xy[2];
    };

    // contructors
    Vector2(); // warning: components NOT initialised ( for performance )
    Vector2( const Real a ); 
    Vector2( const Real x, const Real y ); 
    
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector2< double >
    //the pointer 'values' must be a valid 2 component c array of the resp. type
    Vector2( const float* values );
    Vector2( const double* values );
    
    ~Vector2();

    void set( Real x, Real y );
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector2< double >
    //the pointer 'values' must be a valid 2 component c array of the resp. type
    void set( const float* values );
    void set( const double* values );

    const Vector2& operator=( Real a ); 
    const Vector2& operator=( const Vector2& a ); 

    Real& operator[]( size_t position);
    const Real& operator[]( size_t position) const;

    // vector/scalar operations
    Vector2 operator+( const Real  a ) const;
    Vector2 operator-( const Real  a ) const; 
    Vector2 operator*( const Real  a ) const;
    Vector2 operator/( Real  a ) const;
     
    const Vector2& operator+=( Real  a );
    const Vector2& operator-=( Real  a );
    const Vector2& operator*=( Real  a );
    const Vector2& operator/=( Real  a ); 

    // vector/vector operations
    Vector2 operator+( const Vector2& a ) const; 
    Vector2 operator-( const Vector2& a ) const;
    Vector2 operator*( const Vector2& a ) const; 
    Vector2 operator/( const Vector2& a ) const; 
    Vector2 operator-() const;

    const Vector2& operator+=( const Vector2& a ); 
    const Vector2& operator-=( const Vector2& a ); 
    const Vector2& operator*=( const Vector2& a ); 
    const Vector2& operator/=( const Vector2& a ); 

    bool operator==( const Vector2 &a ) const;
    bool operator!=(const Vector2 &a ) const;

    Real length() const;
    Real lengthSquared() const;

    Real normalise();
    Real normalise( const float* source );
    void scale( Real scale_factor );
    void invert();

    Real min();
    Real max();

    friend std::ostream& operator << ( std::ostream& os, const Vector2& v )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << "[" << std::setw(10) << v.x << " " << std::setw(10) << v.y << "]";
        os.precision( prec );
        os.setf( flags );
        return os;
    };        
};

#ifndef VMMLIB_DISABLE_TYPEDEFS
    typedef Vector2<int>    Vector2i;
    typedef Vector2<float>  Vector2f;
    typedef Vector2<double> Vector2d;
#endif

// - implementation - //
       
template < class Real > 
Vector2< Real >::Vector2() 
{} 

template < class Real > 
Vector2< Real >::Vector2( const Real a )
    : x(a)
    , y(a)
{} 

template < class Real > 
Vector2< Real >::Vector2( const Real i, const Real j )
    : x(i)
    , y(j)
{} 

template < class Real > 
Vector2< Real >::Vector2( const float* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
}

template < class Real > 
Vector2< Real >::Vector2( const double* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
}


template < class Real > 
Vector2< Real >::~Vector2()
{}

template < class Real > 
void Vector2< Real >::set( Real xx, Real yy )
{ 
    x = xx; 
    y = yy; 
}

template < class Real > 
void Vector2< Real >::set( const float* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
}

template < class Real > 
void Vector2< Real >::set( const double* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< Real > ( values[0] );
    y = static_cast< Real > ( values[1] );
}

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator=( Real a )
{ 
    x = y = a; 
    return *this; 
} 

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator=( const Vector2& a ) 
{ 
    x = a.x; 
    y = a.y; 
    return *this;
} 


template < class Real > 
Real& Vector2< Real >::operator[]( size_t position) 
{ 
    assert( position < 2 && "Vector2::operator[] Invalid component index!" ); 
    return xy[position];
}
         
template < class Real > 
const Real& Vector2< Real >::operator[]( size_t position) const
{ 
    assert( position < 2 && "Vector2::operator[] Invalid component index!" ); 
    return xy[position]; 
} 
	
template < class Real > 
Real  Vector2< Real >::length() const 
{ 
    Real l = lengthSquared();
    return sqrt( l ); 
} 

template < class Real > 
Real  Vector2< Real >::lengthSquared() const 
{ 
    return x * x + y * y; 
} 

template < class Real > 
Real  Vector2< Real >::normalise()
{ 
    Real l = length(); 
    if ( l == 0 ) 
        return 0; 
    l = 1.0f / l; 
    x *= l; 
    y *= l; 
    return l; 
} 

template < class Real > 
Real  Vector2< Real >::normalise( const float* source )
{
    Vector2< Real > a ( source );
    Real l = a.length();
    if ( l == 0 ) 
        return 0;
    
    l = 1.0f / l;
    source[0] *= l;
    source[1] *= l;
    return l;
}


template < class Real >
void Vector2< Real >::scale( Real scale_factor )
{
    operator*=( scale_factor );
}

template < class Real > 
Vector2< Real > Vector2< Real >::operator+( const Real  a ) const 
{ 
    return Vector2( x + a, y + a ); 
} 

template < class Real > 
Vector2< Real > Vector2< Real >::operator-( const Real  a ) const 
{ 
    return Vector2( x - a, y - a ); 
}
 
template < class Real > 
Vector2< Real > Vector2< Real >::operator*( const Real  a ) const 
{ 
    return Vector2( x * a, y * a ); 
}

template < class Real > 
Vector2< Real > Vector2< Real >::operator/( Real  a ) const 
{ 
    assert( a != 0.0f ); 
    a = 1.0f / a; 
    return Vector2( x * a, y * a ); 
}

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator+=( Real  a ) 
{ 
    x += a; 
    y += a; 
    return *this; 
} 

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator-=( Real  a ) 
{ 
    x -= a; 
    y -= a; 
    return *this; 
} 

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator*=( Real  a ) 
{ 
    x *= a; 
    y *= a; 
    return *this; 
}
 
template < class Real > 
const Vector2< Real >& Vector2< Real >::operator/=( Real  a ) 
{ 
    a = 1.0f / a; 
    x *= a; 
    y *= a; 
    return *this; 
} 

// vector/vector operations
template < class Real > 
Vector2< Real > Vector2< Real >::operator+( const Vector2 &a ) const 
{ 
    return Vector2( x + a.x, y + a.y ); 
}
 
template < class Real > 
Vector2< Real > Vector2< Real >::operator-( const Vector2 &a ) const 
{ 
    return Vector2( x - a.x, y - a.y ); 
}

template < class Real > 
Vector2< Real > Vector2< Real >::operator*( const Vector2 &a ) const 
{ 
    return Vector2( x * a.x, y * a.y ); 
} 

template < class Real > 
Vector2< Real > Vector2< Real >::operator/( const Vector2 &a ) const 
{ 
    return Vector2( x / a.x, y / a.y ); 
} 

template < class Real > 
Vector2< Real > Vector2< Real >::operator-() const 
{ 
    return Vector2( -x, -y );
}

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator+=( const Vector2 &a ) 
{ 
    x += a.x; 
    y += a.y; 
    return *this; 
} 

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator-=( const Vector2 &a ) 
{ 
    x -= a.x; 
    y -= a.y; 
    return *this; 
}

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator*=( const Vector2 &a ) 
{ 
    x *= a.x; 
    y *= a.y; 
    return *this; 
}

template < class Real > 
const Vector2< Real >& Vector2< Real >::operator/=( const Vector2 &a ) 
{ 
    x /= a.x; 
    y /= a.y; 
    return *this; 
}

template < class Real > 
bool Vector2< Real >::operator==( const Vector2 &a ) const 
{ 
    return ( x == a.x && y == a.y ); 
}

template < class Real > 
bool Vector2< Real >::operator!=(const Vector2 &a ) const 
{ 
    return ( x != a.x || y != a.y ); 
}

template < class Real > 
void Vector2< Real >::invert() 
{	
    x = -x; 
    y = -y; 
}

template < class Real > 
Real Vector2< Real >::max() 
{ 
    return std::max( x,y ); 
}

template < class Real > 
Real Vector2< Real >::min() 
{ 
    return std::min( x,y ); 
} 
}	
#endif
