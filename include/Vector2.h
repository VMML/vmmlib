/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Stefan Eilemann
*
* @license BSD license, check LICENSE
*/ 

#ifndef _vmml_Vector2_H_
#define _vmml_Vector2_H_

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cassert>

// - declaration -

namespace vmml
{

template< class T > 
class Vector2
{
public:
    union
    {
        struct
        {
            T  x, y;
        };
        T xy[2];
    };

    // contructors
    Vector2(); // warning: components NOT initialised ( for performance )
    Vector2( const T a ); 
    Vector2( const T x, const T y ); 
    
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector2< double >
    //the pointer 'values' must be a valid 2 component c array of the resp. type
    Vector2( const float* values );
    Vector2( const double* values );
    
    ~Vector2();

    void set( T x, T y );
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector2< double >
    //the pointer 'values' must be a valid 2 component c array of the resp. type
    void set( const float* values );
    void set( const double* values );

    const Vector2& operator=( T a ); 
    const Vector2& operator=( const Vector2& a ); 

    T& operator[]( size_t position);
    const T& operator[]( size_t position) const;

    // vector/scalar operations
    Vector2 operator+( const T  a ) const;
    Vector2 operator-( const T  a ) const; 
    Vector2 operator*( const T  a ) const;
    Vector2 operator/( T  a ) const;
     
    const Vector2& operator+=( T  a );
    const Vector2& operator-=( T  a );
    const Vector2& operator*=( T  a );
    const Vector2& operator/=( T  a ); 

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

    T length() const;
    T lengthSquared() const;

    T normalise();
    T normalise( const float* source );
    void scale( T scale_factor );
    void invert();

    T getMinComponent() const;
    T getMaxComponent() const;
    T getArea() const;

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

    static const Vector2 ZERO;
};

#ifndef VMMLIB_DISABLE_TYPEDEFS
    typedef Vector2<int>    Vector2i;
    typedef Vector2<float>  Vector2f;
    typedef Vector2<double> Vector2d;
#endif

// - implementation - //
       
template< typename T > 
const Vector2< T > Vector2< T >::ZERO( 0, 0 );

template < class T > 
Vector2< T >::Vector2() 
{} 

template < class T > 
Vector2< T >::Vector2( const T a )
    : x(a)
    , y(a)
{} 

template < class T > 
Vector2< T >::Vector2( const T i, const T j )
    : x(i)
    , y(j)
{} 

template < class T > 
Vector2< T >::Vector2( const float* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
}

template < class T > 
Vector2< T >::Vector2( const double* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
}


template < class T > 
Vector2< T >::~Vector2()
{}

template < class T > 
void Vector2< T >::set( T xx, T yy )
{ 
    x = xx; 
    y = yy; 
}

template < class T > 
void Vector2< T >::set( const float* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
}

template < class T > 
void Vector2< T >::set( const double* values )
{
    assert( values && "Vector2: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
}

template < class T > 
const Vector2< T >& Vector2< T >::operator=( T a )
{ 
    x = y = a; 
    return *this; 
} 

template < class T > 
const Vector2< T >& Vector2< T >::operator=( const Vector2& a ) 
{ 
    x = a.x; 
    y = a.y; 
    return *this;
} 


template < class T > 
T& Vector2< T >::operator[]( size_t position) 
{ 
    assert( position < 2 && "Vector2::operator[] Invalid component index!" ); 
    return xy[position];
}
         
template < class T > 
const T& Vector2< T >::operator[]( size_t position) const
{ 
    assert( position < 2 && "Vector2::operator[] Invalid component index!" ); 
    return xy[position]; 
} 
	
template < class T > 
T  Vector2< T >::length() const 
{ 
    T l = lengthSquared();
    return sqrt( l ); 
} 

template < class T > 
T  Vector2< T >::lengthSquared() const 
{ 
    return x * x + y * y; 
} 

template < class T > 
T  Vector2< T >::normalise()
{ 
    T l = length(); 
    if ( l == 0 ) 
        return 0; 
    l = 1.0f / l; 
    x *= l; 
    y *= l; 
    return l; 
} 

template < class T > 
T  Vector2< T >::normalise( const float* source )
{
    Vector2< T > a ( source );
    T l = a.length();
    if ( l == 0 ) 
        return 0;
    
    l = 1.0f / l;
    source[0] *= l;
    source[1] *= l;
    return l;
}


template < class T >
void Vector2< T >::scale( T scale_factor )
{
    operator*=( scale_factor );
}

template < class T > 
Vector2< T > Vector2< T >::operator+( const T  a ) const 
{ 
    return Vector2( x + a, y + a ); 
} 

template < class T > 
Vector2< T > Vector2< T >::operator-( const T  a ) const 
{ 
    return Vector2( x - a, y - a ); 
}
 
template < class T > 
Vector2< T > Vector2< T >::operator*( const T  a ) const 
{ 
    return Vector2( x * a, y * a ); 
}

template < class T > 
Vector2< T > Vector2< T >::operator/( T  a ) const 
{ 
    assert( a != 0.0f ); 
    a = 1.0f / a; 
    return Vector2( x * a, y * a ); 
}

template < class T > 
const Vector2< T >& Vector2< T >::operator+=( T  a ) 
{ 
    x += a; 
    y += a; 
    return *this; 
} 

template < class T > 
const Vector2< T >& Vector2< T >::operator-=( T  a ) 
{ 
    x -= a; 
    y -= a; 
    return *this; 
} 

template < class T > 
const Vector2< T >& Vector2< T >::operator*=( T  a ) 
{ 
    x *= a; 
    y *= a; 
    return *this; 
}
 
template < class T > 
const Vector2< T >& Vector2< T >::operator/=( T  a ) 
{ 
    a = 1.0f / a; 
    x *= a; 
    y *= a; 
    return *this; 
} 

// vector/vector operations
template < class T > 
Vector2< T > Vector2< T >::operator+( const Vector2 &a ) const 
{ 
    return Vector2( x + a.x, y + a.y ); 
}
 
template < class T > 
Vector2< T > Vector2< T >::operator-( const Vector2 &a ) const 
{ 
    return Vector2( x - a.x, y - a.y ); 
}

template < class T > 
Vector2< T > Vector2< T >::operator*( const Vector2 &a ) const 
{ 
    return Vector2( x * a.x, y * a.y ); 
} 

template < class T > 
Vector2< T > Vector2< T >::operator/( const Vector2 &a ) const 
{ 
    return Vector2( x / a.x, y / a.y ); 
} 

template < class T > 
Vector2< T > Vector2< T >::operator-() const 
{ 
    return Vector2( -x, -y );
}

template < class T > 
const Vector2< T >& Vector2< T >::operator+=( const Vector2 &a ) 
{ 
    x += a.x; 
    y += a.y; 
    return *this; 
} 

template < class T > 
const Vector2< T >& Vector2< T >::operator-=( const Vector2 &a ) 
{ 
    x -= a.x; 
    y -= a.y; 
    return *this; 
}

template < class T > 
const Vector2< T >& Vector2< T >::operator*=( const Vector2 &a ) 
{ 
    x *= a.x; 
    y *= a.y; 
    return *this; 
}

template < class T > 
const Vector2< T >& Vector2< T >::operator/=( const Vector2 &a ) 
{ 
    x /= a.x; 
    y /= a.y; 
    return *this; 
}

template < class T > 
bool Vector2< T >::operator==( const Vector2 &a ) const 
{ 
    return ( x == a.x && y == a.y ); 
}

template < class T > 
bool Vector2< T >::operator!=(const Vector2 &a ) const 
{ 
    return ( x != a.x || y != a.y ); 
}

template < class T > 
void Vector2< T >::invert() 
{	
    x = -x; 
    y = -y; 
}

template < class T > 
T Vector2< T >::getMaxComponent() const
{ 
    return std::max( x,y ); 
}

template < class T > 
T Vector2< T >::getMinComponent() const
{ 
    return std::min( x,y ); 
} 

template < class T > 
T Vector2< T >::getArea() const
{ 
    return x * y; 
} 
}	
#endif
