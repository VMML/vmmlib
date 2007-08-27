/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Jonas Boesch
* @author Stefan Eilemann
* @author Renato Pajarola
* @author David H. Eberly ( Wild Magic )
* @author Andrew Willmott ( VL )
*
* @license revised BSD license, check LICENSE
*
* parts of the source code of VMMLib were inspired by David Eberly's 
* Wild Magic and Andrew Willmott's VL.
* 
*/ 

#ifndef __VMML__VECTOR3__H__
#define __VMML__VECTOR3__H__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

// - declaration -

namespace vmml
{
template< typename T > class Vector2;
template< typename T > class Vector4;

template< typename T > 
class Vector3
{
public:
    union
    {
        struct
        {
            T  x, y, z;
        };
        struct
        {
            T  r, g, b;
        };
        T xyz[3];
        T rgb[3];
    };

    // contructors
    Vector3(); // warning: components NOT initialised ( for performance )
    Vector3( const T  a ); 
    Vector3( const T  i, const T  j, const T  k ); 
    Vector3( const Vector2< T >& xy, const T z ); 
    Vector3( const Vector4< T >& from ); 

    // type conversion constructor
    template< typename U >
    Vector3( const Vector3< U >& orig );
    
    // dangerous, but implemented to allow easy conversion between 
    // Vector< float > and Vector3< double >
    //the pointer 'values must be a valid 3 component c array of the resp. type
    Vector3( const float* values );
    Vector3( const double* values );
    
    ~Vector3();

    void set( T xx, T yy, T zz );
    // the pointer 'values' must be a valid 3 component c array of the resp. type
    void set( const float* values );
    void set( const double* values );

    const Vector3& operator=( T a ); 
    const Vector3& operator=( const Vector3& a ); 
    template< typename U >
    const Vector3& operator=( const Vector3< U >& a );

    T& operator[]( size_t position);
    const T& operator[]( size_t position) const;

    // vector/scalar operations
    Vector3 operator+( const T  a ) const;
    Vector3 operator-( const T  a ) const; 
    Vector3 operator*( const T  a ) const;
    Vector3 operator/( T  a ) const;
     
    const Vector3& operator+=( T  a );
    const Vector3& operator-=( T  a );
    const Vector3& operator*=( T  a );
    const Vector3& operator/=( T  a ); 

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
    bool isAkin(  const Vector3& a, 
                  const T& delta = std::numeric_limits<T>::epsilon() );

    // component-component compare
    // returns a size_t with a bitmask of the component comparison results
    // -> if this->xy[k] is smaller than a[k], the kth bit will be enabled;
    const size_t smaller( const Vector3& a ) const;
    const size_t smaller( const Vector3& a, const size_t axis ) const;
    // -> if this->xy[k] is smaller than a[k], the kth bit will be enabled;
    const size_t greater( const Vector3& a ) const;
    const size_t greater( const Vector3& a, const size_t axis ) const;
    

    T length() const;
    T lengthSquared() const;

    T distance( const Vector3& other ) const;
    T distanceSquared( const Vector3& other ) const;

    T normalize();
    static T normalize( float* source );
    // deprecated
    T normalise();
    static T normalise( float* source );

    void scale( T scale_factor );

    // result = vec1.cross( vec2 ) => vec1 x vec2
    Vector3 cross( const Vector3& a ) const;
    // result.cross( vec1, vec2 ) => vec1 x vec2
    void cross( const Vector3 &a, const Vector3 &b);
    T dot( const Vector3& a) const;
    static T dot( const Vector3& a, const Vector3& b);

    void invert(); 

    // *this is the result
    void normal( const Vector3& aa, const Vector3& bb, const Vector3& cc );
    //returns the normal of *this and the two argument vectors
    Vector3 normal( const Vector3& aa, const Vector3& bb );

    Vector3 rotate( T theta, T  rx, T  ry, T  rz ); 

    T getMinComponent();
    T getMaxComponent();

    friend std::ostream& operator << ( std::ostream& os, const Vector3& v )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << "[" << std::setw(10) << v.x << " " << std::setw(10) << v.y 
           << " " << std::setw(10) << v.z << "]";
        os.precision( prec );
        os.setf( flags );
        return os;
    }; 

    static const Vector3 ZERO;

}; // class vector3
    

typedef Vector3< float >            Vector3f;
typedef Vector3< double >           Vector3d;
typedef Vector3< int >              Vector3i;
typedef Vector3< unsigned char >    Vector3ub;    

typedef Vector3< float >            vec3f;
typedef Vector3< double >           vec3d;
typedef Vector3< int >              vec3i;
typedef Vector3< unsigned char >    vec3ub;    

} // namespace vmml

// - implementation - //
#include <vmmlib/vector2.h>
#include <vmmlib/vector4.h>
 
namespace vmml
{

template< typename T > 
const Vector3< T > Vector3< T >::ZERO( 0, 0, 0 );

template < typename T > 
Vector3< T >::Vector3() 
{} 

template < typename T > 
Vector3< T >::Vector3( const T  a )
    : x(a)
    , y(a)
    , z(a) 
{} 

template < typename T > 
Vector3< T >::Vector3( const T  i, const T  j, const T  k )
    : x(i)
    , y(j)
    , z(k) 
{} 



template < typename T > 
Vector3< T >::Vector3( const Vector4<T>& from )
{
    const T wInv = 1./from.w;
    x = from.x * wInv;
    y = from.y * wInv;
    z = from.z * wInv;
} 



template < typename T > 
Vector3< T >::Vector3( const Vector2< T >& xy, const T z_ )
    : x( xy.x )
    , y( xy.y )
    , z( z_ )
{} 



template < typename T > 
template < typename U >
Vector3< T >::Vector3( const Vector3< U >& a )
    : x( static_cast< T >( a.x ) )
    , y( static_cast< T >( a.y ) )
    , z( static_cast< T >( a.z ) )
{}



template < typename T > 
Vector3< T >::Vector3( const float* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
}



template < typename T > 
Vector3< T >::Vector3( const double* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
}



template < typename T > 
Vector3< T >::~Vector3()
{}



template < typename T > 
void Vector3< T >::set( T xx, T yy, T zz )
{ 
    x = xx; 
    y = yy; 
    z = zz; 
}



template < typename T > 
void Vector3< T >::set( const float* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
}



template < typename T > 
void Vector3< T >::set( const double* values )
{
    assert( values && "Vector3: Nullpointer argument as source for initialisation!" );
    x = static_cast< T > ( values[0] );
    y = static_cast< T > ( values[1] );
    z = static_cast< T > ( values[2] );
}



template < typename T > 
const Vector3< T >& Vector3< T >::operator=( T a )
{ 
    x = y = z = a; 
    return *this; 
} 



template < typename T > 
const Vector3< T >& Vector3< T >::operator=( const Vector3& a ) 
{ 
    x = a.x; 
    y = a.y; 
    z = a.z; 
    return *this;
} 



template < typename T > 
template < typename U >
const Vector3< T >& Vector3< T >::operator=( const Vector3< U >& orig )
{
    x = static_cast< T > ( orig.x );
    y = static_cast< T > ( orig.y );
    z = static_cast< T > ( orig.z );
    return *this;
}



template < typename T > 
T& Vector3< T >::operator[]( size_t index ) 
{ 
    assert( index < 3 && "Vector3::operator[] Invalid component index!" ); 
    return xyz[ index ]; 
}
       



template < typename T > 
const T& Vector3< T >::operator[]( size_t index ) const
{ 
    assert( index < 3 && "Vector3::operator[] Invalid component index!" ); 
    return xyz[ index ]; 
} 


	
template < typename T > 
T  Vector3< T >::length() const 
{ 
    return sqrt( lengthSquared( ) ); 
} 

template <> 
inline float Vector3< float >::length() const 
{ 
    return sqrtf( lengthSquared( )); 
} 



template < typename T > 
T  Vector3< T >::lengthSquared() const 
{ 
    return x * x + y * y + z * z; 
} 



template < typename T > 
T 
Vector3< T >::distance( const Vector3& other ) const
{
    return (*this - other).length();
}



template < typename T > 
T 
Vector3< T >::distanceSquared( const Vector3& other ) const
{
    return (*this - other).lengthSquared();
}



template < typename T > 
T  Vector3< T >::normalise()
{ 
    const T l = length(); 
    if ( l == 0 ) 
        return 0; 
    x /= l; 
    y /= l; 
    z /= l; 
    return l; 
} 



// PRECONDITION: float* source is a valid 3-float array
template < typename T > 
T  Vector3< T >::normalise( float* source )
{
    Vector3< float >* a = ( Vector3< float >* ) source;
    const T l = a->length();
    if ( l == 0 ) 
        return 0;

    source[0] /= l;
    source[1] /= l;
    source[2] /= l;
    return l;
}



template < typename T > 
T  Vector3< T >::normalize()
{ 
    const T l = length();
    if ( l == 0 ) 
        return 0; 
    x /= l; 
    y /= l; 
    z /= l; 
    return l; 
} 



// PRECONDITION: float* source is a valid 3-float array
template < typename T > 
T  Vector3< T >::normalize( float* source )
{
    Vector3< float >* a = ( Vector3< float >* ) source;
    const T l = a->length();
    if ( l == 0 ) 
        return 0;
    
    source[0] /= l;
    source[1] /= l;
    source[2] /= l;
    return l;
}



template < typename T >
void Vector3< T >::scale( T scale_factor )
{
    operator*=( scale_factor );
}



template < typename T > 
Vector3< T > Vector3< T >::operator+( const T  a ) const 
{ 
    return Vector3( x + a, y + a, z + a ); 
} 



template < typename T > 
Vector3< T > Vector3< T >::operator-( const T  a ) const 
{ 
    return Vector3( x - a, y - a, z - a ); 
}



template < typename T > 
Vector3< T > Vector3< T >::operator*( const T  a ) const 
{ 
    return Vector3( x * a, y * a, z * a ); 
}



template < typename T > 
Vector3< T > Vector3< T >::operator/( T  a ) const 
{ 
    assert( a != 0.0f ); 
    a = 1.0f / a; 
    return Vector3( x * a, y * a, z * a ); 
}



template < typename T > 
const Vector3< T >& Vector3< T >::operator+=( T  a ) 
{ 
    x += a; 
    y += a; 
    z += a; 
    return *this; 
} 



template < typename T > 
const Vector3< T >& Vector3< T >::operator-=( T  a ) 
{ 
    x -= a; 
    y -= a; 
    z -= a; 
    return *this; 
} 



template < typename T > 
const Vector3< T >& Vector3< T >::operator*=( T  a ) 
{ 
    x *= a; 
    y *= a; 
    z *= a; 
    return *this; 
}



template < typename T > 
const Vector3< T >& Vector3< T >::operator/=( T  a ) 
{ 
    a = 1.0f / a; 
    x *= a; 
    y *= a; 
    z *= a; 
    return *this; 
} 



// vector/vector operations
template < typename T > 
Vector3< T > Vector3< T >::operator+( const Vector3 &a ) const 
{ 
    return Vector3( x + a.x, y + a.y, z + a.z ); 
}



template < typename T > 
Vector3< T > Vector3< T >::operator-( const Vector3 &a ) const 
{ 
    return Vector3( x - a.x, y - a.y, z - a.z ); 
}



template < typename T > 
Vector3< T > Vector3< T >::operator*( const Vector3 &a ) const 
{ 
    return Vector3( x * a.x, y * a.y, z * a.z ); 
} 



template < typename T > 
Vector3< T > Vector3< T >::operator/( const Vector3 &a ) const 
{ 
    return Vector3( x / a.x, y / a.y, z / a.z ); 
} 



template < typename T > 
Vector3< T > Vector3< T >::operator-() const 
{ 
    return Vector3( -x, -y, -z );
}



template < typename T > 
const Vector3< T >& Vector3< T >::operator+=( const Vector3 &a ) 
{ 
    x += a.x; 
    y += a.y; 
    z += a.z; 
    return *this; 
} 



template < typename T > 
const Vector3< T >& Vector3< T >::operator-=( const Vector3 &a ) 
{ 
    x -= a.x; 
    y -= a.y; 
    z -= a.z; 
    return *this; 
}



template < typename T > 
const Vector3< T >& Vector3< T >::operator*=( const Vector3 &a ) 
{ 
    x *= a.x; 
    y *= a.y; 
    z *= a.z; 
    return *this; 
}



template < typename T > 
const Vector3< T >& Vector3< T >::operator/=( const Vector3 &a ) 
{ 
    x /= a.x; 
    y /= a.y; 
    z /= a.z; 
    return *this; 
}



template < typename T > 
const size_t 
Vector3< T >::smaller( const Vector3< T >& a ) const
{
    size_t result = 0;
    if ( x < a.x )
        result |= 1;
    if ( y < a.y )
        result |= 2;
    if ( z < a.z )
        result |= 4;
    return result;
}



template < typename T > 
const size_t 
Vector3< T >::smaller( const Vector3< T >& a, const size_t axis  ) const
{
    return ( xyz[ axis ] < a.xyz[ axis ] ) ? 1 << axis : 0;
}



template < typename T > 
const size_t 
Vector3< T >::greater( const Vector3< T >& a ) const
{
    size_t result = 0;
    if ( x > a.x )
        result |= 1;
    if ( y > a.y )
        result |= 2;
    if ( z > a.z )
        result |= 4;
    return result;
}



template < typename T > 
const size_t 
Vector3< T >::greater( const Vector3< T >& a, const size_t axis  ) const
{
    return ( xyz[ axis ] > a.xyz[ axis ] ) ? 1 << axis : 0;
}



// result = vec1.cross( vec2 ) => vec1 x vec2
template < typename T > 
Vector3< T > Vector3< T >::cross( const Vector3& a ) const
{ 
    return Vector3( y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x ); 
}



// result.cross( vec1, vec2 ) => vec1 x vec2
template < typename T > 
void Vector3< T >::cross( const Vector3 &a, const Vector3 &b)
{ 
    x = a.y * b.z - a.z * b.y; 
    y = a.z * b.x - a.x * b.z; 
    z = a.x * b.y - a.y * b.x; 
}



template < typename T > 
T Vector3< T >::dot( const Vector3& a) const 
{ 
    return x * a.x + y * a.y + z * a.z; 
}



template < typename T > 
T Vector3< T >::dot( const Vector3& a, const Vector3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}



template < typename T > 
bool Vector3< T >::operator==( const Vector3& a ) const 
{ 
    return ( x == a.x && y == a.y && z == a.z ); 
}



template < typename T > 
bool Vector3< T >::operator!=(const Vector3& a ) const 
{ 
    return ( x != a.x || y != a.y || z != a.z ); 
}



template < typename T > 
bool Vector3< T >::isAkin( const Vector3& a, const T& delta )
{
    if( fabs( x-a.x ) > delta || fabs( y-a.y ) > delta || 
        fabs( z-a.z ) > delta )

        return false;
    return true;
}



template< > 
inline bool Vector3< float >::isAkin( const Vector3& a, const float& delta )
{
    if( fabsf( x-a.x ) > delta || fabsf( y-a.y ) > delta || 
        fabsf( z-a.z ) > delta )

        return false;
    return true;
}



template < typename T > 
void Vector3< T >::invert() 
{	
    x = -x; 
    y = -y; 
    z = -z; 
}



// *this is the result
template< typename T >
void Vector3< T >::normal( const Vector3< T >& aa, 
                              const Vector3< T >& bb,
                              const Vector3< T >& cc )
{
    Vector3< T > u,v;

    // right hand system, CCW triangle
    u = bb - aa;
    v = cc - aa;
    cross( u, v );
    normalize();
}



//returns the normal of *this and the two argument vectors
template< typename T >
Vector3< T > Vector3< T >::normal( const Vector3< T >& aa, 
                                         const Vector3< T >& bb )
{
    Vector3< T > tmp;
    tmp.normal( *this, aa, bb);
    return tmp;
}




template < typename T > 
Vector3< T > Vector3< T >::rotate( T theta, T rx, T ry, 
                                         T rz )      
{   
    Vector3 q( 0, 0, 0 ), r( rx, ry, rz );
    r.normalize();
    const T costheta = ( T ) cos( theta );
    const T sintheta = ( T ) sin( theta );
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




template<>
inline Vector3< float > Vector3< float >::rotate( float theta, float rx, 
                                                  float ry, float rz )      
{   
    Vector3 q( 0, 0, 0 ), r( rx, ry, rz );
    r.normalize();
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




template < typename T > 
T Vector3< T >::getMaxComponent() 
{ 
    T m = std::max( x,y ); 
    m = std::max( m, z); 
    return m; 
}




template < typename T > 
T Vector3< T >::getMinComponent() 
{ 
    T m = std::min( x,y ); 
    m = std::min( m, z); 
    return m; 
} 



} //namespace vmml
	
#endif
