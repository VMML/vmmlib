#ifndef __VMML__VECTOR__HPP__
#define __VMML__VECTOR__HPP__

#include <vmmlib/exception.hpp>
#include <vmmlib/vmmlib_config.hpp>
#include <vmmlib/details.hpp>

#include <iostream>
#include <iomanip>
#include <vector>

namespace vmml
{

template< size_t M, typename float_t = double >
class vector
{
public:
    // accessors 
    inline float_t& operator()( size_t index );
    inline const float_t& operator()( size_t index ) const;
    inline float_t& at( size_t index );
    inline const float_t& at( size_t index ) const;

    inline float_t& operator[]( size_t index );
    inline const float_t& operator[]( size_t index ) const;

    // element accessors for M <= 4;
    inline float_t& x();
    inline float_t& y();
    inline float_t& z();
    inline float_t& w();
    inline const float_t& x() const;
    inline const float_t& y() const;
    inline const float_t& z() const;
    inline const float_t& w() const;

    bool operator==( const vector< M, float_t >& other ) const;
    bool operator!=( const vector< M, float_t >& other ) const;
    bool isEqualTo( const vector< M, float_t >& other, float_t tolerance = 1e-15 ) const;
   
    const vector< M, float_t >& operator=( const float_t* c_array );
    float_t operator=( float_t filler );

    const vector< M, float_t >& operator=( const vector< M, float_t >& other );
    
    // returns void to avoid 'silent' loss of precision when chaining
    template< typename other_float_t >
    void operator=( const vector< M, other_float_t >& other );
    
    vector< M, float_t > operator*( const vector< M, float_t >& other ) const;
    vector< M, float_t > operator/( const vector< M, float_t >& other ) const;    
    vector< M, float_t > operator+( const vector< M, float_t >& other ) const; 
    vector< M, float_t > operator-( const vector< M, float_t >& other ) const;

    void operator*=( const vector< M, float_t >& other );
    void operator/=( const vector< M, float_t >& other );    
    void operator+=( const vector< M, float_t >& other ); 
    void operator-=( const vector< M, float_t >& other );

    vector< M, float_t > operator*( const float_t other ) const;
    vector< M, float_t > operator/( const float_t other ) const;    
    vector< M, float_t > operator+( const float_t other ) const; 
    vector< M, float_t > operator-( const float_t other ) const;

    void operator*=( const float_t other );
    void operator/=( const float_t other );    
    void operator+=( const float_t other ); 
    void operator-=( const float_t other );

    vector< M, float_t > operator-() const;

    void invert();

    // constructors 
    vector() {}; // std ctor - WARNING: NO INITIALIZATION
    vector( float_t a ); // sets all components to a;

    // WARNING: the following constructors will not work for vectors where
    // M != number of arguments. Instead, the compiler will throw an error such as
    // 'no matching function for call to 'number_of_parameters_must_be_M()'.
    vector( float_t x, float_t y );
    vector( float_t x, float_t y, float_t z );
    vector( float_t x, float_t y, float_t z, float_t w );
    
    void set( float_t a ); // sets all components to a;
    void set( const vector< M-1, float_t >& v, float_t a );

    // WARNING: the following set functions will not work for vectors where
    // M != number of arguments. Instead, the compiler will throw an error such as
    // 'no matching function for call to 'number_of_parameters_must_be_M()'.
    void set( float_t x, float_t y );
    void set( float_t x, float_t y, float_t z ); 
    void set( float_t x, float_t y, float_t z, float_t w );
    
    // compute the cross product of two vectors
    // note: there's also a free function:
    // vector<> cross( const vector<>, const vector<> )

    // result = vec1.cross( vec2 ) => retval result = vec1 x vec2
    inline vector< M, float_t > cross( const vector< M, float_t >& rhs ) const;

    // result.cross( vec1, vec2 ) => (this) = vec1 x vec2
    void cross( const vector< M, float_t >& a, const vector< M, float_t >& b );
  
    // compute the dot product of two vectors
    // note: there's also a free function:
    // float_t dot( const vector<>, const vector<> );
    inline float_t dot( const vector< M, float_t >& other ) const;

    // normalize the vector
    // note: there's also a free function:
    // vector<> normalize( const vector<> );
    inline void normalize();
    vector< M, float_t > getNormalized() const;
    
    // L2 norm
    inline float_t norm() const;
    inline float_t normSquared() const;
    inline float_t length() const;
    inline float_t lengthSquared() const;
    
    inline float_t distance( const vector< M, float_t >& other ) const;
    inline float_t distanceSquared( const vector< M, float_t >& other ) const;
    
    // (*this) = normal of v0,v1,v2
    void computeNormal( 
        const vector< M, float_t >& v0, 
        const vector< M, float_t >& v1, 
        const vector< M, float_t >& v2
        );
        
    // retval = normal of (this), v1, v2
    vector< M, float_t > computeNormal(
        const vector< M, float_t >& v1, 
        const vector< M, float_t >& v2
        ) const;
    
    
	// sphere functions - sphere layout: center xyz, radius w
	inline vector< 3, float_t > projectPointOntoSphere( const vector< 3, float_t >& point ) const;
	// returns a negative distance if the point lies in the sphere
	inline float_t getDistanceToSphere( const vector< 3, float_t >& point ) const;

    inline vector< 3, float_t >& getSphereCenter();
    inline const vector< 3, float_t >& getSphereCenter() const;
    

	// plane functions - plane layout; normal xyz, distance w
    inline vector< 3, float_t > projectPointOntoPlane( const vector< 3, float_t >& point ) const;
	inline float_t getDistanceToPlane( const vector< 3, float_t >& point ) const;
    // this normalizes the vector3 consisting of the first three 
    // coordinates, and leaves the fourth as is.
    void normalizePlane();
    
    
    size_t  getSmallestComponentIndex() const;
    size_t  getLargestComponentIndex() const;

    float_t&        getSmallestComponent();
    float_t&        getLargestComponent();
    const float_t&  getSmallestComponent() const;
    const float_t&  getLargestComponent() const;
    
    
    void copyFrom1DimCArray( const float_t* c_array );

    template< typename different_float_t >
    void copyFrom1DimCArray( const different_float_t* c_array );
    
    size_t getM() const;
    
    friend std::ostream& operator << ( std::ostream& os, 
        const vector< M, float_t >& vector_ )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        
        os << " |";
        for( size_t index = 0; index < M; ++index )
        {
            os << std::setw(10) << vector_.at( index ) << " ";
        }
        os << "| ";
        os << std::endl;
        os.precision( prec );
        os.setf( flags );
        return os;
    };  


    // storage
    float_t array[ M ]
    #ifndef VMMLIB_DONT_FORCE_ALIGNMENT
        #ifdef _GCC
            __attribute__((aligned(16)))
        #endif
    #endif
    ;

}; // class vector


//
//  some free functions for convenience
//

// allows float * vector, not only vector * float 
template< size_t M, typename float_t >
static vector< M, float_t >
operator* ( float_t factor, const vector< M, float_t >& vector_ )
{
    return vector_ * factor;
}


template< size_t M, typename float_t >
inline float_t
dot( const vector< M, float_t >& first, const vector< M, float_t >& second ) 
{
    float_t tmp = 0.0;
    for( size_t index = 0; index < M; ++index )
    {
        tmp += first.at( index ) * second.at( index );
    }
    return tmp;
}


template< size_t M, typename float_t >
inline vector< M, float_t >
cross( const vector< M, float_t >& a, const vector< M, float_t >& b )
{
    return a.cross( b );
}


template< size_t M, typename float_t >
inline vector< M, float_t >
normalize( const vector< M, float_t >& vector_ )
{
    return vector_.getNormalized();
}



#if 0
typedef vector< 2, float >  vec2f;
typedef vector< 2, double > vec2d;
typedef vector< 3, float >  vec3f;
typedef vector< 3, double > vec3d;
typedef vector< 4, float >  vec4d;
typedef vector< 4, double > vec4d;
#endif 


template< size_t M, typename float_t >
vector< M, float_t >::vector( float_t a )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) = a;
    }

}



template< size_t M, typename float_t >
vector< M, float_t >::vector( float_t x, float_t y )
{
    details::number_of_parameters_must_be_M< 2, M, vector< M, float_t > >();

    array[ 0 ] = x;
    array[ 1 ] = y;
}


template< size_t M, typename float_t >
vector< M, float_t >::vector( float_t x, float_t y, float_t z )
{
    details::number_of_parameters_must_be_M< 3, M, vector< M, float_t > >();

    array[ 0 ] = x;
    array[ 1 ] = y;
    array[ 2 ] = z;
}



template< size_t M, typename float_t >
vector< M, float_t >::vector( float_t x, float_t y, float_t z, float_t w )
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();
    
    array[ 0 ] = x;
    array[ 1 ] = y;
    array[ 2 ] = z;
    array[ 3 ] = w;

}



template< size_t M, typename float_t >
void
vector< M, float_t >::set( float_t a )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) = a;
    }

}



template< size_t M, typename float_t >
void
vector< M, float_t >::set( const vector< M-1, float_t >& v, float_t a )
{
    memcpy( array, v.array, sizeof( float_t ) * (M-1) );
    at( M-1 ) = a;
}




template< size_t M, typename float_t >
void
vector< M, float_t >::set( float_t x, float_t y )
{
    details::number_of_parameters_must_be_M< 2, M, vector< M, float_t > >();

    array[ 0 ] = x;
    array[ 1 ] = y;
}


template< size_t M, typename float_t >
void
vector< M, float_t >::set( float_t x, float_t y, float_t z )
{
    details::number_of_parameters_must_be_M< 3, M, vector< M, float_t > >();

    array[ 0 ] = x;
    array[ 1 ] = y;
    array[ 2 ] = z;
}



template< size_t M, typename float_t >
void
vector< M, float_t >::set( float_t x, float_t y, float_t z, float_t w )
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();
    
    array[ 0 ] = x;
    array[ 1 ] = y;
    array[ 2 ] = z;
    array[ 3 ] = w;

}


template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::operator()( size_t index )
{
	return at( index );
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::operator()( size_t index ) const
{
	return at( index );
}



template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::at( size_t index )
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( index >= M )
    {
        VMMLIB_ERROR( "at() - index out of bounds", VMMLIB_HERE );
    }
    #endif
    return array[ index ];
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::at( size_t index ) const
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( index >= M )
    {
        VMMLIB_ERROR( "at() - index out of bounds", VMMLIB_HERE );
    }
    #endif
    return array[ index ];
}



template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::operator[]( size_t index )
{
    return at( index );
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::operator[]( size_t index ) const
{
    return at( index );
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator*( const vector< M, float_t >& other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) * other.at( index );
    }
    return result;
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator/( const vector< M, float_t >& other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) / other.at( index );
    }
    return result;
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator+( const vector< M, float_t >& other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) + other.at( index );
    }
    return result;
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator-( const vector< M, float_t >& other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) - other.at( index );
    }
    return result;
}




template< size_t M, typename float_t >
void
vector< M, float_t >::operator*=( const vector< M, float_t >& other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) *= other.at( index );
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::operator/=( const vector< M, float_t >& other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) /= other.at( index );
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::operator+=( const vector< M, float_t >& other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) += other.at( index );
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::operator-=( const vector< M, float_t >& other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) -= other.at( index );
    }
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator*( const float_t other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) * other;
    }
    return result;
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator/( const float_t other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) / other;
    }
    return result;
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator+( const float_t other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) + other;
    }
    return result;
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator-( const float_t other ) const
{
    vector< M, float_t > result;
    for( size_t index = 0; index < M; ++index )
    {
        result.at( index ) = at( index ) - other;
    }
    return result;
}




template< size_t M, typename float_t >
void
vector< M, float_t >::operator*=( const float_t other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) *= other;
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::operator/=( const float_t other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) /= other;
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::operator+=( const float_t other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) += other;
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::operator-=( const float_t other )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) -= other;
    }
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::operator-() const
{
    vector< M, float_t > v;
    for( size_t index = 0; index < M; ++index )
    {
        v.array[ index ] = -array[ index ];
    }
    return v;
}



template< size_t M, typename float_t >
void
vector< M, float_t >::invert()
{
    for( size_t index = 0; index < M; ++index )
    {
        array[ index ] = -array[ index ];
    }
}



template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::x()
{
    details::number_of_parameters_must_be_at_least_M< 1, M, vector< M, float_t > >();
    return array[ 0 ];
}



template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::y()
{
    details::number_of_parameters_must_be_at_least_M< 2, M, vector< M, float_t > >();
    return array[ 1 ];
}



template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::z()
{
    details::number_of_parameters_must_be_at_least_M< 3, M, vector< M, float_t > >();
    return array[ 2 ];
}



template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::w()
{
    details::number_of_parameters_must_be_at_least_M< 4, M, vector< M, float_t > >();
    return array[ 3 ];
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::x() const
{
    details::number_of_parameters_must_be_at_least_M< 1, M, vector< M, float_t > >();
    return array[ 0 ];
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::y() const
{
    details::number_of_parameters_must_be_at_least_M< 2, M, vector< M, float_t > >();
    return array[ 1 ];
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::z() const
{
    details::number_of_parameters_must_be_at_least_M< 3, M, vector< M, float_t > >();
    return array[ 2 ];
}



template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::w() const
{
    details::number_of_parameters_must_be_at_least_M< 4, M, vector< M, float_t > >();
    return array[ 3 ];
}



// result = vec1.cross( vec2 ) => result = vec1 x vec2
template< size_t M, typename float_t >
inline vector< M, float_t >
vector< M, float_t >::cross( const vector< M, float_t >& rhs ) const
{
    vector< M, float_t > result;
    result.cross( *this, rhs );
    return result;
}



// result.cross( vec1, vec2 ) => (this) = vec1 x vec2
template< size_t M, typename float_t >
void
vector< M, float_t >::
cross( const vector< M, float_t >& aa, const vector< M, float_t >& bb )
{ 
    details::number_of_parameters_must_be_M< 3, M, vector< M, float_t > >();

    array[ 0 ] = aa.y() * bb.z() - aa.z() * bb.y(); 
    array[ 1 ] = aa.z() * bb.x() - aa.x() * bb.z(); 
    array[ 2 ] = aa.x() * bb.y() - aa.y() * bb.x(); 
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::dot( const vector< M, float_t >& other ) const
{
    float_t tmp = 0.0;
    for( size_t index = 0; index < M; ++index )
    {
        tmp += at( index ) * other.at( index );
    }
    return tmp;
}


template< size_t M, typename float_t >
inline void
vector< M, float_t >::normalize()
{
    float_t norm_reciprocal = 1.0 / norm();
    this->operator*=( norm_reciprocal );
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::getNormalized() const
{
    vector< M, float_t > n( *this );
    n.normalize();
    return n;
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::norm() const
{
    return details::getSquareRoot( normSquared() );
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::normSquared() const
{
    float_t tmp = 0.0;
    for( size_t index = 0; index < M; ++index )
    {
        tmp += at( index ) * at( index );
    }
    return tmp;
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::length() const
{
    return norm();
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::lengthSquared() const
{
    return normSquared();
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::distance( const vector< M, float_t >& other ) const
{
    return details::getSquareRoot( distanceSquared() );
}



template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::distanceSquared( const vector< M, float_t >& other ) const
{
    vector< M, float_t > tmp( *this );
    tmp -= other;
    return tmp.lengthSquared();
}



template< size_t M, typename float_t >
void
vector< M, float_t >::computeNormal(
    const vector< M, float_t >& aa, 
    const vector< M, float_t >& bb, 
    const vector< M, float_t >& cc
    )
{
    vector< M, float_t > u,v;
    // right hand system, CCW triangle
    u = bb - aa;
    v = cc - aa;
    cross( u, v );
    normalize();
}



template< size_t M, typename float_t >
vector< M, float_t >
vector< M, float_t >::computeNormal(
    const vector< M, float_t >& bb, 
    const vector< M, float_t >& cc
    ) const
{
    vector< M, float_t > tmp;
    tmp.computeNormal( *this, bb, cc);
    return tmp;
}



// sphere layout: center xyz, radius w
template< size_t M, typename float_t >
inline vector< 3, float_t >
vector< M, float_t >::
projectPointOntoSphere( const vector< 3, float_t >& point ) const
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();

    const vector< 3, float_t >& center = reinterpret_cast< const vector< 3, float_t >& >( *this );
    vector< 3, float_t > projPoint( point );
    projPoint -= center;
    projPoint.normalize();
    projPoint *= w();
    return center + projPoint;
}



// sphere layout: center xyz, radius w
template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::
getDistanceToSphere( const vector< 3, float_t >& point ) const
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();

    const vector< 3, float_t >& center_ = reinterpret_cast< const vector< 3, float_t >& >( *this );
	return ( point - center_ ).length() - w();
}



// sphere layout: center xyz, radius w
template< size_t M, typename float_t >
inline vector< 3, float_t >&
vector< M, float_t >::getSphereCenter()
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();
    return reinterpret_cast< vector< 3, float_t >& >( *this );
}



// sphere layout: center xyz, radius w
template< size_t M, typename float_t >
inline const vector< 3, float_t >&
vector< M, float_t >::getSphereCenter() const
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();
    return reinterpret_cast< vector< 3, float_t >& >( *this );
}



// plane: normal xyz, distance w
template< size_t M, typename float_t >
inline float_t
vector< M, float_t >::getDistanceToPlane( const vector< 3, float_t >& point ) const
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();

    const vector< 3, float_t >& normal = reinterpret_cast< const vector< 3, float_t >& >( *this );
    return normal.dot( point ) + w();
}



// plane: normal xyz, distance w
template< size_t M, typename float_t >
vector< 3, float_t >
vector< M, float_t >::projectPointOntoPlane( const vector< 3, float_t >& point ) const
{
    details::number_of_parameters_must_be_M< 4, M, vector< M, float_t > >();

    const vector< 3, float_t >& normal = reinterpret_cast< const vector< 3, float_t >& >( *this );
    return point - ( normal * getDistanceToPlane( point ) );
}



template< size_t M, typename float_t >
bool
vector< M, float_t >::operator==( const vector< M, float_t >& other ) const
{
    bool ok = true;
    for( size_t index = 0; ok && index < M; ++index )
    {
        ok = at( index ) == other.at( index );
    }
    return ok;
}


template< size_t M, typename float_t >
bool
vector< M, float_t >::operator!=( const vector< M, float_t >& other ) const
{
    return ! this->operator==( other );
}


template< size_t M, typename float_t >
bool
vector< M, float_t >::
isEqualTo( const vector< M, float_t >& other, float_t tolerance ) const
{
    bool ok = true;
    for( size_t index = 0; ok && index < M; ++index )
    {
        if ( at( index ) > other.at( index ) )
            ok = abs( at( index ) - other.at( index ) ) < tolerance;
        else
            ok = abs( other.at( index ) - at( index ) ) < tolerance;
    }
    return ok;

}



template< size_t M, typename float_t >
const vector< M, float_t >&
vector< M, float_t >::operator=( const float_t* c_array )
{
    copyFrom1DimCArray( c_array );
    return *this;
}



template< size_t M, typename float_t >
float_t
vector< M, float_t >::operator=( float_t filler_value )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) = filler_value;
    }
    return filler_value;
}




template< size_t M, typename float_t >
const vector< M, float_t >&
vector< M, float_t >::operator=( const vector< M, float_t >& other )
{
    memcpy( array, other.array, M * sizeof( float_t ) );
    return *this;
}



// returns void to avoid 'silent' loss of precision when chaining
template< size_t M, typename float_t >
template< typename other_float_t >
void
vector< M, float_t >::operator=( const vector< M, other_float_t >& other )
{
    for( size_t index = 0; index < M; ++index )
    {
        array[ index ] = static_cast< float_t >( other.array[ index ] );
    }
}



template< size_t M, typename float_t >
void
vector< M, float_t >::copyFrom1DimCArray( const float_t* c_array )
{
    memcpy( array, c_array, M * sizeof( float_t ) );
}



template< size_t M, typename float_t >
template< typename different_float_t >
void
vector< M, float_t >::copyFrom1DimCArray( const different_float_t* c_array )
{
    for( size_t index = 0; index < M; ++index )
    {
        at( index ) = static_cast< float_t >( c_array[ index ] );
    }
}



template< size_t M, typename float_t >
size_t
vector< M, float_t >::getM() const
{
    return M;
}



template< size_t M, typename float_t >
size_t
vector< M, float_t >::getSmallestComponentIndex() const
{
    size_t smallest_index = 0;
    for( size_t index = 1; index != M; ++index )
    {
        if ( array[ index ] < array[ smallest_index ] )
            smallest_index = index;
    }   
    return smallest_index;
}



template< size_t M, typename float_t >
size_t
vector< M, float_t >::getLargestComponentIndex() const
{
    size_t largest_index = 0;
    for( size_t index = 1; index != M; ++index )
    {
        if ( array[ index ] > array[ largest_index ] )
            largest_index = index;
    }   
    return largest_index;
}



template< size_t M, typename float_t >
float_t&
vector< M, float_t >::getSmallestComponent()
{
    return at( getSmallestComponentIndex() );
}



template< size_t M, typename float_t >
const float_t&
vector< M, float_t >::getSmallestComponent() const
{
    return at( getSmallestComponentIndex() );
}



template< size_t M, typename float_t >
float_t&
vector< M, float_t >::getLargestComponent()
{
    return at( getLargestComponentIndex() );
}



template< size_t M, typename float_t >
const float_t&
vector< M, float_t >::getLargestComponent() const
{
    return at( getLargestComponentIndex() );
}


} // namespace vmml

#endif

