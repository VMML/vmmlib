#ifndef __VMML__VECTOR__HPP__
#define __VMML__VECTOR__HPP__

#include <iostream>
#include <iomanip>
#include <vector>

#include <vmmlib/exception.hpp>
#include <vmmlib/vmmlib_config.hpp>

namespace vmml
{

template< size_t M, typename float_t = double >
class vector
{
public:

    inline float_t& operator()( size_t index );
    inline const float_t& operator()( size_t index ) const;
    inline float_t& at( size_t index );
    inline const float_t& at( size_t index ) const;

    inline float_t& operator[]( size_t index );
    inline const float_t& operator[]( size_t index ) const;

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
    
    inline float_t dot( const vector< M, float_t >& other ) const;
    static inline float_t dot( const vector< M, float_t >& v0, 
        const vector< M, float_t >& v1 );

    inline void normalize();
    
    inline float_t norm() const;
    inline float_t normSquared() const;
    
    void copyFrom1DimCArray( const float_t* c_array );

    template< typename different_float_t >
    void copyFrom1DimCArray( const different_float_t* c_array );
    
    float_t array[ M ];

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

}; // class vector


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
inline float_t&
vector< M, float_t >::x()
{
    return at( 0 );
}


template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::y()
{
    return at( 1 );
}


template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::z()
{
    return at( 2 );
}


template< size_t M, typename float_t >
inline float_t&
vector< M, float_t >::w()
{
    return at( 3 );
}


template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::x() const
{
    return at( 0 );
}


template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::y() const
{
    return at( 1 );
}


template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::z() const
{
    return at( 2 );
}


template< size_t M, typename float_t >
inline const float_t&
vector< M, float_t >::w() const
{
    return at( 3 );
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
inline float_t
vector< M, float_t >::dot( 
    const vector< M, float_t >& first, 
    const vector< M, float_t >& second ) 
{
    float_t tmp = 0.0;
    for( size_t index = 0; index < M; ++index )
    {
        tmp += first.at( index ) * second.at( index );
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
inline float_t
vector< M, float_t >::norm() const
{
    return sqrt( normSquared() );
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

} // namespace vmml

#endif

