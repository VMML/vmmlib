#ifndef __VMML__VECTOR__HPP__
#define __VMML__VECTOR__HPP__

#include <iostream>
#include <iomanip>
#include <vector>
#include <vmmlib/vmmlib_config.hpp>

namespace vmml
{

template< size_t M, typename float_t = double >
class vector
{
public:

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
vector< M, float_t >::at( size_t index )
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( index >= M )
    {
        throw "FIXME.";
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
        throw "FIXME.";
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

