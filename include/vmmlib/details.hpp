#ifndef __VMML__DETAILS__HPP__
#define __VMML__DETAILS__HPP__

#include <cmath>
namespace vmml
{

namespace details
{

// helper structs for SFINAE
template < typename T, typename U, bool b = T::value >
struct enable_if
{};


template < typename T, typename U >
struct enable_if <T, U, true > 
{
  typedef U type;
};


template< size_t M, size_t N >
struct is_square
{
	enum { value = M == N };
};


template< size_t M, size_t N >
struct is_4x4
{
	enum { value = ( M == 4 && N == 4 ) };
};



template< size_t M, size_t N, typename T >
void
matrix_is_square( typename enable_if< is_square< M, N >, T >::type* dummy = 0 )
{
    // intentionally left empty.
}



template< size_t M, size_t N, typename T >
void
matrix_is_4x4( typename enable_if< is_4x4< M, N >, T >::type* dummy = 0 )
{
    // intentionally left empty.
}


// helpers for certain cmath functions
template< typename float_t >
inline float_t
getSine( const float_t& angleInRadians )
{
    return sin( angleInRadians );
}

template<>
inline float
getSine( const float& angleInRadians )
{
    return sinf( angleInRadians );
}


template< typename float_t >
inline float_t
getCosine( const float_t& angleInRadians )
{
    return cos( angleInRadians );
}

template<>
inline float
getCosine( const float& angleInRadians )
{
    return cosf( angleInRadians );
}


} // namespace details

} // namespace vmml

#endif

