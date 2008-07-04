#ifndef __VMML__ARRAY__H__
#define __VMML__ARRAY__H__ 

#include <cstdlib>
#include <iterator>

/**
*   @brief stl-style wrapper for c-style arrays
*
*   @author jonas boesch
*/

namespace vmml
{

template < typename T, size_t N >
class array 
{
public:
	typedef T value_type;

	typedef value_type*         pointer;
	typedef const value_type*   const_pointer;

	typedef value_type&         reference;
	typedef const value_type&   const_reference;

	typedef ptrdiff_t           difference_type;
	typedef size_t              size_type;

	typedef pointer             iterator;
	typedef const_pointer       const_iterator;

    typedef typename std::reverse_iterator< iterator > 
        reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator>
        const_reverse_iterator;	

    array(); // uninitialized ctor
    array( const T& element );

	inline iterator begin();
	inline iterator end();

	inline const_iterator begin() const;
	inline const_iterator end() const;

	inline reverse_iterator rbegin();
	inline reverse_iterator rend();

	inline const_reverse_iterator rbegin() const;
	inline const_reverse_iterator rend() const;

	inline reference operator[] ( size_type index );
	inline const_reference operator[] ( size_type index ) const;

	inline size_type size() const;

protected:
	T _data[ N ];
};


template< typename T, size_t N >
array< T, N >::array()
{
    // NO INITIALISATION
}


template< typename T, size_t N >
array< T, N >::array( const T& element )
{
    for ( size_t i = 0; i < N; ++i )
    {
        memcpy( &_data[ i ], &element, sizeof( T ) ); 
    }
}



template< typename T, size_t N >
inline typename array< T, N >::iterator
array< T, N >::begin()
{ 
	return _data; 
}



template< typename T, size_t N >
inline typename array< T, N >::iterator 
array< T, N >::end()
{ 
	return _data + N ;
}



template< typename T, size_t N >
inline typename array< T, N >::const_iterator 
array< T, N >::begin() const
{ 
	return _data;
}



template< typename T, size_t N >
inline typename array< T, N >::const_iterator 
array< T, N >::end() const
{
	return _data + N; 
}



template< typename T, size_t N >
inline typename array< T, N >::reverse_iterator
array< T, N >::rbegin()
{ 
	return reverse_iterator( &_data[ N ] ); 
}



template< typename T, size_t N >
inline typename array< T, N >::reverse_iterator 
array< T, N >::rend()
{ 
	return reverse_iterator( _data );
}



template< typename T, size_t N >
inline typename array< T, N >::const_reverse_iterator 
array< T, N >::rbegin() const
{ 
	return const_reverse_iterator( &_data[ N ] );
}



template< typename T, size_t N >
inline typename array< T, N >::const_reverse_iterator 
array< T, N >::rend() const
{
	return const_reverse_iterator( _data ); 
}



template< typename T, size_t N >
inline typename array< T, N >::reference 
array< T, N >::operator[]( size_type n )
{ 
	return _data[ n ];
}



template< typename T, size_t N >
inline typename array< T, N >::const_reference 
array< T, N >::operator[]( size_type n ) const 
{ 
	return _data[ n ]; 
}



template< typename T, size_t N >
inline typename array< T, N >::size_type 
array< T, N >::size() const
{ 
	return N; 
}


} // namespace vmmltl

#endif

