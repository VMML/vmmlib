/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 * @author Jonas Boesch
 *
 * 
 */

#ifndef __VMML__TENSOR3_ITERATOR__HPP__
#define __VMML__TENSOR3_ITERATOR__HPP__

#include <iostream>
#include <assert.h>
#include <vmmlib/matrix.hpp>


namespace vmml
{
	

template< typename T >
class tensor3_iterator
{
public:
	
	typedef typename T::value_type		value_type;
	typedef typename T::pointer			pointer;
	typedef typename T::reference		reference;
	typedef typename matrix< T::ROWS, T::COLS, typename T::value_type >::iterator matrix_iterator;
	typedef typename matrix< T::ROWS, T::COLS, typename T::value_type >::const_iterator matrix_const_iterator;
	
	tensor3_iterator() : _tensor3( 0 ), _matrix_index( 0 ) {};
	
	tensor3_iterator( T& t_ ) : _tensor3( &t_ ), _matrix_index( 0 ) {};
	
	tensor3_iterator end() 
	{ 
		return tensor3_iterator( *this + T::SIZE - 1);
	};
	
	
	tensor3_iterator begin()
	{
		return tensor3_iterator( *this );
	}
	
	
	value_type operator*()
	{
		return *_matrix_it;
	}
	
	const value_type operator*() const
	{
		return *_matrix_it;
	}
	
	void operator++() 
	{
		assert( _tensor3 && "singular iterator"); // SHOULD BE EXCEPTION
		if ( _matrix_it == _matrix_it_end )
		{
			++_matrix_index;
			if ( _matrix_index == T::SLICES )
			{
				// don't do anything, we finished iterating
			}
			{
				_matrix_it = _tensor3.get_frontal_slice( _matrix_index ).begin();
				_matrix_it_end	= _tensor3.get_frontal_slice( _matrix_index ).end();
			}
		}
	}
	
	
protected:
	matrix_iterator		_matrix_it;
	matrix_iterator		_matrix_it_end;
	T*					_tensor3;
	size_t				_matrix_index;
	

}; //end tensor3_iterator class
}// end vmml namespace

#endif
