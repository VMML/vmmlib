#ifndef __VMML__TENSOR3_DATA_STORE__HPP__
#define __VMML__TENSOR3_DATA_STORE__HPP__

#include <vmmlib/matrix.hpp>

namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, typename T >
	class tensor3_data_store
	{
	public:
		typedef matrix< I1, I2, T >		matrix_type;
		
		tensor3_data_store() 
		{
			for( size_t index = 0; index < I3; ++index )
			{
				_array[ index ] = new matrix_type( matrix_type::ZERO );
			}
		}
		
		tensor3_data_store( const tensor3_data_store& source_ ) 
		{
			for( size_t index = 0; index < I3; ++index )
			{
				_array[ index ] = new matrix_type( source_[ index ] );
			}
		}
		
		~tensor3_data_store()
		{
			for( size_t index = 0; index < I3; ++index )
			{
				delete _array[ index ];
			}
		}
		
		matrix_type& operator[]( size_t index )
		{
			assert( index < I3 );
			return *_array[ index % I3 ];
		}
		
		const matrix_type& operator[]( size_t index ) const
		{
			assert( index < I3 );
			return *_array[ index % I3 ];
		}
		
		void operator=( const tensor3_data_store& other_ )
		{
			for( size_t index = 0; index < I3; ++index )
			{
				*(_array[ index ]) = other_[ index ];
			}
		}
		
		
		
	protected:
		matrix< I1, I2, T >*	_array[I3];
		
	};
	
} // namespace vmml


#endif