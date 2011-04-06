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
        
        static const size_t matrix_size = I1 * I2;
        
		tensor3_data_store() 
		{
            _allocate_memory();
		}
		
		tensor3_data_store( const tensor3_data_store& source_ ) 
		{
            _allocate_memory();
            (*this) = source_;
		}
		
		~tensor3_data_store()
		{
            _deallocate_memory();
		}
		
		inline matrix_type& operator[]( size_t index )
		{
			assert( index < I3 );
            return *reinterpret_cast< matrix_type* >( array + matrix_size * index );
		}
		
		inline const matrix_type& operator[]( size_t index ) const
		{
			assert( index < I3 );
            return *reinterpret_cast< const matrix_type* >( array + matrix_size * index );
		}
		
		inline void operator=( const tensor3_data_store& other_ )
		{
            memcpy( array, other_.array, I1 * I2 * I3 * sizeof( T ) );
		}
        
        inline void _allocate_memory()
        {
            array = new T[ I1 * I2 * I3];
        }
        
        inline void _deallocate_memory()
        {
            delete[] array;
        }
		
		
		
	protected:
		T*          array;
		
	};
	
} // namespace vmml


#endif