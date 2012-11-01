/*
 * VMMLib - Tensor Classes
 *
 * @author Susanne Suter
 * @author Jonas Boesch
 * @author David Klaper
 *
 * a tensor is a generalization of a multidimensional array
 * a tensor4 is a tensor data structure with four modes I1, I2, I3 and I4
 */

#ifndef __VMML__TENSOR4__HPP__
#define __VMML__TENSOR4__HPP__

#include <fstream>   // file I/O
//#include <vmmlib/tensor4_iterator.hpp>
#include <vmmlib/enable_if.hpp>
//#include <vmmlib/blas_dot.hpp>
//#include <sys/mman.h>
//#include <fcntl.h>
//#include <omp.h>
//#include <limits>
//#undef min
//#undef max
#include "tensor3.hpp"


namespace vmml
{
	
	// tensor with four modes, containing a series I4 x tensor4 as (I3 of I1 x I2 vmml matrices)
	//I1 is number of rows, I2 is number of columns and I3 is number of tubes, I4 is number of tensor4s
	
	template< size_t I1, size_t I2, size_t I3, size_t I4, typename T = float >
	class tensor4
	{
	public:
		typedef T                                       value_type;
		typedef T*                                      pointer;
		typedef T&                                      reference;
		
		typedef tensor3< I1, I2, I3, T> tensor3_t;
		
		//TODO: maybe tensor4 iterator
		//TODO: unfolding along all modes
		//TODO: accessors to tensor3 (along all modes)
        
        // Average all values along 1 axis
		tensor3_t& average_I4(tensor3_t& t3) const;
        
		inline void get_tensor3( const size_t i4_, tensor3_t& t3_data_ ) const; //TODO: DK done
		
		inline tensor3_t& get_tensor3( size_t i4_ );
		inline const tensor3_t& get_tensor3( size_t i4_ ) const;
		
		static const size_t ROWS	       = I1;
		static const size_t COLS	       = I2;
		static const size_t SLICES	       = I3;
		static const size_t T3S			   = I4;
		static const size_t MATRIX_SIZE    = I1 * I2;
		static const size_t T3_SIZE        = I1 * I2 * I3;
		static const size_t SIZE           = I1 * I2 * I3 * I4;
		
		static size_t get_array_size_in_bytes();
		
	    // WARNING: dangerous. Use before destruction if you want to prevent
		// a delete call for the assigned T* _array in the destructor.
		void clear_array_pointer();
		
		// accessors
		inline T& operator()( size_t i1, size_t i2, size_t i3, size_t i4 ); //TODO: DK done
		inline const T& operator()( size_t i1, size_t i2, size_t i3, size_t i4 ) const;  //TODO: DK done
		
		inline T& at( size_t i1, size_t i2, size_t i3, size_t i4 );  //TODO: DK done
		inline const T& at( size_t i1, size_t i2, size_t i3, size_t i4 ) const;  //TODO: DK done
		
		
		// ctors
		tensor4();
		
		explicit tensor4(void* memory) ;
		
		tensor4( const tensor4& source );  //TODO: DK done
		
		template< typename U >
		tensor4( const tensor4< I1, I2, I3, I4, U >& source_ );  //TODO: DK done
		
		template< size_t J1, size_t J2, size_t J3, size_t J4 >
		tensor4( const tensor4< J1, J2, J3, J4, T >& source_ );  //TODO: DK done
		
		~tensor4();
		
		size_t size() const; // return I1 * I2 * I3 * I4;    //TODO: DK done
		
		// sets all elements to fill_value
		void operator=( T fill_value ); //@SUS: todo  //TODO: DK done
		void fill( T fill_value ); //special case of set method (all values are set to the same value!)  //TODO: DK done
		
		//sets all tensor values with random values
		//if seed is negative, srand( seed ) should have been set outside fill_random
		//if seed is 0 or greater srand( seed ) will be called with the given seed
		void fill_random( int seed = -1 );  //TODO: DK done
		void fill_random_signed( int seed = -1 ); //TODO: DK done
		void fill_increasing_values( ); //TODO: DK done
		
		const tensor4& operator=( const tensor4& source_ ); //TODO: DK done
		
		
		// note: this function copies elements until either the matrix is full or
		// the iterator equals end_.
		template< typename input_iterator_t >
		void set( input_iterator_t begin_, input_iterator_t end_,
				 bool row_major_layout = true );
		void zero(); //TODO: DK done
		
		T get_min() const;
		T get_max() const;
		T get_abs_min() const;
		T get_abs_max() const;
		
		//returns number of non-zeros
		size_t nnz() const;
		size_t nnz( const T& threshold_ ) const;
		void threshold( const T& threshold_value_ );
		
		//error computation
		double frobenius_norm() const;
		double frobenius_norm( const tensor4< I1, I2, I3, I4, T >& other ) const;
		double avg_frobenius_norm() const;
		double mse( const tensor4< I1, I2, I3, I4, T >& other ) const; // mean-squared error
		double rmse( const tensor4< I1, I2, I3, I4, T >& other ) const; //root mean-squared error
		double compute_psnr( const tensor4< I1, I2, I3, I4, T >& other, const T& max_value_ ) const; //peak signal-to-noise ratio
		void mean( T& mean_ ) const;
		double mean() const;
		double variance() const;
		double stdev() const;
		
		template< typename TT >
		void cast_from( const tensor4< I1, I2, I3, I4, TT >& other ); //TODO: DK (2) done
		
		template< size_t J1, size_t J2, size_t  J3, size_t J4, typename TT >
		void cast_from( const tensor4< J1, J2, J3, J4, TT >& other, const long& slice_idx_start_ = 0 ); //TODO: DK (2) done
		
		
		template< typename TT >
		void float_t_to_uint_t( const tensor4< I1, I2, I3, I4, TT >& other ); //TODO: DK (2) done
		
	    //check if corresponding tensor values are equal or not
		bool operator==( const tensor4& other ) const; //TODO: DK done
		bool operator!=( const tensor4& other ) const; //TODO: DK done
		
		// due to limited precision, two 'idential' tensor4 might seem different.
		// this function allows to specify a tolerance when comparing matrices.
		bool equals( const tensor4& other, T tolerance ) const;
		// this version takes a comparison functor to compare the components of
		// the two tensor4 data structures
		template< typename compare_t >
		bool equals( const tensor4& other, compare_t& cmp ) const;
		
		
		inline tensor4 operator+( T scalar ) const; //TODO: DK (2) done
		inline tensor4 operator-( T scalar ) const; //TODO: DK (2) done
		
		void operator+=( T scalar ); //TODO: DK (2) done
		void operator-=( T scalar ); //TODO: DK (2) done
		
		inline tensor4 operator+( const tensor4& other ) const; //TODO: DK (2) done
		inline tensor4 operator-( const tensor4& other ) const; //TODO: DK (2) done
		
		template< size_t J1, size_t J2, size_t J3, size_t J4 >
		typename enable_if< J1 < I1 && J2 < I2 && J3 < I3 && J4 < I4 >::type*
		operator+=( const tensor4< J1, J2, J3, J4, T>& other ); //TODO: DK (2) done
		
		void operator+=( const tensor4& other );
		void operator-=( const tensor4& other );
		
		//
		// tensor4-scalar operations / scaling
		//
		tensor4 operator*( T scalar ); //TODO: DK (2) done
		void operator*=( T scalar ); //TODO: DK (2) done
		
		tensor4 operator/( T scalar ); //TODO: DK (2) done
		void operator/=( T scalar ); //TODO: DK (2) done
		
		
		inline tensor4< I1, I2, I3, I4, T > operator-() const; //TODO: DK (2) done
		tensor4< I1, I2, I3, I4, T > negate() const; //TODO: DK (2) done
		
		friend std::ostream& operator << ( std::ostream& os, const tensor4< I1, I2, I3, I4, T >& t4 ) //TODO: DK (tensor3) done
		{
			//FIXME: to this directly with tensors
			//sth like:
			//os << t4.get_tensor3( i ) << " xxx " << std::endl;
			for(size_t i4 = 0; i4 < I4; ++i4)
			{
				for(size_t i3 = 0; i3 < I3; ++i3)
				{
					for(size_t i1 = 0; i1 < I1; ++i1)
					{
						os << "(";
						
						for(size_t i2 = 0; i2 < I2; ++i2)
						{
							if (i2 == (I2 - 1) )
							{
								os << t4.at( i1, i2, i3, i4 ) ;
							} else {
								os << t4.at( i1, i2, i3, i4 ) << ", " ;
							}
							
						}
						os << ")" << std::endl;
					}
					if ( i3 < (I3 -1 ))
					{
						os << " xxx " << std::endl;
					}
				}
				if ( i4 + 1 < I4 )
				{
					os << "---- " << std::endl;
				}
			}
			return os;
		}
		
		
		// static members
		static void     tensor4_allocate_data( T*& array_ );
		static void     tensor4_deallocate_data( T*& array_ );
		
		static const tensor4< I1, I2, I3, I4, T > ZERO;
		
		T*          get_array_ptr();
		const T*    get_array_ptr() const;
		
		// computes the array index for direct access
		inline size_t compute_index( size_t i1, size_t i2, size_t i3, size_t i4 ) const;
		
		
	protected:
		tensor3_t&                   _get_tensor3( size_t index_ );
		const tensor3_t&             _get_tensor3( size_t index_ ) const;
		
		T*   _array;
		
		}; // class tensor4
		
#define VMML_TEMPLATE_STRING    template< size_t I1, size_t I2, size_t I3, size_t I4, typename T >
#define VMML_TEMPLATE_CLASSNAME tensor4< I1, I2, I3, I4, T >
		
		
		
		// WARNING: make sure that memory is a pointer to a memory block of 
		// sufficient size (that is, is at least get_array_size_in_bytes())
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME::tensor4( void* memory ) 
		: _array( reinterpret_cast<T*>( memory ) )
		{
			assert( _array );
		}	
		
		
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME::tensor4()
		: _array()
		{
			tensor4_allocate_data( _array );
		}
		
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME::~tensor4()
		{
			tensor4_deallocate_data( _array );
		}
		
		
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::
		tensor4_allocate_data( T*& array_ )
		{
			array_ = new T[ I1 * I2 * I3 * I4];
		}
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::
		tensor4_deallocate_data( T*& array_ )
		{
			if (array_)
			{
				delete[] array_;
			}
		}
		
		
		
		VMML_TEMPLATE_STRING
		T*
		VMML_TEMPLATE_CLASSNAME::get_array_ptr()
		{
			return _array;
		}
		
		
		
		VMML_TEMPLATE_STRING
		const T*
		VMML_TEMPLATE_CLASSNAME::get_array_ptr() const
		{
			return _array;
		}
		
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME::tensor4( const tensor4& source_ )
		{
			tensor4_allocate_data( _array );
			(*this) = source_;
			
		}
		
		VMML_TEMPLATE_STRING
		template< typename U >
		VMML_TEMPLATE_CLASSNAME::tensor4( const tensor4< I1, I2, I3, I4, U >& source_ )
		{
			const U* s_array = source_.get_array_ptr();
			tensor4_allocate_data( _array );
			for (size_t index = 0; index < SIZE; ++index)
			{
				_array[ index ] = static_cast< T >( s_array[ index ] );
			}
		}
		
		VMML_TEMPLATE_STRING
		template< size_t J1, size_t J2, size_t J3, size_t J4 >
		VMML_TEMPLATE_CLASSNAME::tensor4( const tensor4< J1, J2, J3, J4, T >& source_ )
		{
			tensor4_allocate_data( _array );
			size_t min1 = J1 < I1 ? J1 : I1;
			size_t min2 = J2 < I2 ? J2 : I2;
			size_t min3 = J3 < I3 ? J3 : I3;
			size_t min4 = J4 < I4 ? J4 : I4;
			
			zero();
			
			for(size_t i4 = 0; i4 < min4; ++i4 )
			{
				for(size_t i3 = 0; i3 < min3; ++i3)
				{
					for(size_t i2 = 0; i2 < min2; ++i2)
					{
						for(size_t i1 = 0; i1 < min1; ++i1)
						{
							at(i1, i2, i3, i4) = source_(i1, i2, i3, i4);
						}
					}
				}
			}
			
		} 
		
		
		
		VMML_TEMPLATE_STRING
		size_t 
		VMML_TEMPLATE_CLASSNAME::size() const
		{
			return SIZE;
		}
		
		VMML_TEMPLATE_STRING
		bool 
		VMML_TEMPLATE_CLASSNAME::operator==( const tensor4< I1, I2, I3, I4, T >& other ) const
		{
			const T* other_array = other.get_array_ptr();
			for(size_t index = 0; index < SIZE; ++index)
			{
				if(_array[ index ] != other_array[ index] )
				{
					return false;
				}
			}
			
			return true;
		}
		
		
		VMML_TEMPLATE_STRING
		bool 
		VMML_TEMPLATE_CLASSNAME::operator!=( const tensor4< I1, I2, I3, I4, T >& other ) const
		{
			return ! operator==( other );
		}	
		
		
		VMML_TEMPLATE_STRING
		typename VMML_TEMPLATE_CLASSNAME::tensor3_t&
		VMML_TEMPLATE_CLASSNAME::
		_get_tensor3( size_t index_ )
		{
			tensor3<I1, I2, I3, T>* tens3 = new tensor3<I1, I2, I3, T>( (void *)(_array + I1 * I2 * I3 * index_));
			return *tens3;
		}
		
		
		VMML_TEMPLATE_STRING
		const typename VMML_TEMPLATE_CLASSNAME::tensor3_t&
		VMML_TEMPLATE_CLASSNAME::
		_get_tensor3( size_t index_ ) const
		{
			const tensor3<I1, I2, I3, T>* tens3 = new tensor3<I1, I2, I3, T>( (void *)(_array + I1 * I2 * I3 * index_));
			return *tens3;
		}
		
		
		VMML_TEMPLATE_STRING
		inline void
		VMML_TEMPLATE_CLASSNAME::
		get_tensor3( const size_t i4_, tensor3_t& t3_data_ ) const
		{
#ifdef VMMLIB_SAFE_ACCESSORS
			if ( i4_ >= I4 )
				VMMLIB_ERROR( "get_tensor3() - index out of bounds.", VMMLIB_HERE );
#endif
			
			t3_data_ = _get_tensor3( i4_ );
		}
		
		
		
		VMML_TEMPLATE_STRING
		inline typename VMML_TEMPLATE_CLASSNAME::tensor3_t&
		VMML_TEMPLATE_CLASSNAME::
		get_tensor3( size_t i4_ )
		{
#ifdef VMMLIB_SAFE_ACCESSORS
			if ( i4_ >= I4 )
				VMMLIB_ERROR( "get_tensor3() - index out of bounds.", VMMLIB_HERE );
#endif
			return _get_tensor3( i4_ );
		}
		
		
		VMML_TEMPLATE_STRING
		inline const typename VMML_TEMPLATE_CLASSNAME::tensor3_t&
		VMML_TEMPLATE_CLASSNAME::
		get_tensor3( size_t i4_ ) const
		{
#ifdef VMMLIB_SAFE_ACCESSORS
			if ( i4_ >= I4 )
				VMMLIB_ERROR( "get_tensor3() - index out of bounds.", VMMLIB_HERE );
#endif
			return _get_tensor3( i4_ );
		}
		
		
		//fill
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::
		fill( T fillValue )
		{
			for( size_t i4 = 0; i4 < I4; ++i4 )
			{
				for( size_t i3 = 0; i3 < I3; ++i3 )
				{
					for( size_t i1 = 0; i1 < I1; ++i1 )
					{
						for( size_t i2 = 0; i2 < I2; ++i2 )
						{
							at(i1, i2, i3, i4) = fillValue;
						}
					}
				}
			}
			//FIXME:
			//_get_tensor3( i4 ).fill( fillValue );
			
		}
		
		
		VMML_TEMPLATE_STRING
		inline T&
		VMML_TEMPLATE_CLASSNAME::at( size_t i1, size_t i2, size_t i3, size_t i4 )
		{
#ifdef VMMLIB_SAFE_ACCESSORS
			if ( i1 >= I1 || i2 >= I2 || i3 >= I3 || i4 >= I4 )
				VMMLIB_ERROR( "at( i1, i2, i3, i4 ) - index out of bounds", VMMLIB_HERE );
#endif
			return _array[ i4 * T3_SIZE + i3 * MATRIX_SIZE + i2 * ROWS + i1 ];
			
		}
		
		
		
		VMML_TEMPLATE_STRING
		const inline T&
		VMML_TEMPLATE_CLASSNAME::at( size_t i1, size_t i2, size_t i3, size_t i4 ) const
		{
#ifdef VMMLIB_SAFE_ACCESSORS
			if ( i1 >= I1 || i2 >= I2 || i3 >= I3 || i4 >= I4 )
				VMMLIB_ERROR( "at( i1, i2, i3, i4 ) - index out of bounds", VMMLIB_HERE );
#endif
			return _array[ i4 * T3_SIZE + i3 * MATRIX_SIZE + i2 * ROWS + i1 ];
		}
		
		
		
		VMML_TEMPLATE_STRING
		inline T&
		VMML_TEMPLATE_CLASSNAME::operator()( size_t i1, size_t i2, size_t i3, size_t i4 )
		{
			return at( i1, i2, i3, i4 );
		}
		
		
		
		VMML_TEMPLATE_STRING
		const inline T&
		VMML_TEMPLATE_CLASSNAME::operator()(  size_t i1, size_t i2, size_t i3, size_t i4 ) const
		{     
			return at( i1, i2, i3, i4 );
		}
		
		/*VMML_TEMPLATE_STRING
		 void 
		 VMML_TEMPLATE_CLASSNAME::fill( T fill_value )
		 {
		 for(size_t index = 0; index < SIZE; ++index)
		 {
		 _array[ index ] = fill_value;
		 }
		 }*/
		
		VMML_TEMPLATE_STRING
		void 
		VMML_TEMPLATE_CLASSNAME::operator=( T fill_value )
		{
			fill( fill_value );
		}
		
		VMML_TEMPLATE_STRING
		void 
		VMML_TEMPLATE_CLASSNAME::zero()
		{
			fill( static_cast< T >( 0.0 ) );
		}
		
		VMML_TEMPLATE_STRING
		void 
		VMML_TEMPLATE_CLASSNAME::fill_increasing_values()
		{
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[ index ] = index;
			}
		}
		
		VMML_TEMPLATE_STRING
		void 
		VMML_TEMPLATE_CLASSNAME::fill_random(int seed)
		{
			if ( seed >= 0 )
				srand( seed );
			
			double fillValue = 0.0f;
			for( size_t index = 0; index < SIZE; ++index )
			{
				fillValue = rand();
				fillValue /= RAND_MAX;
				fillValue *= std::numeric_limits< T >::max();
				_array[ index ] = static_cast< T >( fillValue )  ;
			}
		}
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::fill_random_signed(int seed)
		{
			if ( seed >= 0 )
				srand( seed );
			
			double fillValue = 0.0f;
			for( size_t index = 0; index < SIZE; ++index )
			{
				
				fillValue = rand();
				fillValue /= RAND_MAX;
				fillValue *= std::numeric_limits< T >::max();
				T fillValue2 = static_cast< T >(fillValue) % std::numeric_limits< T >::max();
				fillValue2 -= std::numeric_limits< T >::max()/2;
				_array[ index ] = fillValue2  ;
				// test if ever > max/2 or < -max/2
				
			}
		}
		
		VMML_TEMPLATE_STRING
		const VMML_TEMPLATE_CLASSNAME&
		VMML_TEMPLATE_CLASSNAME::operator=( const VMML_TEMPLATE_CLASSNAME& source_ )
		{
			if(this != &source_) // avoid self assignment
			{
				memcpy( _array, source_._array, I1 * I2 * I3 * I4 * sizeof( T ) );
			}
			return *this;
		}
		
		VMML_TEMPLATE_STRING
		template< typename TT >
		void 
		VMML_TEMPLATE_CLASSNAME::cast_from( const tensor4< I1, I2, I3, I4, TT >& other ) 
		{
#pragma omp parallel for
			for(long tensor_idx = 0; tensor_idx < (long)I4; ++tensor_idx)
			{
#pragma omp parallel for
				for( long slice_idx = 0; slice_idx < (long)I3; ++slice_idx )
				{
#pragma omp parallel for
					for( long row_index = 0; row_index < (long)I1; ++row_index )
					{
#pragma omp parallel for
						for( long col_index = 0; col_index < (long)I2; ++col_index )
						{
							at( row_index, col_index, slice_idx, tensor_idx ) =  static_cast< T >(other.at( row_index, col_index, slice_idx, tensor_idx ));
						}
					}
				}
			}
			
		}
		
		VMML_TEMPLATE_STRING
		template< size_t J1, size_t J2, size_t J3, size_t J4, typename TT >
		void
		VMML_TEMPLATE_CLASSNAME::cast_from( const tensor4< J1, J2, J3, J4, TT >& other, const long& slice_idx_start_ )
		{
			this->zero();
			int j4 = I4 < J4? I4 : J4;
			int j3 = I3 < J3? I3 : J3;
			int j2 = I2 < J2? I2 : J2;
			int j1 = I1 < J1? I1 : J1;
#pragma omp parallel for
			for(long tensor_idx = 0; tensor_idx < (long)j4; ++tensor_idx)
			{
#pragma omp parallel for
				for( long slice_idx = 0; slice_idx < (long)j3; ++slice_idx )
				{
#pragma omp parallel for
					for( long row_index = 0; row_index < (long)j1; ++row_index )
					{
#pragma omp parallel for
						for( long col_index = 0; col_index < (long)j2; ++col_index )
						{
							at( row_index, col_index, slice_idx, tensor_idx ) =  static_cast< T >(other.at( row_index, col_index, slice_idx, tensor_idx ));
						}
					}
				}
			}
			
		}
		
		VMML_TEMPLATE_STRING
		template< typename TT >
		void
		VMML_TEMPLATE_CLASSNAME::float_t_to_uint_t( const tensor4< I1, I2, I3, I4, TT >& other )
		{
			typedef tensor4< I1, I2, I3, I4, TT > t3_tt_type ;
			
			if( sizeof(T) == 1 || sizeof(T) == 2) {
				const TT* otherdata = other.get_array_ptr();
				for( size_t index = 0;index < SIZE; ++index )
				{
					_array[index] = T( std::min( std::max(int(0), int( otherdata[index] + 0.5)), int(std::numeric_limits< T >::max()) ));
				}
			}
			else {
				this->cast_from( other );
				return;
			}
		}
		
		VMML_TEMPLATE_STRING
		inline VMML_TEMPLATE_CLASSNAME
		VMML_TEMPLATE_CLASSNAME::operator+( T scalar ) const
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			result += scalar;
			return result;
			
		}
		
		VMML_TEMPLATE_STRING
		inline VMML_TEMPLATE_CLASSNAME
		VMML_TEMPLATE_CLASSNAME::operator-( T scalar ) const
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			result -= scalar;
			return result;
			
		}
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::operator+=( T scalar )
		{
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[index] += scalar;
			}
			
		}
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::operator-=( T scalar )
		{
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[index] -= scalar;
			}
			
		}
		
		VMML_TEMPLATE_STRING
		inline VMML_TEMPLATE_CLASSNAME
		VMML_TEMPLATE_CLASSNAME::operator+( const tensor4& other ) const
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			result += other;
			return result;
		}
		
		
		VMML_TEMPLATE_STRING
		inline VMML_TEMPLATE_CLASSNAME
		VMML_TEMPLATE_CLASSNAME::operator-( const tensor4& other ) const
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			result -= other;
			return result;
		}
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::operator+=( const tensor4& other )
		{
			const T* dataptr = other.get_array_ptr();
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[index] += dataptr[index];
			}
		}
		
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::operator-=( const tensor4& other )
		{
			const T* dataptr = other.get_array_ptr();
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[index] -= dataptr[index];
			}
		}
		
		VMML_TEMPLATE_STRING
		template< size_t J1, size_t J2, size_t J3, size_t J4 >
		typename enable_if< J1 < I1 && J2 < I2 && J3 < I3 && J4 < I4 >::type*
		VMML_TEMPLATE_CLASSNAME::operator+=( const tensor4< J1, J2, J3, J4, T>& other )
		{
			for(size_t i4 = 0; i4 < J4; ++i4)
			{
				for( size_t i3 = 0; i3 < J3; ++i3 )
				{
					for( size_t i1 = 0; i1 < J1; ++i1 )
					{
						for( size_t i2 = 0; i2 < J2; ++i2 )
						{
							at( i1, i2, i3, i4 ) += other.at(i1, i2, i3, i4);
						}
					}
				}
			}
			return 0;
		}
		
		
		
		//
		// tensor4-scalar operations / scaling
		//
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME 
		VMML_TEMPLATE_CLASSNAME::operator*( T scalar )
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			result *= scalar;
			return result;
		}
		
		VMML_TEMPLATE_STRING
		void 
		VMML_TEMPLATE_CLASSNAME::operator*=( T scalar )
		{
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[index] *= scalar;
			}
		}
		
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME 
		VMML_TEMPLATE_CLASSNAME::operator/( T scalar )
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			result /= scalar;
			return result;
		}
		
		VMML_TEMPLATE_STRING
		void 
		VMML_TEMPLATE_CLASSNAME::operator/=( T scalar )
		{
			for(size_t index = 0; index < SIZE; ++index)
			{
				_array[index] /= scalar;
			}
		}
		
		
		VMML_TEMPLATE_STRING
		inline VMML_TEMPLATE_CLASSNAME 
		VMML_TEMPLATE_CLASSNAME::operator-() const
		{
			return negate();   
		}
		
		VMML_TEMPLATE_STRING
		VMML_TEMPLATE_CLASSNAME 
		VMML_TEMPLATE_CLASSNAME::negate() const
		{
			vmml::tensor4<I1, I2, I3, I4, T> result(*this);
			T* dataptr = result.get_array_ptr();
			
			for(size_t index = 0; index < SIZE; ++index)
			{
				dataptr[index] = _array[index] * -1;
			}
			return result;
			
		}
		
		VMML_TEMPLATE_STRING
		typename VMML_TEMPLATE_CLASSNAME::tensor3_t&
		VMML_TEMPLATE_CLASSNAME::average_I4(tensor3_t& t3) const
		{
			
			for(size_t i3 = 0; i3 < I3; ++i3)
			{
				for( size_t i1 = 0; i1 < I1; ++i1 )
				{
					for( size_t i2 = 0; i2 < I2; ++i2 )
					{
						T sum = 0;
						for( size_t i4 = 0; i4 < I4; ++i4 )
						{
							sum += at(i1,i2,i3,i4);
						}
						t3.at(i1,i2,i3) = sum/I4;
					}
				}
			}
			return t3;
			
		}
		
		VMML_TEMPLATE_STRING
		size_t
		VMML_TEMPLATE_CLASSNAME::get_array_size_in_bytes()
		{
			return (sizeof(T) * SIZE);
		}
		
		
		
		// WARNING: dangerous. Use before destruction if you want to prevent
		// a delete call for the assigned T* _array in the destructor.
		VMML_TEMPLATE_STRING
		void
		VMML_TEMPLATE_CLASSNAME::clear_array_pointer() 
		{
			_array = 0;
		}
		
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
		
		} // namespace vmml
		
#endif
