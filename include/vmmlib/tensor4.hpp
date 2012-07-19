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
		typedef float T_blas;
		
		
		
		//TODO: unfolding along all modes
		//TODO: accessors to tensor3 (along all modes)
		//TODO: maybe tensor4 iterator
		
		
		static const size_t ROWS	       = I1;
		static const size_t COLS	       = I2;
		static const size_t SLICES	       = I3;
		static const size_t T3S			   = I4;
		static const size_t MATRIX_SIZE    = I1 * I2;
		static const size_t T3_SIZE        = I1 * I2 * I3;
		static const size_t SIZE           = I1 * I2 * I3 * I4;
		
		
		
		// accessors
		inline T& operator()( size_t i1, size_t i2, size_t i3, size_t i4 );
		inline const T& operator()( size_t i1, size_t i2, size_t i3, size_t i4 ) const;
		
		inline T& at( size_t i1, size_t i2, size_t i3, size_t i4 );
		inline const T& at( size_t i1, size_t i2, size_t i3, size_t i4 ) const;
		
		
		// ctors
		tensor4();
        
		tensor4( const tensor4& source );
		
		template< typename U >
		tensor4( const tensor4< I1, I2, I3, I4, U >& source_ );
		
		template< size_t J1, size_t J2, size_t J3, size_t J4 >
		tensor4( const tensor4< J1, J2, J3, J4, T >& source_ );
		
		~tensor4();
		
		size_t size() const; // return I1 * I2 * I3 * I4;   
		
		// sets all elements to fill_value
		void operator=( T fill_value ); //@SUS: todo
		void fill( T fill_value ); //special case of set method (all values are set to the same value!)
		
		//sets all tensor values with random values
		//set srand(time(NULL)) or srand( seed )
		//if seed is set to -1, srand( seed ) was set outside set_random
		//otherwise srand( seed ) will be called with the given seed
		void fill_random( int seed = -1 );
		void fill_random_signed( int seed = -1 );
		void fill_increasing_values( );
		
		const tensor4& operator=( const tensor4& source_ );
		
		
		// note: this function copies elements until either the matrix is full or
		// the iterator equals end_.
		template< typename input_iterator_t >
		void set( input_iterator_t begin_, input_iterator_t end_, 
				 bool row_major_layout = true );	
		void zero();
		
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
		void cast_from( const tensor4< I1, I2, I3, I4, TT >& other );
		
		template< size_t J1, size_t J2, size_t  J3, typename TT >
		void cast_from( const tensor4< J1, J2, J3, I4, TT >& other, const long& slice_idx_start_ = 0 );
		
		
		template< typename TT >
		void float_t_to_uint_t( const tensor4< I1, I2, I3, I4, TT >& other );
		
		
		//TODO should move to other file. 
		void write_to_raw( const std::string& dir_, const std::string& filename_ ) const;
		void read_from_raw( const std::string& dir_, const std::string& filename_ ) ;
		void write_datfile( const std::string& dir_, const std::string& filename_ ) const;
		void write_to_csv( const std::string& dir_, const std::string& filename_ ) const;
		
	    
		bool operator==( const tensor4& other ) const;
		bool operator!=( const tensor4& other ) const;
		
		// due to limited precision, two 'idential' tensor4 might seem different.
		// this function allows to specify a tolerance when comparing matrices.
		bool equals( const tensor4& other, T tolerance ) const;
		// this version takes a comparison functor to compare the components of
		// the two tensor4 data structures
		template< typename compare_t >
		bool equals( const tensor4& other, compare_t& cmp ) const;
		

		inline tensor4 operator+( T scalar ) const;
		inline tensor4 operator-( T scalar ) const;
		
		void operator+=( T scalar );
		void operator-=( T scalar );
		
		inline tensor4 operator+( const tensor4& other ) const;
		inline tensor4 operator-( const tensor4& other ) const;
		
		template< size_t J1, size_t J2, size_t J3, size_t J4 >
		typename enable_if< J1 < I1 && J2 < I2 && J3 < I3 && J4 < I4 >::type*
		operator+=( const tensor4< J1, J2, J3, J4, T>& other );
		
		void operator+=( const tensor4& other );
		void operator-=( const tensor4& other );
		
		// 
		// tensor4-scalar operations / scaling 
		// 
		tensor4 operator*( T scalar );
		void operator*=( T scalar );
		
		tensor4 operator/( T scalar );
		void operator/=( T scalar );
		
		
		inline tensor4< I1, I2, I3, I4, T > operator-() const;
		tensor4< I1, I2, I3, I4, T > negate() const;
		
		friend std::ostream& operator << ( std::ostream& os, const tensor4< I1, I2, I3, I4, T >& t4 )
		{
			for(size_t i = 0; i < I4; ++i)
			{
				//TODO for tensor4
				//os << t3.get_frontal_slice_fwd( i ) << " *** " << std::endl;
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
		
		T*   _array;
		
}; // class tensor4
		
#define VMML_TEMPLATE_STRING    template< size_t I1, size_t I2, size_t I3, size_t I4, typename T >
#define VMML_TEMPLATE_CLASSNAME tensor4< I1, I2, I3, I4, T >
		


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
	delete[] array_;
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
		
				   
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
				   
} // namespace vmml
				   
#endif