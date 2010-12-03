/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 * @author Jonas Boesch
 * 
 * a tensor is a generalization of a multidimensional array
 * a tensor3 is a tensor data structure with three modes I1, I2 and I3
 */
 
#ifndef __VMML__TENSOR3__HPP__
#define __VMML__TENSOR3__HPP__

#include <fstream>   // file I/O
#include <vmmlib/tensor3_iterator.hpp>
#include <vmmlib/tensor3_data_store.hpp>
#include <vmmlib/enable_if.hpp>

namespace vmml
{
	
// tensor with three modes, containing a series I3 of I1 x I2 vmml matrices
//I1 is number of rows, I2 is number of columns and I3 is number of tubes
template< size_t I1, size_t I2, size_t I3, typename T = float >
class tensor3
{
public:
    typedef T                                       value_type;
    typedef T*                                      pointer;
    typedef T&                                      reference;
    
    typedef typename matrix< I1, I2, T>::iterator                                matrix_iterator;
		
    typedef typename vmml::tensor3_iterator< tensor3< I1, I2, I3, T > >	         iterator;
    typedef typename vmml::tensor3_const_iterator< tensor3< I1, I2, I3, T > >    const_iterator;
	
    typedef typename vmml::tensor3_iterator< tensor3< I1, I2, I3, T > >          reverse_iterator;
    typedef typename vmml::tensor3_iterator< tensor3< I1, I2, I3, T > >          const_reverse_iterator;
    
    typedef tensor3_data_store< I1, I2, I3, T >                                  store_type;
	
    typedef matrix< I1, I2, T >        front_slice_type; //fwd: forward cylcling (after kiers et al., 2000)
    typedef matrix< I3, I1, T >        lat_slice_type;
    typedef matrix< I2, I3, T >        horiz_slice_type;

    typedef matrix< I2, I1, T >        bwd_front_slice_type; //bwd: backward cylcling (after lathauwer et al., 2000a)
    typedef matrix< I1, I3, T >        bwd_lat_slice_type;
    typedef matrix< I3, I2, T >        bwd_horiz_slice_type;

    typedef matrix< I1, I2*I3, T >     bwd_lat_unfolding_type;
    typedef matrix< I2, I1*I3, T >     bwd_front_unfolding_type;
    typedef matrix< I3, I1*I2, T >     bwd_horiz_unfolding_type;
    
    static const size_t ROWS	       = I1;
    static const size_t COLS	       = I2;
    static const size_t SLICES	       = I3;
    static const size_t SIZE           = I1 * I2 * I3;
 
    // accessors
    inline T& operator()( size_t i1, size_t i2, size_t i3 );
    inline const T& operator()( size_t i1, size_t i2, size_t i3 ) const;

    inline T& at( size_t i1, size_t i2, size_t i3 );
    inline const T& at( size_t i1, size_t i2, size_t i3 ) const;
    
    // element iterators - NOTE: column-major order
    iterator                begin();
    iterator                end();
    
    const_iterator          begin() const;
    const_iterator          end() const;
	
#if 0
    reverse_iterator        rbegin();
    reverse_iterator        rend();
    const_reverse_iterator  rbegin() const;
    const_reverse_iterator  rend() const;
#endif
    
    // ctors
    tensor3();
        
    tensor3( const tensor3& source );
    
    template< typename U >
    tensor3( const tensor3< I1, I2, I3, U >& source_ );
	
    template< size_t J1, size_t J2, size_t J3>
    tensor3( const tensor3< J1, J2, J3, T >& source_ );
	
    size_t size() const; // return I1 * I2 * I3;   
	
	template< size_t J1, size_t J2, size_t J3 >
	tensor3<J1, J2, J3, T>
    get_sub_tensor3( size_t row_offset, size_t col_offset, size_t slice_offset = 0,
					typename enable_if< J1 <= I1 && J2 <= I2 && J3 <= I3 >::type* = 0 ) const;
	
	template< size_t J1, size_t J2, size_t J3 >
	typename enable_if< J1 <= I1 && J2 <= I2 && J3 <= I3 >::type*
    get_sub_tensor3( tensor3<J1, J2, J3, T >& result, 
					size_t row_offset = 0, size_t col_offset = 0, size_t slice_offset = 0 ) const;

    inline void get_I1_vector( size_t i2, size_t i3, vector< I1, T >& data ) const; // I1_vector is a column vector with all values i1 at i2 and i3
    inline void get_I2_vector( size_t i1, size_t i3, vector< I2, T >& data ) const; // I2_vector is a row vector with all values i2 at i1 and i3
    inline void get_I3_vector( size_t i1, size_t i2, vector< I3, T >& data ) const; // I3_vector is a vector with all values i3 at a given i1 and i2

    inline void get_row( size_t i1, size_t i3, vector< I2, T >& data ) const; // same as get_I2_vector
    inline void get_column( size_t i2, size_t i3, vector< I1, T >& data ) const; // same as get_I1_vector
    inline void get_tube( size_t i1, size_t i2, vector< I3, T >& data ) const; // same as get_I3_vector
	
    inline void set_I1_vector( size_t i2, size_t i3, const vector< I1, T >& data ); // I1_vector is a column vector with all values i1 at i2 and i3
    inline void set_I2_vector( size_t i1, size_t i3, const vector< I2, T >& data ); // I2_vector is a row vector with all values i2 at i1 and i3
    inline void set_I3_vector( size_t i1, size_t i2, const vector< I3, T >& data ); // I3_vector is a vector with all values i3 at a given i1 and i2

    inline void set_row( size_t i1, size_t i3, const vector< I2, T >& data ); // same as set_I2_vector
    inline void set_column( size_t i2, size_t i3, const vector< I1, T >& data ); // same as set_I1_vector
    inline void set_tube( size_t i1, size_t i2, const vector< I3, T >& data ); // same as set_I3_vector

    inline void get_frontal_slice_fwd( size_t i3, front_slice_type& data ) const;  
    inline void get_lateral_slice_bwd( size_t i2, bwd_lat_slice_type& data ) const; 
    inline void get_horizontal_slice_fwd( size_t i1, horiz_slice_type& data ) const;
	
    inline void get_frontal_slice_bwd( size_t i3, bwd_front_slice_type& data ) const;  
    inline void get_lateral_slice_fwd( size_t i2, lat_slice_type& data ) const; 
    inline void get_horizontal_slice_bwd( size_t i1, bwd_horiz_slice_type& data ) const;
    
    inline void set_frontal_slice_fwd( size_t i3, const front_slice_type& data ); 
    inline void set_lateral_slice_bwd( size_t i2, const bwd_lat_slice_type& data ); 
    inline void set_horizontal_slice_fwd( size_t i1, const horiz_slice_type& data );
	
    inline void set_frontal_slice_bwd( size_t i3, const bwd_front_slice_type& data ); 
    inline void set_lateral_slice_fwd( size_t i2, const lat_slice_type& data ); 
    inline void set_horizontal_slice_bwd( size_t i1, const bwd_horiz_slice_type& data );
	
    inline front_slice_type& get_frontal_slice_fwd( size_t index );
    inline const front_slice_type& get_frontal_slice_fwd( size_t index ) const;
	
	// sets all elements to fill_value
    void operator=( T fill_value ); //@SUS: todo
    void fill( T fill_value ); //special case of set method (all values are set to the same value!)
    void fill_random( );
    void fill_random_signed( );
    void fill_increasing_values( );
	
    template< size_t R >
    typename enable_if< R == I1 && R == I2 && R == I3 >::type* 
	 diag( const vector< R, T >& diag_values_ );
	
    void range_threshold(tensor3< I1, I2, I3, T >& other_, const T& start_value, const T& end_value) const;
    
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
	
	template< typename TT  >
		void quantize( tensor3< I1, I2, I3, TT >& quantized_, T& min_value, T& max_value ) const;
	template< typename TT  >
		void quantize_to( tensor3< I1, I2, I3, TT >& quantized_, const T& min_value, const T& max_value ) const;
	template< typename TT  >
	void quantize_log( tensor3< I1, I2, I3, TT >& quantized_, T& min_value, T& max_value ) const;
	template< typename TT  >
		void dequantize( tensor3< I1, I2, I3, TT >& dequantized_, const TT& min_value, const TT& max_value ) const;
	
    bool operator==( const tensor3& other ) const;
    bool operator!=( const tensor3& other ) const;
	
    // due to limited precision, two 'idential' tensor3 might seem different.
    // this function allows to specify a tolerance when comparing matrices.
    bool equals( const tensor3& other, T tolerance ) const;
    // this version takes a comparison functor to compare the components of
    // the two tensor3 data structures
    template< typename compare_t >
    bool equals( const tensor3& other, compare_t& cmp ) const;
	 
	
    //tensor times matrix multiplication along different modes
    template< size_t J1, size_t J2, size_t J3 > 
    void multiply_horizontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I3, J3, T >& other_slice_ ); //output: tensor3< J1, J2, I3, T >  
	
    template< size_t J1, size_t J2, size_t J3 > 
    void multiply_lateral_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I1, J1, T >& other_slice_ ); //output: tensor3< I1, J2, J3, T > 
	
    template< size_t J1, size_t J2, size_t J3 > 
    void multiply_frontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I2, J2, T >& other_slice_ ); //output: tensor3< J1, I2, J3, T >
    
    //backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
    template< size_t J1, size_t J2, size_t J3 > 
    void full_tensor3_matrix_multiplication( const tensor3< J1, J2, J3, T >& core, const matrix< I1, J1, T >& U1, const matrix< I2, J2, T >& U2, const matrix< I3, J3, T >& U3 );
    
    void horizontal_unfolding_bwd( bwd_horiz_unfolding_type& unfolding) const;
    void lateral_unfolding_bwd( bwd_lat_unfolding_type& unfolding) const;
    void frontal_unfolding_bwd( bwd_front_unfolding_type& unfolding) const;
    
    
    //error computation 
    double frobenius_norm( ) const;
    double rmse( const tensor3< I1, I2, I3, T >& other ) const; //root mean-squared error
    
    template< typename TT >
    void cast_from( const tensor3< I1, I2, I3, TT >& other );
	
    template< typename TT >
    void float_t_to_uint_t( const tensor3< I1, I2, I3, TT >& other );
	
    void export_to( std::vector< T >& data_ ) const ;
    void import_from( const std::vector< T >& data_ ) ;
	void write_to_raw( const std::string& dir_, const std::string& filename_ ) const;
	void read_from_raw( const std::string& dir_, const std::string& filename_ ) ;
	void write_datfile( const std::string& dir_, const std::string& filename_ ) const;
	    
    inline tensor3 operator+( T scalar ) const;
    inline tensor3 operator-( T scalar ) const;
    
    void operator+=( T scalar );
    void operator-=( T scalar );
    
    inline tensor3 operator+( const tensor3& other ) const;
    inline tensor3 operator-( const tensor3& other ) const;
    
    void operator+=( const tensor3& other );
    void operator-=( const tensor3& other );

    // 
    // tensor3-scalar operations / scaling 
    // 
    tensor3 operator*( T scalar );
    void operator*=( T scalar );
	
    //
    // matrix-vector operations
    //
    // transform column vector by matrix ( vec = matrix * vec )
    vector< I1, T > operator*( const vector< I2, T >& other ) const;
	
    // transform column vector by matrix ( vec = matrix * vec )
    // assume homogenous coords, e.g. vec3 = mat4x4 * vec3, with w = 1.0
    template< size_t O >
    vector< O, T > operator*( const vector< O, T >& vector_ ) const;
    
    inline tensor3< I1, I2, I3, T > operator-() const;
    tensor3< I1, I2, I3, T > negate() const;
	
    friend std::ostream& operator << ( std::ostream& os, const tensor3< I1, I2, I3, T >& t3 )
    {
        for(size_t i = 0; i < I3; ++i)
        {
                os << t3.array[ i ] << "***" << std::endl;			
        }	
        return os;
    }
	
	
    // static members
    static const tensor3< I1, I2, I3, T > ZERO;


protected:
	store_type			array;
	
}; // class tensor3

#define VMML_TEMPLATE_STRING    template< size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME tensor3< I1, I2, I3, T >


VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tensor3()
	: array()
{}

	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tensor3( const tensor3& source )
    : array( source.array )
{}


VMML_TEMPLATE_STRING
template< typename U >
VMML_TEMPLATE_CLASSNAME::tensor3( const tensor3< I1, I2, I3, U >& source_ )
{
	(*this) = source_;
}
	

VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 >
VMML_TEMPLATE_CLASSNAME::tensor3( const tensor3< J1, J2, J3, T >& source_ )
{
	const size_t minL =  J1 < I1 ? J1 : I1;
	const size_t minC =  J2 < I2 ? J2 : I2;
	const size_t minS =  J3 < I3 ? J3 : I3;
	
	zero();
	
	for ( size_t i = 0 ; i < minL ; i++ ) {
		for ( size_t j = 0 ; j < minC ; j++ ) {
			for ( size_t k = 0 ; k < minS ; k++ )
			{
				at( i,j, k ) = source_( i, j, k ); 
			}
		}
	}
}
	
	
VMML_TEMPLATE_STRING
inline T&
VMML_TEMPLATE_CLASSNAME::at( size_t i1, size_t i2, size_t i3 )
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( i3 >= I3 )
        VMMLIB_ERROR( "at( i1, i2, i3 ) - index out of bounds", VMMLIB_HERE );
    #endif
    return array[ i3 ].at( i1, i2 );
}



VMML_TEMPLATE_STRING
const inline T&
VMML_TEMPLATE_CLASSNAME::at( size_t i1, size_t i2, size_t i3 ) const
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( i3 >= I3 )
        VMMLIB_ERROR( "at( i1, i2, i3 ) - i3 index out of bounds", VMMLIB_HERE );
    #endif
    return array[ i3 ].at( i1, i2 );
}



VMML_TEMPLATE_STRING
inline T&
VMML_TEMPLATE_CLASSNAME::operator()( size_t i1, size_t i2, size_t i3 )
{
    return at( i1, i2, i3 );
}



VMML_TEMPLATE_STRING
const inline T&
VMML_TEMPLATE_CLASSNAME::operator()(  size_t i1, size_t i2, size_t i3 ) const
{
    return at( i1, i2, i3 );
}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_I2_vector( size_t i1, size_t i3, vector< I2, T >& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
	
    if ( i3 >= I3 )
        VMMLIB_ERROR( "get_I1_vector() - i3 index out of bounds.", VMMLIB_HERE );
	
#endif
    
	array[ i3 ].get_row( i1, data );	
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_I1_vector( size_t i2, size_t i3, vector< I1, T >& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
	
    if ( i3 >= I3 )
        VMMLIB_ERROR( "get_I2_vector() - i3 index out of bounds.", VMMLIB_HERE );
	
#endif
    
	array[ i3 ].get_column( i2, data );	
	
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_I3_vector( size_t i1, size_t i2, vector< I3, T >& data ) const
{
	for (size_t i3 = 0; i3 < I3; ++i3)
	{
		data[ i3 ] = array[ i3 ].at( i1, i2 );		
	}

}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_row( size_t i1, size_t i3, vector< I2, T >& data ) const
{
    get_I2_vector( i1, i3, data );
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_column( size_t i2, size_t i3, vector< I1, T >& data ) const
{
    get_I1_vector( i2, i3, data );
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_tube( size_t i1, size_t i2, vector< I3, T >& data ) const
{
	get_I3_vector( i1, i2, data );
}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_I2_vector( size_t i1, size_t i3, const vector< I2, T >& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
	
    if ( i3 >= I3 )
        VMMLIB_ERROR( "set_I1_vector() - i3 index out of bounds.", VMMLIB_HERE );
	
#endif
    
	array[ i3 ].set_row( i1, data );	//@SUS: bug fix
}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_I1_vector( size_t i2, size_t i3, const vector< I1, T >& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
	
    if ( i3 >= I3 )
        VMMLIB_ERROR( "set_I2_vector() - i3 index out of bounds.", VMMLIB_HERE );
	
#endif
    
	array[ i3 ].set_column( i2, data );	
	
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_I3_vector( size_t i1, size_t i2, const vector< I3, T >& data )
{
    for (size_t i3 = 0; i3 < I3; ++i3)
    {
            array[ i3 ].at( i1, i2 ) = data[ i3 ];
    }
	
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_row( size_t i1, size_t i3, const vector< I2, T >& data )
{
	set_I2_vector( i1, i3, data );
}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_column( size_t i2, size_t i3, const vector< I1, T >& data )
{
	set_I1_vector( i2, i3, data );
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_tube( size_t i1, size_t i2, const vector< I3, T >& data )
{
	set_I3_vector( i1, i2, data );
}

VMML_TEMPLATE_STRING
inline typename VMML_TEMPLATE_CLASSNAME::front_slice_type& 
VMML_TEMPLATE_CLASSNAME::
get_frontal_slice_fwd( size_t index )
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( index >= I3 )
        VMMLIB_ERROR( "get_frontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	return array[ index ];
}


VMML_TEMPLATE_STRING
inline const typename VMML_TEMPLATE_CLASSNAME::front_slice_type& 
VMML_TEMPLATE_CLASSNAME::
get_frontal_slice_fwd( size_t index ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( index >= I3 )
        VMMLIB_ERROR( "get_frontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	return array[ index ];
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_frontal_slice_fwd( size_t i3, front_slice_type& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( i3 >= I3 )
        VMMLIB_ERROR( "get_frontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	
	data = array[ i3 ];
}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_lateral_slice_bwd( size_t i2, bwd_lat_slice_type& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( i2 >= I2 )
        VMMLIB_ERROR( "get_lateral_slice_bwd() - index out of bounds.", VMMLIB_HERE );
#endif
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		data.set_column( i3, array[ i3 ].get_column( i2 ));
	}
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_horizontal_slice_fwd( size_t i1, horiz_slice_type& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( i1 >= I1 )
        VMMLIB_ERROR( "get_horizontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		data.set_column( i3, array[ i3 ].get_row( i1 )); //or for every i2 get/set column
	}
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_frontal_slice_bwd( size_t i3, bwd_front_slice_type& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( i3 >= I3 )
		VMMLIB_ERROR( "get_frontal_slice_bwd() - index out of bounds.", VMMLIB_HERE );
#endif
	
	front_slice_type* data_t = new front_slice_type();
	*data_t = array[ i3 ];
	data = transpose( *data_t );
	delete data_t;
}


VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_lateral_slice_fwd( size_t i2, lat_slice_type& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( i2 >= I2 )
		VMMLIB_ERROR( "get_lateral_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	bwd_lat_slice_type* data_t = new bwd_lat_slice_type();
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		data_t->set_column( i3, array[ i3 ].get_column( i2 ));
	}
	data = transpose( *data_t );
	delete data_t;
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
get_horizontal_slice_bwd( size_t i1, bwd_horiz_slice_type& data ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( i1 >= I1 )
		VMMLIB_ERROR( "get_horizontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	horiz_slice_type* data_t = new horiz_slice_type();
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		data_t->set_column( i3, array[ i3 ].get_row( i1 )); //or for every i2 get/set column
	}
	data = transpose( *data_t );
	delete data_t;
}
	
	
	
//setter
VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_frontal_slice_fwd( size_t i3, const front_slice_type& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( i3 >= I3 )
        VMMLIB_ERROR( "set_frontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	
	array[ i3 ] = data;
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_lateral_slice_bwd( size_t i2, const bwd_lat_slice_type& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( i2 >= I2 )
        VMMLIB_ERROR( "set_lateral_slice_bwd() - index out of bounds.", VMMLIB_HERE );
#endif
	
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		array[ i3 ].set_column(i2, data.get_column( i3 ));
	}
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_horizontal_slice_fwd( size_t i1, const horiz_slice_type& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
    if ( i1 >= I1 )
        VMMLIB_ERROR( "set_horizontal_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		array[ i3 ].set_row( i1, data.get_column( i3 ));
	}
	
}
	
VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_frontal_slice_bwd( size_t i3, const bwd_front_slice_type& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( i3 >= I3 )
		VMMLIB_ERROR( "set_frontal_slice_bwd() - index out of bounds.", VMMLIB_HERE );
#endif
	front_slice_type* data_t  = new front_slice_type();
	*data_t = transpose( data );
	array[ i3 ] = *data_t;
	delete data_t;
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_lateral_slice_fwd( size_t i2, const lat_slice_type& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( i2 >= I2 )
		VMMLIB_ERROR( "set_lateral_slice_fwd() - index out of bounds.", VMMLIB_HERE );
#endif
	bwd_lat_slice_type* data_t = new bwd_lat_slice_type();
	*data_t = transpose( data );
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		array[ i3 ].set_column(i2, data_t->get_column( i3 ));
	}
	
	delete data_t;
}

VMML_TEMPLATE_STRING
inline void 
VMML_TEMPLATE_CLASSNAME::
set_horizontal_slice_bwd( size_t i1, const bwd_horiz_slice_type& data )
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( i1 >= I1 )
		VMMLIB_ERROR( "set_horizontal_slice_bwd() - index out of bounds.", VMMLIB_HERE );
#endif
	horiz_slice_type* data_t = new horiz_slice_type();
	*data_t = transpose( data );
	
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		array[ i3 ].set_row( i1, data_t->get_column( i3 ));
	}
	delete data_t;
}


	
//fill
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::
fill( T fillValue )
{
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		array[ i3 ].fill( fillValue );
	}
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::
fill_random( )
{
	double fillValue = 0.0f;
	srand(time(NULL));
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		for( size_t i1 = 0; i1 < I1; ++i1 )
		{
			for( size_t i2 = 0; i2 < I2; ++i2 )
			{
				fillValue = rand();
				fillValue /= RAND_MAX;
				fillValue *= std::numeric_limits< T >::max();
				at( i1, i2, i3 ) = static_cast< T >( fillValue )  ;
			}
		}
	}
}

VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::
fill_random_signed( )
{
	double fillValue = 0.0f;
	srand(time(NULL));
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		for( size_t i1 = 0; i1 < I1; ++i1 )
		{
			for( size_t i2 = 0; i2 < I2; ++i2 )
			{
				fillValue = rand();
				fillValue /= RAND_MAX;
				fillValue *= std::numeric_limits< T >::max();
				T fillValue2 = static_cast< T >(fillValue) % std::numeric_limits< T >::max();
				fillValue2 -= std::numeric_limits< T >::max()/2;
				at( i1, i2, i3 ) = fillValue2  ;
			}
		}
	}
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::
fill_increasing_values( )
{
	double fillValue = 0.0f;
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		for( size_t i1 = 0; i1 < I1; ++i1 )
		{
			for( size_t i2 = 0; i2 < I2; ++i2 )
			{
				at( i1, i2, i3 ) = static_cast< T >( fillValue );
				fillValue++;
			}
		}
	}
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::range_threshold(tensor3<I1, I2, I3, T>& other_, const T& start_value, const T& end_value) const
{
	
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		for( size_t i1 = 0; i1 < I1; ++i1 )
		{
			for( size_t i2 = 0; i2 < I2; ++i2 )
			{
				T value = at( i1, i2, i3 );
				if (value >= start_value && value <= end_value)
					other_.at(i1, i2, i3 ) = static_cast< T >( value );
			}
		}
	}
}

	
VMML_TEMPLATE_STRING
template< size_t R>
typename enable_if< R == I1 && R == I2 && R == I3>::type* 
VMML_TEMPLATE_CLASSNAME::diag( const vector< R, T >& diag_values_ )
{
	zero();
	for( size_t r = 0; r < R; ++r )
	{
		at(r, r, r) = static_cast< T >( diag_values_.at(r) );
	}
}
	
	
	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::zero()
{
    fill( static_cast< T >( 0.0 ) );
}


VMML_TEMPLATE_STRING
bool
VMML_TEMPLATE_CLASSNAME::operator==( const tensor3< I1, I2, I3, T >& other ) const
{
    bool ok = true;
    for( size_t i3 = 0; ok && i3 < I3; ++i3 )
	{
        ok = array[ i3 ] == other.array[ i3 ];
	}
	
    return ok;
}


VMML_TEMPLATE_STRING
bool
VMML_TEMPLATE_CLASSNAME::operator!=( const tensor3< I1, I2, I3, T >& other ) const
{
    return ! operator==( other );
}


VMML_TEMPLATE_STRING
bool equals( const tensor3< I1, I2, I3, T >& t3_0, const tensor3< I1, I2, I3, T >& t3_1, T tolerance )
{
    return t3_0.equals( t3_1, tolerance );
}

VMML_TEMPLATE_STRING
bool
VMML_TEMPLATE_CLASSNAME::equals( const tensor3< I1, I2, I3, T >& other, T tolerance ) const
{
    bool ok = true;
	for (size_t i3 = 0; ok && i3 < I3; ++i3 ) 
	{
		ok = array[ i3 ].equals(other.array[ i3 ], tolerance);
	}
    return ok;	
}





VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::size() const
{
    return I1 * I2 * I3;
}

VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 >
tensor3<J1, J2, J3, T>
VMML_TEMPLATE_CLASSNAME::
get_sub_tensor3( size_t row_offset, size_t col_offset, size_t slice_offset,
				typename enable_if< J1 <= I1 && J2 <= I2 && J3 <= I3 >::type* ) const
{
	tensor3< J1, J2, J3, T > result;
	get_sub_tensor3( result, row_offset, col_offset, slice_offset );
	return result;
}


VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 >
typename enable_if< J1 <= I1 && J2 <= I2 && J3 <= I3 >::type*
VMML_TEMPLATE_CLASSNAME::
get_sub_tensor3( tensor3<J1, J2, J3, T >& result,  
				size_t row_offset, size_t col_offset, size_t slice_offset ) const
{
#ifdef VMMLIB_SAFE_ACCESSORS
	if ( J1 + row_offset > I1 || J2 + col_offset > I2 || J3 + slice_offset > I3 )
		VMMLIB_ERROR( "index out of bounds.", VMMLIB_HERE );
#endif
	
	for( size_t slice = 0; slice < J3; ++slice )
	{
		for( size_t row = 0; row < J1; ++row )
		{
			for( size_t col = 0; col < J2; ++col )
			{
				result.at( row, col, slice ) 
				= at( row_offset + row, col_offset + col, slice_offset + slice );
			}
		}
	}
	return 0;
}
	


VMML_TEMPLATE_STRING
typename VMML_TEMPLATE_CLASSNAME::iterator
VMML_TEMPLATE_CLASSNAME::begin()
{
    return iterator( *this, true );
}




VMML_TEMPLATE_STRING
typename VMML_TEMPLATE_CLASSNAME::iterator
VMML_TEMPLATE_CLASSNAME::end()
{
    return iterator( *this, false );
}



VMML_TEMPLATE_STRING
typename VMML_TEMPLATE_CLASSNAME::const_iterator
VMML_TEMPLATE_CLASSNAME::begin() const
{
    return const_iterator( *this, true );
}



VMML_TEMPLATE_STRING
typename VMML_TEMPLATE_CLASSNAME::const_iterator
VMML_TEMPLATE_CLASSNAME::end() const
{
    return const_iterator( *this, false );
}


VMML_TEMPLATE_STRING
template< typename input_iterator_t >
void 
VMML_TEMPLATE_CLASSNAME::set( input_iterator_t begin_, input_iterator_t end_, bool row_major_layout )
{
	input_iterator_t it( begin_ );
    if( row_major_layout )
    {
        for ( size_t i3 = 0; i3 < I3; ++i3 )
		{
			for( size_t i1 = 0; i1 < I1; ++i1 )
			{
				for( size_t i2 = 0; i2 < I2; ++i2, ++it )
				{
					if ( it == end_ )
						return;
					at( i1, i2, i3 ) = static_cast< T >( *it );
				}
			}
		}
    }
    else
    {
        std::copy( it, it + ( I1 * I2 * I3 ), begin() );
    }
}





VMML_TEMPLATE_STRING
inline VMML_TEMPLATE_CLASSNAME 
VMML_TEMPLATE_CLASSNAME::operator+( const tensor3< I1, I2, I3, T >& other ) const
{
	tensor3< I1, I2, I3, T > result( *this );
	result += other;
	return result;
}



VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::operator+=( const tensor3< I1, I2, I3, T >& other )
{
    iterator it = begin(), it_end = end(); 
    const_iterator other_it = other.begin();
    for( ; it != it_end; ++it, ++other_it )
    {
        *it += *other_it;
    }
}



VMML_TEMPLATE_STRING
inline VMML_TEMPLATE_CLASSNAME 
VMML_TEMPLATE_CLASSNAME::operator-( const tensor3< I1, I2, I3, T >& other ) const
{
	tensor3< I1, I2, I3, T > result( *this );
	result -= other;
	return result;
}



VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::operator-=( const tensor3< I1, I2, I3, T >& other )
{
    iterator it = begin(), it_end = end(); 
    const_iterator other_it = other.begin();
    for( ; it != it_end; ++it, ++other_it )
    {
        *it -= *other_it;
    }
}


//sum with scalar


VMML_TEMPLATE_STRING
inline VMML_TEMPLATE_CLASSNAME 
VMML_TEMPLATE_CLASSNAME::operator+( T scalar ) const
{
	tensor3< I1, I2, I3, T > result( *this );
	result += scalar;
	return result;
}



VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::operator+=( T scalar )
{
    iterator it = begin(), it_end = end(); 
    for( ; it != it_end; ++it )
    {
        *it += scalar;
    }
}



VMML_TEMPLATE_STRING
inline VMML_TEMPLATE_CLASSNAME 
VMML_TEMPLATE_CLASSNAME::operator-( T scalar ) const
{
	tensor3< I1, I2, I3, T > result( *this );
	result -= scalar;
	return result;
}



VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::operator-=( T scalar )
{
    iterator it = begin(), it_end = end(); 
    for( ; it != it_end; ++it )
    {
        *it -= scalar;
    }
}




//tensor matrix multiplications

VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 > 
void
VMML_TEMPLATE_CLASSNAME::multiply_horizontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I3, J3, T >& other_slice_ )
{
	matrix< J3, J2, T>* slice = new matrix< J3, J2, T>();
	matrix< I3, J2, T>* slice_new = new matrix< I3, J2, T>();
	for (size_t i1 = 0; i1 < J1; ++i1)
	{
		other.get_horizontal_slice_bwd( i1, *slice );
		slice_new->multiply( other_slice_, *slice );
		set_horizontal_slice_bwd( i1, *slice_new );		
	}
	
	delete slice;
	delete slice_new;
}

VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 > 
void
VMML_TEMPLATE_CLASSNAME::multiply_lateral_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I1, J1, T >& other_slice_ )
{
	matrix< J1, J3, T>* slice = new matrix< J1, J3, T>();
	matrix< I1, J3, T>* slice_new = new matrix< I1, J3, T>();
	for ( size_t i2 = 0; i2 < J2; ++i2 )
	{
		other.get_lateral_slice_bwd( i2, *slice );
		slice_new->multiply( other_slice_, *slice );
		set_lateral_slice_bwd( i2, *slice_new );		
	}
	delete slice;
	delete slice_new;
 }
 
 
VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 > 
void
VMML_TEMPLATE_CLASSNAME::multiply_frontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I2, J2, T >& other_slice_ )
{
	 matrix< J2, J1, T>* slice = new matrix< J2, J1, T>();
	 matrix< I2, J1, T>* slice_new = new matrix< I2, J1, T>();
	 for ( size_t i3 = 0; i3 < J3; ++i3 )
	 {
		 other.get_frontal_slice_bwd( i3, *slice ); 
		 slice_new->multiply( other_slice_, *slice );
		 set_frontal_slice_bwd( i3, *slice_new );		
	 }
	delete slice;
	delete slice_new;
}
 
 
 
VMML_TEMPLATE_STRING
template< size_t J1, size_t J2, size_t J3 > 
void
VMML_TEMPLATE_CLASSNAME::full_tensor3_matrix_multiplication(  const tensor3< J1, J2, J3, T >& core, 
														   const matrix< I1, J1, T >& U1, 
														   const matrix< I2, J2, T >& U2, 
														   const matrix< I3, J3, T >& U3 )
{
 tensor3< I1, J2, J3, T>* t3_result_1 = new  tensor3< I1, J2, J3, T>();
 tensor3< I1, I2, J3, T>* t3_result_2 = new tensor3< I1, I2, J3, T>();
 
 //backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
 t3_result_1->multiply_lateral_bwd( core, U1 );
 t3_result_2->multiply_frontal_bwd( *t3_result_1, U2 );
 multiply_horizontal_bwd( *t3_result_2, U3 );
 
	delete t3_result_1;
	delete t3_result_2;
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::horizontal_unfolding_bwd( bwd_horiz_unfolding_type& unfolding) const
{
	bwd_horiz_slice_type* horizontal_slice = new bwd_horiz_slice_type();
	for( size_t i = 0; i < I1; ++i )
	{
		get_horizontal_slice_bwd(i, *horizontal_slice );
		for( size_t col = 0; col < I2; ++col )
		{
			unfolding.set_column( i*I2+col,  horizontal_slice->get_column(col));
		} 
	}
	delete horizontal_slice;
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::lateral_unfolding_bwd( bwd_lat_unfolding_type& unfolding) const
{
	bwd_lat_slice_type* lateral_slice = new bwd_lat_slice_type();
	for( size_t i = 0; i < I2; ++i )
	{
		get_lateral_slice_bwd(i, *lateral_slice );
		for( size_t col = 0; col < I3; ++col )
		{
			unfolding.set_column( i*I3+col,  lateral_slice->get_column(col));
		} 
	}
	delete lateral_slice;
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::frontal_unfolding_bwd( bwd_front_unfolding_type& unfolding) const
{
	bwd_front_slice_type* frontal_slice = new bwd_front_slice_type();
	for( size_t i = 0; i < I3; ++i )
	{
		get_frontal_slice_bwd(i, *frontal_slice );
		for( size_t col = 0; col < I1; ++col )
		{
			unfolding.set_column( i*I1+col, frontal_slice->get_column(col));
		} 
	}
	delete frontal_slice;
}




VMML_TEMPLATE_STRING
tensor3< I1, I2, I3, T >
VMML_TEMPLATE_CLASSNAME::operator*( T scalar )
{
    tensor3< I1, I2, I3, T >* result = (*this);
    
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		result.array[ i3 ] = array[ i3 ] * scalar;
	}

    return *result;
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::operator*=( T scalar )
{
	for( size_t i3 = 0; i3 < I3; ++i3 )
	{
		array[ i3 ] *= scalar;
	}
	
}



VMML_TEMPLATE_STRING
inline tensor3< I1, I2, I3, T >
VMML_TEMPLATE_CLASSNAME::operator-() const
{
    return negate();
}


VMML_TEMPLATE_STRING
tensor3< I1, I2, I3, T >
VMML_TEMPLATE_CLASSNAME::negate() const
{
    tensor3< I1, I2, I3, T > result;
    result *= -1.0;
    return result;
}


VMML_TEMPLATE_STRING
double 
VMML_TEMPLATE_CLASSNAME::frobenius_norm( ) const
{
    double f_norm = 0.0;
    const_iterator it = begin(), it_end = end(); 
    for( ; it != it_end; ++it )
         f_norm += *it * *it;

    return sqrt(f_norm);
}

VMML_TEMPLATE_STRING
double 
VMML_TEMPLATE_CLASSNAME::rmse( const tensor3< I1, I2, I3, T >& other ) const
{
	double mse = 0.0;
	double diff = 0.0;
	const_iterator it = begin(), it_end = end(); 
	const_iterator other_it = other.begin(), other_it_end = other.end(); 
	for( ; it != it_end; ++it, ++other_it ){
		diff = abs( *it ) - abs( *other_it );
		mse += diff * diff;
	}
	
	return sqrt(mse/size());
}	
	
VMML_TEMPLATE_STRING
template< typename TT >
void
VMML_TEMPLATE_CLASSNAME::cast_from( const tensor3< I1, I2, I3, TT >& other )
{
    typedef tensor3< I1, I2, I3, TT > t3_tt_type ;
    typedef typename t3_tt_type::const_iterator tt_const_iterator;
    
    iterator it = begin(), it_end = end();
    tt_const_iterator other_it = other.begin();
    for( ; it != it_end; ++it, ++other_it )
    {
        *it = static_cast< T >( *other_it );
    }
}

VMML_TEMPLATE_STRING
template< typename TT >
void
VMML_TEMPLATE_CLASSNAME::float_t_to_uint_t( const tensor3< I1, I2, I3, TT >& other )
{
	typedef tensor3< I1, I2, I3, TT > t3_tt_type ;
	typedef typename t3_tt_type::const_iterator tt_const_iterator;
	
	if( sizeof(T) == 1 || sizeof(T) == 2) {
		iterator it = begin(), it_end = end();
		tt_const_iterator other_it = other.begin();
		for( ; it != it_end; ++it, ++other_it )
		{
			*it = T( std::min( std::max(int(0), int( *other_it + 0.5)), int(std::numeric_limits< T >::max()) ));
		}
	} 
	else {
		//std::cout << "Warning: use a different type as target (uint8 or uint16). No converstion done.\n" << std::endl;
		this->cast_from( other );
		return;
	}
}	

VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::export_to( std::vector< T >& data_ ) const
{	
    data_.clear();
    const_iterator  it = begin(),
    it_end = end();
    for( ; it != it_end; ++it )
    {
        data_.push_back( *it );
    }
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_from( const std::vector< T >& data_ )
{
    size_t i = 0; //iterator over data_
    size_t input_size = data_.size();
    
    iterator  it = begin(),
    it_end = end();
    for( ; it != it_end; ++it, ++i )
    {		
        if ( i >= input_size) 
                *it = static_cast<T> (0);
        else 
                *it = data_.at( i );
    }
}	

	
VMML_TEMPLATE_STRING
T
VMML_TEMPLATE_CLASSNAME::get_min() const
{
	T tensor3_min = static_cast<T>(std::numeric_limits<T>::max());
	
	const_iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( *it < tensor3_min ) {
			tensor3_min = *it;
		}
	}
	return tensor3_min;
}	
	
VMML_TEMPLATE_STRING
T
VMML_TEMPLATE_CLASSNAME::get_max() const
{
	T tensor3_max = static_cast<T>(0);
	
	const_iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( *it > tensor3_max ) {
			tensor3_max = *it;
		}
	}
	return tensor3_max;
}

VMML_TEMPLATE_STRING
T
VMML_TEMPLATE_CLASSNAME::get_abs_min() const
{
	T tensor3_min = static_cast<T>(std::numeric_limits<T>::max());
	
	const_iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( fabs(*it) < fabs(tensor3_min) ) {
			tensor3_min = fabs(*it);
		}
	}
	return tensor3_min;
}	

VMML_TEMPLATE_STRING
T
VMML_TEMPLATE_CLASSNAME::get_abs_max() const
{
	T tensor3_max = static_cast<T>(0);
	
	const_iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( fabs(*it) > fabs(tensor3_max) ) {
			tensor3_max = fabs(*it);
		}
	}
	return tensor3_max;
}
	
	
VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::nnz() const
{
	size_t counter = 0;
 	
	const_iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( *it != 0 ) {
			++counter;
		}
	}
	
	return counter;
}	
	
VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::nnz( const T& threshold_ ) const
{
	size_t counter = 0;
	
	const_iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( fabs(*it) > threshold_ ) {
			++counter;
		}
	}
	
	return counter;
}		
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::threshold( const T& threshold_value_ ) 
{
	iterator  it = begin(),
	it_end = end();
	for( ; it != it_end; ++it)
	{		
		if ( fabs(*it) <= threshold_value_ ) {
			*it = static_cast<T> (0);
		}
	}
}	
	
VMML_TEMPLATE_STRING
template< typename TT  >
void
VMML_TEMPLATE_CLASSNAME::quantize_to( tensor3< I1, I2, I3, TT >& quantized_, const T& min_value, const T& max_value ) const
{
	long max_tt_range = long(std::numeric_limits< TT >::max());
	long min_tt_range = long(std::numeric_limits< TT >::min());
	long tt_range = (max_tt_range - min_tt_range) + 1;
	
	T t_range = max_value - min_value;
	
	typedef tensor3< I1, I2, I3, TT > t3_tt_type ;
	typedef typename t3_tt_type::iterator tt_iterator;
	tt_iterator it_quant = quantized_.begin();
	const_iterator it = begin(), it_end = end();
	
	for( ; it != it_end; ++it, ++it_quant )
	{
		if (std::numeric_limits<TT>::is_signed ) {
			*it_quant = TT( std::min( std::max( min_tt_range, long(( *it * tt_range / t_range ) + 0.5)), max_tt_range ));
		} else {
			*it_quant = TT( std::min( std::max( min_tt_range, long(((*it - min_value) * tt_range / t_range) + 0.5)), max_tt_range ));
		}
	}
}
	
	
VMML_TEMPLATE_STRING
template< typename TT  >
void
VMML_TEMPLATE_CLASSNAME::quantize( tensor3< I1, I2, I3, TT >& quantized_, T& min_value, T& max_value ) const
{
	min_value = get_min();
	max_value = get_max();
	
	quantize_to( quantized_, min_value, max_value ); 
}	
	
	
VMML_TEMPLATE_STRING
template< typename TT  >
void
VMML_TEMPLATE_CLASSNAME::quantize_log( tensor3< I1, I2, I3, TT >& quantized_, T& min_value, T& max_value ) const
{
	long max_tt_range = long(std::numeric_limits< TT >::max());
	long min_tt_range = long(std::numeric_limits< TT >::min());
	long tt_range = max_tt_range - min_tt_range;
	
	min_value = get_min();
	max_value = get_max();
	T t_range = max_value - min_value;
	
	typedef tensor3< I1, I2, I3, TT > t3_tt_type ;
	typedef typename t3_tt_type::iterator tt_iterator;
	tt_iterator it_quant = quantized_.begin();
	const_iterator it = begin(), it_end = end();
	
	for( ; it != it_end; ++it, ++it_quant )
	{
		if (std::numeric_limits<TT>::is_signed ) {
			*it_quant = TT( std::min( std::max( min_tt_range, long(( *it * tt_range / t_range ) + 0.5)), max_tt_range ));
			std::cout << "quant value: " << float(log(*it) * (tt_range) / (t_range)) << std::endl;
		} else {
			*it_quant = TT( std::min( std::max( min_tt_range, long(((*it - min_value) * tt_range / t_range) + 0.5)), max_tt_range ));
			std::cout << "quant value: " << float((log(*it) - min_value) * (tt_range) / (t_range)) << std::endl;
		}
	}
}		

	
VMML_TEMPLATE_STRING
template< typename TT  >
void
VMML_TEMPLATE_CLASSNAME::dequantize( tensor3< I1, I2, I3, TT >& dequantized_, const TT& min_value, const TT& max_value ) const
{
	T max_t_range = get_max();
	T min_t_range = get_min();
	long t_range = long(max_t_range) - long(min_t_range) +1;
	
	TT tt_range = max_value - min_value;

	typedef tensor3< I1, I2, I3, TT > t3_tt_type ;
	typedef typename t3_tt_type::iterator tt_iterator;
	tt_iterator it_dequant = dequantized_.begin();
	const_iterator it = begin(), it_end = end();
	for( ; it != it_end; ++it, ++it_dequant )
	{
		if (std::numeric_limits<T>::is_signed ) {
			*it_dequant = std::min( std::max( min_value, TT((TT(*it) / t_range) * tt_range)), max_value );
		} else {
			*it_dequant = std::min( std::max( min_value, TT((((TT(*it) / t_range)) * tt_range ) + min_value)), max_value );
		}
	}
}	

#if 0
	
std::string format_path( const std::string& dir_, const std::string& filename_, const std::string& format_ )
{
	std::string path = dir_;
	int dir_length = dir_.size() -1;
	int last_separator = dir_.find_last_of( "/");
	std::string path = dir_;
	if (last_separator < dir_length ) {
		path.append( "/" );
	}
	path.append( filename_ );
	//check for format
	if( filename_.find( format_, filename_.size() -3) == (-1)) {
		path.append( ".");
		path.append( format_ );
	}
	return path;	
}

#endif	

	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::write_to_raw( const std::string& dir_, const std::string& filename_ ) const
{		
	int dir_length = dir_.size() -1;
	int last_separator = dir_.find_last_of( "/");
	std::string path = dir_;
	if (last_separator < dir_length ) {
		path.append( "/" );
	}
	path.append( filename_ );
	//check for format
	if( filename_.find( "raw", filename_.size() -3) == (-1)) {
		path.append( ".");
		path.append( "raw" );
	}
	std::string path_raw = path;
	
	std::ofstream outfile;	
	outfile.open( path_raw.c_str() );
	if( outfile.is_open() ) {
		const_iterator oit = begin(), oit_end = end();
		for( ; oit != oit_end; ++oit )
		{
			const T& byte = T(*oit);
			outfile << byte;
		}
	} else {
		std::cout << "no file open" << std::endl;
	}
	
	outfile.close();
}	


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::write_datfile( const std::string& dir_, const std::string& filename_ ) const
{	
	int dir_length = dir_.size() -1;
	int last_separator = dir_.find_last_of( "/");
	std::string path = dir_;
	if (last_separator < dir_length ) {
		path.append( "/" );
	}
	path.append( filename_ );
	//check for format
	if( filename_.find( "dat", filename_.size() -3) == (-1)) {
		path.append( ".");
		path.append( "dat" );
	}
	
	std::string path_dat = path;

	std::string filename = "";
	int pos = filename_.size() -4;
	if( ( filename_.find(".raw", pos ) == pos) || ( filename_.find(".dat", pos ) == pos) ) {
		filename = filename_.substr(0, filename_.size() -4);
	} else {
		filename = filename_;
	}
	
	const char* format = (sizeof(T) == sizeof(2)) ? "USHORT": "UCHAR";
	
	FILE* datfile = fopen(path_dat.c_str(), "w");
	fprintf(datfile, "ObjectFileName:\t%s.raw\n", filename.c_str());
	fprintf(datfile, "TaggedFileName:\t---\nResolution:\t%i %i %i\n", I1, I2, I3);
	fprintf(datfile, "SliceThickness:\t1.0 1.0 1.0\n");
	fprintf(datfile, "Format:\t%s\nNbrTags:\t0\n", format);
	fprintf(datfile, "ObjectType:\tTEXTURE_VOLUME_OBJECT\nObjectModel:\tI\nGridType:\tEQUIDISTANT\n");
	fprintf(datfile, "Modality:\tunknown\nTimeStep:\t0\n");
	fclose(datfile);
}	

	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::read_from_raw( const std::string& dir_, const std::string& filename_ ) 
{	
	int dir_length = dir_.size() -1;
	int last_separator = dir_.find_last_of( "/");
	std::string path = dir_;
	if (last_separator < dir_length ) {
		path.append( "/" );
	}
	path.append( filename_ );
	
	char * data;
	data = new char[ SIZE ];
	std::ifstream infile;
	infile.open(path.c_str(), std::ios::in); 
	
	if( infile.is_open())
	{
		iterator  it = begin(),
		it_end = end();
		infile.read( data, (sizeof(T) * SIZE));
		
		size_t counter = 0;
		for( ; it != it_end; ++it, ++counter)
		{
			*it =  data[counter];
		}		
		delete data;
		
	} else {
		std::cout << "no file open" << std::endl;
	}
	
	infile.close();
}	
	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME

} // namespace vmml

#endif

