#ifndef __VMML__MATRIX__HPP__
#define __VMML__MATRIX__HPP__

#include <iostream>
#include <iomanip>
#include <vector>

#include <vmmlib/vmmlib_config.hpp>

#include <vmmlib/exception.hpp>
#include <vmmlib/vector.hpp>
#include <vmmlib/matrix_functors.hpp>

namespace vmml
{

// matrix of type float_t with m rows and n columns
template< size_t M, size_t N, typename float_t = double >
class matrix
{
public:
	static const bool is_square = M == N;
    
    matrix();
    matrix( const matrix_functor< matrix< M, N, float_t > >& functor );
    
    // accessor for matrix elements
    inline float_t& at( size_t rowIndex, size_t colIndex );
    inline const float_t& at( size_t rowIndex, size_t colIndex ) const;

	// legacy/compatibility accessor
	struct row_accessor
	{
		row_accessor( float_t* array_ ) : array( array_ ) {}
		float_t&
		operator[]( size_t colIndex )
		{ 
			#ifdef VMMLIB_SAFE_ACCESSORS
			if ( colIndex >= N ) 
                VMMLIB_ERROR( "column index out of bounds", VMMLIB_HERE );
			#endif
			return array[ colIndex * M ]; 
		}

		const float_t&
		operator[]( size_t colIndex ) const
		{ 
			#ifdef VMMLIB_SAFE_ACCESSORS
			if ( colIndex >= N ) 
                VMMLIB_ERROR( "column index out of bounds", VMMLIB_HERE );
			#endif
			return array[ colIndex * M ]; 
		}
		
		float_t* array;
		private: row_accessor() {} // disallow std ctor
	};
	// this is a hack to allow array-style access to matrix elements
	// usage: matrix< 2, 2, float > m; m[ 1 ][ 0 ] = 37.0f;
	inline row_accessor operator[]( size_t rowIndex )
	{ 
		#ifdef VMMLIB_SAFE_ACCESSORS
		if ( rowIndex > M ) 
            VMMLIB_ERROR( "row index out of bounds", VMMLIB_HERE );
		#endif
		return row_accessor( array + rowIndex );
	}

    // (this) matrix = left matrix_mxp * right matrix_pxn
    template< size_t P >
    void multiply( 
        const matrix< M, P, float_t >& left,
        const matrix< P, N, float_t >& right 
        );
		
    // WARNING: data_array[] must be at least of size M * N - otherwise CRASH!
    // WARNING: assumes row_by_row layout - if this is not the case, 
    // use copyFrom1DimCArray( data_array, false )
    void operator=( const float_t* data_array );
    void operator=( const std::vector< float_t >& data );

    void operator=( float_t fill_value );
    
    bool operator==( const matrix< M, N, float_t >& other );
    bool operator!=( const matrix< M, N, float_t >& other );
    // due to limited precision, two 'idential' matrices might seem different.
    // this function allows to specify a tolerance when comparing matrices.
    bool isEqualTo( const matrix< M, N, float_t >& other, float_t tolerance );
    
    matrix< M, N, float_t > operator*( float_t scalar );
    void operator*=( float_t scalar );

    // returned matrix_mxp = (this) matrix * other matrix_nxp;
    // use multiply(...) for performance reasons, it avoids a copy of the 
	// resulting matrix
    template< size_t P >
    matrix< M, P, float_t > operator*( matrix< N, P, float_t >& other ); 

	matrix< M, N, float_t > operator+( const matrix< M, N, float_t >& other ) const;
	void operator+=( const matrix< M, N, float_t >& other );

	template< size_t P, size_t Q >
	void direct_sum( const matrix< P, Q, float_t >& other, matrix< M + P, N + Q, float_t >& result );
	
	template< size_t P, size_t Q >
	matrix< M + P, N + Q, float_t > direct_sum( const matrix< P, Q, float_t >& other );

	template< size_t Mret, size_t Nret >
	matrix< Mret, Nret, float_t > getSubMatrix( size_t rowOffset, 
		size_t colOffset ) const;

	template< size_t Mret, size_t Nret >
	void getSubMatrix( matrix< Mret, Nret, float_t >& result, 
		size_t rowOffset, size_t colOffset ) const;

    // copies a transposed version of *this into transposedMatrix
    void transposeTo( matrix<N, M, float_t >& transposedMatrix ) const;
    // slower ( one additional copy )
    matrix< N, M, float_t > getTransposed() const;

    /*
    *
    * WARNING: data_array[] must be at least of size M * N - otherwise CRASH!
    *
    * if row_by_row_layout is true, the following layout in the c array is assumed:
    * e.g. for matrix_3x2
    * c_array[ 0 ] = row0 col0,
    * c_array[ 1 ] = row0 col1, 
    * c_array[ 2 ] = row1 col0, 
    * ...
    * 
    * otherwise, col_by_col is assumed: 
    * c_array[ 0 ] = row0 col0,
    * c_array[ 1 ] = row1 col0, 
    * c_array[ 2 ] = row2 col0, 
    * c_array[ 3 ] = row1 col1, 
    * ...    
    **/
    void copyFrom1DimCArray( const float_t* c_array, 
        bool row_by_row_layout = true );

    template< typename different_float_t >
    void copyFrom1DimCArray( const different_float_t* c_array, 
        bool row_by_row_layout = true );


    vector< M, float_t > getColumn( size_t columnNumber ) const;
    void getColumn( size_t columnNumber, vector< M, float_t >& column ) const;
    void setColumn( size_t columnNumber, const vector< M, float_t >& column );

    void getColumn( size_t columnNumber, matrix< M, 1, float_t >& column ) const;
    void setColumn( size_t columnNumber, const matrix< M, 1, float_t >& column );

    vector< N, float_t > getRow( size_t rowNumber ) const;
    void getRow( size_t rowNumber, vector< N, float_t >& row ) const;
    void setRow( size_t rowNumber,  const vector< N, float_t >& row );

    void getRow( size_t rowNumber, matrix< 1, N, float_t >& row ) const;
    void setRow( size_t rowNumber,  const matrix< 1, N, float_t >& row );
       
    size_t getM() const;
    size_t getNumberOfRows() const;
    
    size_t getN() const;
    size_t getNumberOfColumns() const;
    
    void fill( float_t fillValue );

    // square matrices only
    float_t computeDeterminant() const;
    void getAdjugate( matrix< M, M, float_t >& adjugate ) const;
	
    // the return value indicates if the matrix is invertible.
    // we need a tolerance term since the determinant is subject
	// to precision errors.
	bool computeInverse( matrix< M, M, float_t >& inverse, float_t tolerance = 1e-9 ) const;

    friend std::ostream& operator << ( std::ostream& os, 
        const matrix< M, N, float_t >& matrix )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        
        os << "\n";
        for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
        {
            os << "|";
            for( size_t colIndex = 0; colIndex < N; ++colIndex )
            {
                os << std::setw(10) << matrix.at( rowIndex, colIndex ) << " ";
            }
            os << "|\n";
        }
        os << std::endl;
        os.precision( prec );
        os.setf( flags );
        return os;
    };  

    // protected:
    // matrix.array[ row ][ column ];
    //array< array< float_t, N >, M >    _rows;
    // column_by_column
    float_t array[ M * N ];
    

    // static members
    static const matrix< M, N, float_t > IDENTITY;
    static const matrix< M, N, float_t > ZERO;
    static const set_to_identity< matrix< M, N, float_t > > IDENTITY_FUNCTOR;
    static const set_to_zero< matrix< M, N, float_t > > ZERO_FUNCTOR;

}; // class matrix

/*
template< size_t M, size_t N, typename float_t >
const matrix< M, N, float_t > 
matrix< M, N, float_t >::ZERO( set_to_zero< matrix< M, N, float_t > > );
*/

template< size_t M, size_t N, typename float_t >
const set_to_identity< matrix< M, N, float_t > >
matrix< M, N, float_t >::IDENTITY_FUNCTOR;

template< size_t M, size_t N, typename float_t >
const set_to_zero< matrix< M, N, float_t > >
matrix< M, N, float_t >::ZERO_FUNCTOR;

template< size_t M, size_t N, typename float_t >
const matrix< M, N, float_t > 
matrix< M, N, float_t >::IDENTITY( IDENTITY_FUNCTOR );

template< size_t M, size_t N, typename float_t >
const matrix< M, N, float_t > 
matrix< M, N, float_t >::ZERO( ZERO_FUNCTOR );


template< size_t M, size_t N, typename float_t >
matrix< M, N, float_t >::matrix()
{
    // no initialization for performance reasons.
}



template< size_t M, size_t N, typename float_t >
matrix< M, N, float_t >::matrix( const matrix_functor< matrix< M, N, float_t > >& functor )
{
    functor( *this );
}



template< size_t M, size_t N, typename float_t >
inline float_t&
matrix< M, N, float_t >::at( size_t rowIndex, size_t colIndex )
{
	#ifdef VMMLIB_SAFE_ACCESSORS
	if ( rowIndex >= M || colIndex >= N )
        VMMLIB_ERROR( "at( row, col ) - index out of bounds", VMMLIB_HERE );
	#endif
    return array[ colIndex * M + rowIndex ];
}



template< size_t M, size_t N, typename float_t >
const inline float_t&
matrix< M, N, float_t >::at( size_t rowIndex, size_t colIndex ) const
{
	#ifdef VMMLIB_SAFE_ACCESSORS
	if ( rowIndex >= M || colIndex >= N )
        VMMLIB_ERROR( "at( row, col ) - index out of bounds", VMMLIB_HERE );
	#endif
    return array[ colIndex * M + rowIndex ];
}


#if 0
template< size_t M, size_t N, typename float_t >
matrix< M, N, float_t >::matrix()
{
}



template< size_t M, size_t N, typename float_t >
matrix< M, N, float_t >::matrix( const matrix< M, N >& original )
{

}
#endif


template< size_t M, size_t N, typename float_t >
bool
matrix< M, N, float_t >::
operator==( const matrix< M, N, float_t >& other )
{
    bool ok = true;
    for( size_t i = 0; ok && i < M * N; ++i )
    {
        ok = array[ i ] == other.array[ i ];
    }
    return ok;
}



template< size_t M, size_t N, typename float_t >
bool
matrix< M, N, float_t >::
operator!=( const matrix< M, N, float_t >& other )
{
    return ! operator==( other );
}



template< size_t M, size_t N, typename float_t >
bool
matrix< M, N, float_t >::
isEqualTo( const matrix< M, N, float_t >& other, float_t tolerance )
{
    bool ok = true;
    for( size_t rowIndex = 0; ok && rowIndex < M; rowIndex++)
    {
        for( size_t colIndex = 0; ok && colIndex < N; colIndex++) 
        {
            ok = abs( at( rowIndex, colIndex ) - other.at( rowIndex, colIndex ) ) < tolerance;
        }
    }
    return ok;
}


// WARNING: data_array[] must be at least of size M * N - otherwise CRASH!
// WARNING: assumes row_by_row layout - if this is not the case, 
// use copyFrom1DimCArray( data_array, false )
template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::operator=( const float_t* data_array )
{
    copyFrom1DimCArray( data_array, true );
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::operator=( const std::vector< float_t >& data )
{
    //assert( data.size() >= M * N );
    if ( data.size() >= M * N )
        VMMLIB_ERROR( "input data vector too small", VMMLIB_HERE );
    copyFrom1DimCArray( &data[ 0 ], true );
}



template< size_t M, size_t N, typename float_t >
template< size_t P >
void
matrix< M, N, float_t >::multiply( 
    const matrix< M, P, float_t >& left,
    const matrix< P, N, float_t >& right 
    )
{
    float_t tmp;
    for( size_t rowIndex = 0; rowIndex < M; rowIndex++)
    {
        for( size_t colIndex = 0; colIndex < N; colIndex++) 
        {
            tmp = static_cast< float_t >( 0.0 );
            for( size_t p = 0; p < P; p++)
            {
                tmp += left.at( rowIndex, p ) * right.at( p, colIndex );
            }
            at( rowIndex, colIndex ) = tmp;
        }
    }
}



template< size_t M, size_t N, typename float_t >
template< size_t P >
matrix< M, P, float_t >
matrix< M, N, float_t >::operator*( matrix< N, P, float_t >& other )
{
    matrix< M, P, float_t > result;
    result.multiply( *this, other );
    return result;
}


template< size_t M, size_t N, typename float_t >
matrix< M, N, float_t >
matrix< M, N, float_t >::operator*( float_t scalar )
{
    matrix< M, N, float_t > result;
    
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
            result.at( rowIndex, colIndex ) = at( rowIndex, colIndex ) * scalar;
        }
    }
    return result;
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::operator*=( float_t scalar )
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
            at( rowIndex, colIndex ) *= scalar;
        }
    }
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
transposeTo( matrix< N, M, float_t >& tM ) const
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
            tM.at( colIndex, rowIndex ) = at( rowIndex, colIndex );
        }
    }
}



template< size_t M, size_t N, typename float_t >
matrix< N, M, float_t >
matrix< M, N, float_t >::getTransposed() const
{
    matrix< N, M, float_t > result;
    transposeTo( result );
    return result;
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
copyFrom1DimCArray( const float_t* c_array, bool row_by_row_layout )
{
    if ( row_by_row_layout )
    {
        for( size_t index = 0, rowIndex = 0; rowIndex < M; ++rowIndex )
        {
            for( size_t colIndex = 0; colIndex < N; ++colIndex, ++index )
            {
                at( rowIndex, colIndex ) = c_array[ index ];
            }
        }
    }
    else
    {
        memcpy( array, c_array, M * N * sizeof( float_t ) );
    }
}


template< size_t M, size_t N, typename float_t >
template< typename different_float_t >
void
matrix< M, N, float_t >::
copyFrom1DimCArray( const different_float_t* c_array, bool row_by_row_layout )
{
    if ( row_by_row_layout )
    {
        for( size_t index = 0, rowIndex = 0; rowIndex < M; ++rowIndex )
        {
            for( size_t colIndex = 0; colIndex < N; ++colIndex, ++index )
            {
                at( rowIndex, colIndex ) = static_cast< float_t >( c_array[ index ] );
            }
        }
    }
    else
    {
        for( size_t index = 0, colIndex = 0; colIndex < N; ++colIndex )
        {
            for( size_t rowIndex = 0; rowIndex < M; ++rowIndex, ++index  )
            {
                at( rowIndex, colIndex ) = static_cast< float_t >( c_array[ index ] );
            }
        }
    }
}



template< size_t M, size_t N, typename float_t >
vector< M, float_t > 
matrix< M, N, float_t >::
getColumn( size_t columnNumber ) const
{
	vector< M, float_t > column;
	getColumn( columnNumber, column );
	return column;
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
getColumn( size_t columnNumber, vector< M, float_t >& column ) const
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( columnNumber >= N )
        VMMLIB_ERROR( "getColumn() - index out of bounds.", VMMLIB_HERE );
    #endif

    memcpy( &column.array[0], &array[ M * columnNumber ], M * sizeof( float_t ) );
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
setColumn( size_t columnNumber, const vector< M, float_t >& column )
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( columnNumber >= N )
        VMMLIB_ERROR( "setColumn() - index out of bounds.", VMMLIB_HERE );
    #endif

    memcpy( &array[ M * columnNumber ], &column.array[0], M * sizeof( float_t ) );
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
getColumn( size_t columnNumber, matrix< M, 1, float_t >& column ) const
{
    #if 0
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        //column[ rowIndex ] = _rows[ rowIndex ][ columnNumber ];
        column.at( rowIndex, 1 ) = at( rowIndex, columnNumber );
    }
    #endif

    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( columnNumber >= N )
        VMMLIB_ERROR( "getColumn() - index out of bounds.", VMMLIB_HERE );
    #endif

    memcpy( &column.array[0], &array[ M * columnNumber ], M * sizeof( float_t ) );
}


template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
setColumn( size_t columnNumber, const matrix< M, 1, float_t >& column )
{
    #if 0
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        at( rowIndex, columnNumber ) = column.at( rowIndex, 0 );
    }
    #else

    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( columnNumber >= N )
        VMMLIB_ERROR( "setColumn() - index out of bounds.", VMMLIB_HERE );
    #endif

    memcpy( &array[ M * columnNumber ], &column.array[0], M * sizeof( float_t ) );

    #endif
}


template< size_t M, size_t N, typename float_t >
vector< N, float_t > 
matrix< M, N, float_t >::
getRow( size_t rowIndex ) const
{
	vector< N, float_t > row;
	getRow( rowIndex, row );
	return row;
}


template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
getRow( size_t rowIndex, vector< N, float_t >& row ) const
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( rowIndex >= M )
        VMMLIB_ERROR( "getRow() - index out of bounds.", VMMLIB_HERE );
    #endif

    for( size_t colIndex = 0; colIndex < N; ++colIndex )
    {
        row.at( colIndex ) = at( rowIndex, colIndex );
    }
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
setRow( size_t rowIndex,  const vector< N, float_t >& row )
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( rowIndex >= M )
        VMMLIB_ERROR( "setRow() - index out of bounds.", VMMLIB_HERE );
    #endif

    for( size_t colIndex = 0; colIndex < N; ++colIndex )
    {
        at( rowIndex, colIndex ) = row.at( colIndex );
    }
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
getRow( size_t rowIndex, matrix< 1, N, float_t >& row ) const
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( rowIndex >= M )
        VMMLIB_ERROR( "getRow() - index out of bounds.", VMMLIB_HERE );
    #endif

    for( size_t colIndex = 0; colIndex < N; ++colIndex )
    {
        row.at( 0, colIndex ) = at( rowIndex, colIndex );
    }
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
setRow( size_t rowIndex,  const matrix< 1, N, float_t >& row )
{
    #ifdef VMMLIB_SAFE_ACCESSORS
    if ( rowIndex >= M )
        VMMLIB_ERROR( "getRow() - index out of bounds.", VMMLIB_HERE );
    #endif

    for( size_t colIndex = 0; colIndex < N; ++colIndex )
    {
        at( rowIndex, colIndex ) = row.at( 0, colIndex );
    }
}



/*
template< size_t M, size_t N, typename float_t >
matrix< 1, N, float_t >&
matrix< M, N, float_t >::
operator[]( size_t rowIndex )
{
    if ( rowIndex >= M ) throw "FIXME.";
    assert( rowIndex < M );
    return _rows[ rowIndex ];
}



template< size_t M, size_t N, typename float_t >
const matrix< 1, N, float_t >&
matrix< M, N, float_t >::
operator[]( size_t rowIndex ) const
{
    if ( rowIndex >= M ) throw "FIXME.";
    assert( rowIndex < M );
    return _rows[ rowIndex ];
}


*/
template< size_t M, size_t N, typename float_t >
size_t
matrix< M, N, float_t >::
getM() const
{
    return M;
}



template< size_t M, size_t N, typename float_t >
size_t
matrix< M, N, float_t >::
getN() const
{ 
    return N; 
}



template< size_t M, size_t N, typename float_t >
size_t
matrix< M, N, float_t >::
getNumberOfRows() const
{
    return M;
}



template< size_t M, size_t N, typename float_t >
size_t
matrix< M, N, float_t >::
getNumberOfColumns() const
{ 
    return N; 
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
fill( float_t fillValue )
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
			at( rowIndex, colIndex ) = fillValue;
		}
    }
}


template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
operator=( float_t fillValue )
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
			at( rowIndex, colIndex ) = fillValue;
		}
    }
}



template< size_t M, size_t N, typename float_t >
matrix< M, N, float_t > 
matrix< M, N, float_t >::
operator+( const matrix< M, N, float_t >& other ) const
{
	matrix< M, N, float_t > result( *this );
	result += other;
	return result;
}



template< size_t M, size_t N, typename float_t >
void
matrix< M, N, float_t >::
operator+=( const matrix< M, N, float_t >& other )
{
	for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < N; ++colIndex )
		{
			at( rowIndex, colIndex ) += other.at( rowIndex, colIndex );
		}		
	}
}



template< size_t M, size_t N, typename float_t >
template< size_t P, size_t Q >
void
matrix< M, N, float_t >::
direct_sum( const matrix< P, Q, float_t >& other, matrix< M + P, N + Q, float_t >& result )
{
	result.fill( 0.0 );
	
	// copy this into result, upper-left part
	for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < N; ++colIndex )
		{
			result.at( rowIndex, colIndex ) = at( rowIndex, colIndex );
		}
	}
	// copy other into result, lower-right part
	for( size_t rowIndex = 0; rowIndex < P; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < Q; ++colIndex )
		{
			result.at( M + rowIndex, N + colIndex ) = other.at( rowIndex, colIndex );
		}
	
	}

}


template< size_t M, size_t N, typename float_t >
template< size_t P, size_t Q >
matrix< M + P, N + Q, float_t > 
matrix< M, N, float_t >::
direct_sum( const matrix< P, Q, float_t >& other )
{
	matrix< M + P, N + Q, float_t > result;
	direct_sum( other, result );
	return result;
}


template< size_t M, size_t N, typename float_t >
template< size_t Mret, size_t Nret >
matrix< Mret, Nret, float_t >
matrix< M, N, float_t >::getSubMatrix( size_t rowOffset, 
	size_t colOffset ) const
{
	matrix< Mret, Nret, float_t > result;
	getSubMatrix( result, rowOffset, colOffset );
	return result;
}



template< size_t M, size_t N, typename float_t >
template< size_t Mret, size_t Nret >
void
matrix< M, N, float_t >::getSubMatrix( matrix< Mret, Nret, float_t >& result, 
	size_t rowOffset, size_t colOffset ) const
{
	for( size_t rowIndex = 0; rowIndex < Mret; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < Nret; ++colIndex )
		{
			result.at( rowIndex, colIndex ) 
				= at( rowOffset + rowIndex, colOffset + colIndex );
		}
	}
}



template< size_t M, size_t N, typename float_t >
inline float_t
matrix< M, N, float_t >::computeDeterminant() const
{
    if ( ! is_square )
        VMMLIB_ERROR( "determinant is not defined for non-square matrices.", VMMLIB_HERE );
    else
        VMMLIB_ERROR( "not implemented yet.", VMMLIB_HERE );

}


// specializations ... 
template<>
inline float
matrix< 2, 2, float >::computeDeterminant() const
{
        const float& a = at( 0, 0 );
        const float& b = at( 0, 1 );
        const float& c = at( 1, 0 );
        const float& d = at( 1, 1 );
        return a * d - b * c;
}


template<>
inline double
matrix< 2, 2, double >::computeDeterminant() const
{
        const double& a = at( 0, 0 );
        const double& b = at( 0, 1 );
        const double& c = at( 1, 0 );
        const double& d = at( 1, 1 );
        return a * d - b * c;
}



template< size_t M, size_t N, typename float_t >
inline void
matrix< M, N, float_t >::getAdjugate( matrix< M, M, float_t >& adjugate ) const
{
    if ( ! is_square )
        VMMLIB_ERROR( "adjugate matrix is not defined for non-square matrices.", VMMLIB_HERE );
    else
        VMMLIB_ERROR( "not implemented yet.", VMMLIB_HERE );
}



template<>
inline void
matrix< 2, 2, float >::getAdjugate( matrix< 2, 2, float >& adjugate ) const
{
    const float& a = at( 0, 0 );
    const float& b = at( 0, 1 );
    const float& c = at( 1, 0 );
    const float& d = at( 1, 1 );

    adjugate.at( 0, 0 ) = d;
    adjugate.at( 0, 1 ) = -b;
    adjugate.at( 1, 0 ) = -c;
    adjugate.at( 1, 1 ) = a;
}



template<>
inline void
matrix< 2, 2, double >::getAdjugate( matrix< 2, 2, double >& adjugate ) const
{
    const double& a = at( 0, 0 );
    const double& b = at( 0, 1 );
    const double& c = at( 1, 0 );
    const double& d = at( 1, 1 );

    adjugate.at( 0, 0 ) = d;
    adjugate.at( 0, 1 ) = -b;
    adjugate.at( 1, 0 ) = -c;
    adjugate.at( 1, 1 ) = a;
}



template< size_t M, size_t N, typename float_t >
inline bool
matrix< M, N, float_t >::computeInverse( matrix< M, M, float_t >& Minverse, float_t tolerance ) const
{
    if ( ! is_square )
        VMMLIB_ERROR( "inverse of a matrix is not defined for non-square matrices.", VMMLIB_HERE );
    else
        VMMLIB_ERROR( "not implemented yet.", VMMLIB_HERE );
        
    return false;
}



template<>
inline bool
matrix< 2, 2, float >::computeInverse( matrix< 2, 2, float >& Minverse, float tolerance ) const
{
    float det = computeDeterminant();
    if ( det > tolerance )
    {
        return false;
    }
    float reciprocal_of_determinant = 1.0 / det;
    
    // set Minverse to the adjugate of M
    getAdjugate( Minverse );
    
    Minverse *= reciprocal_of_determinant;
    
    return true;
}


template<>
inline bool
matrix< 2, 2, double >::computeInverse( matrix< 2, 2, double >& Minverse, double tolerance ) const
{
    float det = computeDeterminant();
    if ( det > tolerance )
    {
        return false;
    }
    float reciprocal_of_determinant = 1.0 / det;
    
    // set Minverse to the adjugate of M
    getAdjugate( Minverse );
    
    Minverse *= reciprocal_of_determinant;
    
    return true;
}



template<>
inline bool
matrix< 3, 3, double >::computeInverse( matrix< 3, 3, double >& result, double tolerance ) const
{
    // Invert a 3x3 using cofactors.  This is about 8 times faster than
    // the Numerical Recipes code which uses Gaussian elimination.

    result.at( 0, 0 ) = at( 1, 1 ) * at( 2, 2 ) - at( 1, 2 ) * at( 2, 1 );
    result.at( 0, 1 ) = at( 0, 2 ) * at( 2, 1 ) - at( 0, 1 ) * at( 2, 2 );
    result.at( 0, 2 ) = at( 0, 1 ) * at( 1, 2 ) - at( 0, 2 ) * at( 1, 1 );
    result.at( 1, 0 ) = at( 1, 2 ) * at( 2, 0 ) - at( 1, 0 ) * at( 2, 2 );
    result.at( 1, 1 ) = at( 0, 0 ) * at( 2, 2 ) - at( 0, 2 ) * at( 2, 0 );
    result.at( 1, 2 ) = at( 0, 2 ) * at( 1, 0 ) - at( 0, 0 ) * at( 1, 1 );
    result.at( 2, 0 ) = at( 1, 0 ) * at( 2, 1 ) - at( 1, 1 ) * at( 2, 0 );
    result.at( 2, 1 ) = at( 0, 1 ) * at( 2, 0 ) - at( 0, 0 ) * at( 2, 1 );
    result.at( 2, 2 ) = at( 0, 0 ) * at( 1, 1 ) - at( 0, 1 ) * at( 1, 0 );
    
    const double determinant = at( 0, 0 ) * result.at( 0, 0 ) 
        + at( 0, 1 ) * result.at( 1, 0 )
        + at( 0, 2 ) * result.at( 2, 0 );
    
    if ( abs( determinant ) <= tolerance )
        return false; // matrix is not invertible

    const double detinv = 1.0 / determinant;
    
    result.at( 0, 0 ) *= detinv;
    result.at( 0, 1 ) *= detinv;
    result.at( 0, 2 ) *= detinv;
    result.at( 1, 0 ) *= detinv;
    result.at( 1, 1 ) *= detinv;
    result.at( 1, 2 ) *= detinv;
    result.at( 2, 0 ) *= detinv;
    result.at( 2, 1 ) *= detinv;
    result.at( 2, 2 ) *= detinv;

    return true;
}


template<>
inline bool
matrix< 4, 4, double >::computeInverse( matrix< 4, 4, double >& result, double tolerance ) const
{
    // tuned version from Claude Knaus
    /* first set of 2x2 determinants: 12 multiplications, 6 additions */
    const double t1[6] = { array[ 2] * array[ 7] - array[ 6] * array[ 3],
                      array[ 2] * array[11] - array[10] * array[ 3],
                      array[ 2] * array[15] - array[14] * array[ 3],
                      array[ 6] * array[11] - array[10] * array[ 7],
                      array[ 6] * array[15] - array[14] * array[ 7],
                      array[10] * array[15] - array[14] * array[11] };

    /* first half of comatrix: 24 multiplications, 16 additions */
    result.array[0] = array[ 5] * t1[5] - array[ 9] * t1[4] + array[13] * t1[3];
    result.array[1] = array[ 9] * t1[2] - array[13] * t1[1] - array[ 1] * t1[5];
    result.array[2] = array[13] * t1[0] - array[ 5] * t1[2] + array[ 1] * t1[4];
    result.array[3] = array[ 5] * t1[1] - array[ 1] * t1[3] - array[ 9] * t1[0];
    result.array[4] = array[ 8] * t1[4] - array[ 4] * t1[5] - array[12] * t1[3];
    result.array[5] = array[ 0] * t1[5] - array[ 8] * t1[2] + array[12] * t1[1];
    result.array[6] = array[ 4] * t1[2] - array[12] * t1[0] - array[ 0] * t1[4];
    result.array[7] = array[ 0] * t1[3] - array[ 4] * t1[1] + array[ 8] * t1[0];

   /* second set of 2x2 determinants: 12 multiplications, 6 additions */
    const double t2[6] = { array[ 0] * array[ 5] - array[ 4] * array[ 1],
                      array[ 0] * array[ 9] - array[ 8] * array[ 1],
                      array[ 0] * array[13] - array[12] * array[ 1],
                      array[ 4] * array[ 9] - array[ 8] * array[ 5],
                      array[ 4] * array[13] - array[12] * array[ 5],
                      array[ 8] * array[13] - array[12] * array[ 9] };

    /* second half of comatrix: 24 multiplications, 16 additions */
    result.array[8]  = array[ 7] * t2[5] - array[11] * t2[4] + array[15] * t2[3];
    result.array[9]  = array[11] * t2[2] - array[15] * t2[1] - array[ 3] * t2[5];
    result.array[10] = array[15] * t2[0] - array[ 7] * t2[2] + array[ 3] * t2[4];
    result.array[11] = array[ 7] * t2[1] - array[ 3] * t2[3] - array[11] * t2[0];
    result.array[12] = array[10] * t2[4] - array[ 6] * t2[5] - array[14] * t2[3];
    result.array[13] = array[ 2] * t2[5] - array[10] * t2[2] + array[14] * t2[1];
    result.array[14] = array[ 6] * t2[2] - array[14] * t2[0] - array[ 2] * t2[4];
    result.array[15] = array[ 2] * t2[3] - array[ 6] * t2[1] + array[10] * t2[0];

   /* determinant: 4 multiplications, 3 additions */
   const double determinant = array[0] * result.array[0] + array[4] * result.array[1] +
                         array[8] * result.array[2] + array[12] * result.array[3];

   if( abs( determinant ) <= tolerance )
        return false; // matrix is not invertible

   /* division: 16 multiplications, 1 division */
   const double detinv = 1.0 / determinant;
   for( unsigned i = 0; i != 16; ++i )
       result.array[i] *= detinv;
       
    return true;
}




} // namespace vmml

#endif

