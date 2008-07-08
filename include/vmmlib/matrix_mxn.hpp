#ifndef __VMML__MATRIX_MXN__HPP__
#define __VMML__MATRIX_MXN__HPP__

#include <vmmlib/array.hpp>

#include <iostream>
#include <iomanip>
#include <vector>

namespace vmml
{

// matrix of type float_t with m rows and n columns

template< size_t M, size_t N, typename float_t = double >
class matrix_mxn
{
public:

    // (this) matrix_mxn = left matrix_mxp * right matrix_pxn
    template< size_t P >
    void mul( 
        const matrix_mxn< M, P, float_t >& left,
        const matrix_mxn< P, N, float_t >& right 
        );
		
    // returned matrix_mxp = (this) matrix_mxn * other matrix_nxp;
    // use mul for performance reasons, it avoids a copy of the result matrix
    template< size_t P >
    matrix_mxn< M, P, float_t > operator*( matrix_mxn< N, P, float_t >& other ); 

	matrix_mxn< M, N, float_t > operator+( const matrix_mxn< M, N, float_t >& other ) const;
	void operator+=( const matrix_mxn< M, N, float_t >& other );

	template< size_t P, size_t Q >
	void direct_sum( const matrix_mxn< P, Q, float_t >& other, matrix_mxn< M + P, N + Q, float_t >& result );
	
	template< size_t P, size_t Q >
	matrix_mxn< M + P, N + Q, float_t > direct_sum( const matrix_mxn< P, Q, float_t >& other );

    // WARNING: data_array[] must be at least of size M * N - otherwise CRASH!
    // WARNING: assumes row_by_row layout - if this is not the case, 
    // use copyFrom1DimCArray( data_array, false )
    void operator=( const float_t* data_array );
    void operator=( const std::vector< float_t >& data );
    
    bool operator==( const matrix_mxn< M, N, float_t >& other );
    bool operator!=( const matrix_mxn< M, N, float_t >& other );
    // due to limited precision, two 'idential' matrices might seem different.
    // this function allows to specify a tolerance when comparing matrices.
    bool isEqualTo( const matrix_mxn< M, N, float_t >& other, float_t tolerance );
    
    matrix_mxn< M, N, float_t > operator*( float_t scalar );
    void operator*=( float_t scalar );

    // copies a transposed version of *this into transposedMatrix
    void transposeTo( matrix_mxn<N, M, float_t >& transposedMatrix ) const;
    // slower ( one additional copy )
    matrix_mxn< N, M, float_t > getTransposed() const;

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

    void getColumn( size_t columnNumber, array< float_t, N >& column ) const;
    void setColumn( size_t columnNumber, const array< float_t, N >& column );

    void getRow( size_t rowNumber, array< float_t, M >& row ) const;
    void setRow( size_t rowNumber,  const array< float_t, M >& row );
       
    float_t& at( size_t row, size_t column );
    
    array< float_t, N >& operator[]( size_t rowIndex );
    const array< float_t, N >& operator[]( size_t rowIndex ) const;

    size_t getM() const;
    size_t getNumberOfRows() const;
    
    size_t getN() const;
    size_t getNumberOfColumns() const;
    
    void fill( float_t fillValue );

    friend std::ostream& operator << ( std::ostream& os, 
        const matrix_mxn< M, N, float_t >& matrix )
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
                os << std::setw(10) << matrix._rows[ rowIndex ][ colIndex ] << " ";
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
    array< array< float_t, N >, M >    _rows;

}; // class matrix_mxn


typedef matrix_mxn< 3, 3, float > matrix_3x3_f;


template< size_t M, size_t N, typename float_t >
bool
matrix_mxn< M, N, float_t >::
operator==( const matrix_mxn< M, N, float_t >& other )
{
    bool ok = true;
    for( size_t rowIndex = 0; ok && rowIndex < M; rowIndex++)
    {
        for( size_t colIndex = 0; ok && colIndex < N; colIndex++) 
        {
            ok = _rows[ rowIndex ][ colIndex ] == other._rows[ rowIndex ][ colIndex ];
        }
    }
    return ok;
}



template< size_t M, size_t N, typename float_t >
bool
matrix_mxn< M, N, float_t >::
operator!=( const matrix_mxn< M, N, float_t >& other )
{
    return ! operator==( other );
}



template< size_t M, size_t N, typename float_t >
bool
matrix_mxn< M, N, float_t >::
isEqualTo( const matrix_mxn< M, N, float_t >& other, float_t tolerance )
{
    bool ok = true;
    for( size_t rowIndex = 0; ok && rowIndex < M; rowIndex++)
    {
        for( size_t colIndex = 0; ok && colIndex < N; colIndex++) 
        {
            ok = abs( _rows[ rowIndex ][ colIndex ] - other._rows[ rowIndex ][ colIndex ] ) < tolerance;
        }
    }
    return ok;
}


// WARNING: data_array[] must be at least of size M * N - otherwise CRASH!
// WARNING: assumes row_by_row layout - if this is not the case, 
// use copyFrom1DimCArray( data_array, false )
template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::operator=( const float_t* data_array )
{
    copyFrom1DimCArray( data_array, true );
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::operator=( const std::vector< float_t >& data )
{
    //assert( data.size() >= M * N );
    if ( data.size() >= M * N ) throw "FIXME.";
    copyFrom1DimCArray( &data[ 0 ], true );
}



template< size_t M, size_t N, typename float_t >
template< size_t P >
void
matrix_mxn< M, N, float_t >::mul( 
    const matrix_mxn< M, P, float_t >& left,
    const matrix_mxn< P, N, float_t >& right 
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
                tmp += left._rows[rowIndex][p] * right._rows[p][colIndex];
            }
            _rows[rowIndex][colIndex] = tmp;
        }
    }
}



template< size_t M, size_t N, typename float_t >
template< size_t P >
matrix_mxn< M, P, float_t >
matrix_mxn< M, N, float_t >::operator*( matrix_mxn< N, P, float_t >& other )
{
    matrix_mxn< M, P, float_t > result;
    result.mul( *this, other );
    return result;
}


template< size_t M, size_t N, typename float_t >
matrix_mxn< M, N, float_t >
matrix_mxn< M, N, float_t >::operator*( float_t scalar )
{
    matrix_mxn< M, N, float_t > result;
    
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
            result._rows[ rowIndex ][ colIndex ] = _rows[ rowIndex ][ colIndex ] * scalar;
        }
    }
    return result;
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::operator*=( float_t scalar )
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
            _rows[ rowIndex ][ colIndex ] *= scalar;
        }
    }
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
transposeTo( matrix_mxn< N, M, float_t >& tM ) const
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
            tM._rows[ colIndex ][ rowIndex ] = _rows[ rowIndex ][ colIndex ];
        }
    }
}



template< size_t M, size_t N, typename float_t >
matrix_mxn< N, M, float_t >
matrix_mxn< M, N, float_t >::getTransposed() const
{
    matrix_mxn< N, M, float_t > result;
    transposeTo( result );
    return result;
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
copyFrom1DimCArray( const float_t* c_array, bool row_by_row_layout )
{
    if ( row_by_row_layout )
    {
        for( size_t index = 0, rowIndex = 0; rowIndex < M; ++rowIndex )
        {
            for( size_t colIndex = 0; colIndex < N; ++colIndex, ++index )
            {
                _rows[ rowIndex ][ colIndex ] = c_array[ index ];
            }
        }
    }
    else
    {
        for( size_t index = 0, colIndex = 0; colIndex < N; ++colIndex )
        {
            for( size_t rowIndex = 0; rowIndex < M; ++rowIndex, ++index )
            {
                _rows[ rowIndex ][ colIndex ] = c_array[ index ];
            }
        }
    }
}


template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
getColumn( size_t columnNumber, array< float_t, N >& column ) const
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        column[ rowIndex ] = _rows[ rowIndex ][ columnNumber ];
    }
}


template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
setColumn( size_t columnNumber, const array< float_t, N >& column )
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        _rows[ rowIndex ][ columnNumber ] = column[ rowIndex ];
    }
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
getRow( size_t rowNumber, array< float_t, M >& row ) const
{
    row = _rows[ rowNumber ];
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
setRow( size_t rowNumber,  const array< float_t, M >& row )
{
    _rows[ rowNumber ] = row;
}



template< size_t M, size_t N, typename float_t >
float_t&
matrix_mxn< M, N, float_t >::
at( size_t row, size_t column )
{
    assert( row < M );
    assert( column < N );
    return _rows[ row ][ column ];
}



template< size_t M, size_t N, typename float_t >
array< float_t, N >&
matrix_mxn< M, N, float_t >::
operator[]( size_t rowIndex )
{
    if ( rowIndex >= M ) throw "FIXME.";
    assert( rowIndex < M );
    return _rows[ rowIndex ];
}



template< size_t M, size_t N, typename float_t >
const array< float_t, N >&
matrix_mxn< M, N, float_t >::
operator[]( size_t rowIndex ) const
{
    if ( rowIndex >= M ) throw "FIXME.";
    assert( rowIndex < M );
    return _rows[ rowIndex ];
}



template< size_t M, size_t N, typename float_t >
size_t
matrix_mxn< M, N, float_t >::
getM() const
{
    return M;
}



template< size_t M, size_t N, typename float_t >
size_t
matrix_mxn< M, N, float_t >::
getN() const
{ 
    return N; 
}



template< size_t M, size_t N, typename float_t >
size_t
matrix_mxn< M, N, float_t >::
getNumberOfRows() const
{
    return M;
}



template< size_t M, size_t N, typename float_t >
size_t
matrix_mxn< M, N, float_t >::
getNumberOfColumns() const
{ 
    return N; 
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
fill( float_t fillValue )
{
    for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
    {
        for( size_t colIndex = 0; colIndex < N; ++colIndex )
        {
			_rows[ rowIndex ][ colIndex ] = fillValue;
		}
    }
}



template< size_t M, size_t N, typename float_t >
matrix_mxn< M, N, float_t > 
matrix_mxn< M, N, float_t >::
operator+( const matrix_mxn< M, N, float_t >& other ) const
{
	matrix_mxn< M, N, float_t > result( *this );
	result += other;
	return result;
}



template< size_t M, size_t N, typename float_t >
void
matrix_mxn< M, N, float_t >::
operator+=( const matrix_mxn< M, N, float_t >& other )
{
	for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < N; ++colIndex )
		{
			_rows[ rowIndex ][ colIndex ] += other._rows[ rowIndex ][ colIndex ];
		}		
	}
}



template< size_t M, size_t N, typename float_t >
template< size_t P, size_t Q >
void
matrix_mxn< M, N, float_t >::
direct_sum( const matrix_mxn< P, Q, float_t >& other, matrix_mxn< M + P, N + Q, float_t >& result )
{
	result.fill( 0.0 );
	
	// copy this into result, upper-left part
	for( size_t rowIndex = 0; rowIndex < M; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < N; ++colIndex )
		{
			result._rows[ rowIndex ][ colIndex ] = _rows[ rowIndex ][ colIndex ];
		}
	}
	// copy other into result, lower-right part
	for( size_t rowIndex = 0; rowIndex < P; ++rowIndex )
	{
		for( size_t colIndex = 0; colIndex < Q; ++colIndex )
		{
			result._rows[ M + rowIndex ][ N + colIndex ] = other._rows[ rowIndex ][ colIndex ];
		}
	
	}

}


template< size_t M, size_t N, typename float_t >
template< size_t P, size_t Q >
matrix_mxn< M + P, N + Q, float_t > 
matrix_mxn< M, N, float_t >::
direct_sum( const matrix_mxn< P, Q, float_t >& other )
{
	matrix_mxn< M + P, N + Q, float_t > result;
	direct_sum( other, result );
	return result;
}



} // namespace vmml

#endif

