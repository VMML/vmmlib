#ifndef __VMML__MATRIX_MXM__HPP__
#define __VMML__MATRIX_MXM__HPP__

#include "matrix_mxn.hpp"

/*
*
* a subclass of the mxn matrices that adds some functionality that is only
* defined for square matrices.
*
*/

namespace vmml
{
template< size_t M, typename float_t = double >
class matrix_mxm : public matrix_mxn< M, M, float_t >
{
public:
    float_t computeDeterminant() const;
    void getAdjugate( matrix_mxm< M, float_t >& adjugate ) const;
	// we need a tolerance term since the determinant is often != 0 
	// with real numbers because of precision errors.
	void computeInverse( matrix_mxm< M, float_t >& inverse, float_t tolerance = 1e-9 ) const;

protected:

}; // class matrix_mxm



template< size_t M, typename float_t >
float_t
matrix_mxm< M, float_t >::computeDeterminant() const
{
    throw "not implemented yet.";
}


template< size_t M, typename float_t >
void
matrix_mxm< M, float_t >::computeInverse( matrix_mxm< M, float_t >& Minverse, float_t tolerance ) const
{
    throw "not implemented yet.";
}



template< size_t M, typename float_t >
void
matrix_mxm< M, float_t >::getAdjugate( matrix_mxm< M, float_t >& adjugate ) const
{
    throw "not implemented yet.";
}


// specialization for M = 2

template< typename float_t >
class matrix_mxm< 2, float_t > : public matrix_mxn< 2, 2, float_t >
{
public: 
    float_t computeDeterminant() const
    {
        const array< array< float_t, 2 >, 2 >& rows = matrix_mxn< 2, 2, float_t >::_rows;
        const float_t& a = rows[ 0 ][ 0 ];
        const float_t& b = rows[ 0 ][ 1 ];
        const float_t& c = rows[ 1 ][ 0 ];
        const float_t& d = rows[ 1 ][ 1 ];
        return a * d - b * c;
    }

    void computeInverse( matrix_mxm< 2, float_t >& Minverse, float_t tolerance = 1e-9 ) const
    {
		float_t det = computeDeterminant();
		if ( det > tolerance )
		{
			std::cerr 
				<< "error: matrix is not invertible - "
				<< "determinant is " << det << "!"
				<< std::endl;
			throw "matrix is not invertible.";
        }
		float_t reciprocal_of_determinant = 1.0 / det;
        
        // set Minverse to the adjugate of M
        getAdjugate( Minverse );
        Minverse *= reciprocal_of_determinant;
    }
    
    void getAdjugate( matrix_mxm< 2, float_t >& adjugate ) const
    {
        const array< array< float_t, 2 >, 2 >& rows = matrix_mxn< 2, 2, float_t >::_rows;
        const float_t& a = rows[ 0 ][ 0 ];
        const float_t& b = rows[ 0 ][ 1 ];
        const float_t& c = rows[ 1 ][ 0 ];
        const float_t& d = rows[ 1 ][ 1 ];

        adjugate[ 0 ][ 0 ] = d;
        adjugate[ 0 ][ 1 ] = -b;
        adjugate[ 1 ][ 0 ] = -c;
        adjugate[ 1 ][ 1 ] = a;
    }

};


// specialization for M = 3
template< typename float_t >
class matrix_mxm< 3, float_t > : public matrix_mxn< 3, 3, float_t >
{
public: 
    float_t computeDeterminant() const
    {
        const array< array< float_t, 3 >, 3 >& rows = matrix_mxn< 3, 3, float_t >::_rows;
        const float_t& a = rows[ 0 ][ 0 ];
        const float_t& b = rows[ 0 ][ 1 ];
        const float_t& c = rows[ 0 ][ 2 ];
        const float_t& d = rows[ 1 ][ 0 ];
        const float_t& e = rows[ 1 ][ 1 ];
        const float_t& f = rows[ 1 ][ 2 ];
        const float_t& g = rows[ 2 ][ 0 ];
        const float_t& h = rows[ 2 ][ 1 ];
        const float_t& i = rows[ 2 ][ 2 ];

        return ( a * e * i + b * f * g + c * d * h ) - ( g * e * c + h * f * a + i * d * b );
    }
};


} // namespace vmml

#endif

