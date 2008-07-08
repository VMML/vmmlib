#ifndef __VMML__QR_DECOMPOSITION__HPP__
#define __VMML__QR_DECOMPOSITION__HPP__

#include <vmmlib/matrix.hpp>

#include <cmath>

/*
* QR decomposition, using the modified (stabilized) Gram-Schmidt process.
*
* A  -> matrix to be factorized
* Q  -> orthonormal
* Rn -> upper triangular 
*/

namespace vmml
{

template< size_t M, size_t N, typename float_t >
void qr_decompose( 
    const matrix< M, N, float_t >& A_, 
    matrix< M, M, float_t >& Q, 
    matrix< N, N, float_t >& Rn 
    )
{
    Q   = 0.0; 
    Rn  = 0.0;

    matrix< M, N, float_t > A( A_ );

    matrix< M, 1, float_t > column;
    matrix< M, 1, float_t > q_column;
    for( size_t k = 0; k < N; ++k )
    {
        A.getColumn( k, column );
        float_t rkk = 0.0;
        for( size_t i = 0; i < M; ++i )
        {
            rkk += column.at( i, 0 ) * column.at( i, 0 );
        }
        rkk = sqrt( rkk );
        Rn.at( k, k ) = rkk;
        
        if ( rkk == 0.0 )
            break;
        
        for( size_t i = 0; i < M; ++i )
        {
            Q.at( i, k ) = A.at( i, k ) / rkk;
        }
        
        for( size_t j = k+1; j < N; ++j )
        {
            Q.getColumn( k, q_column );
            A.getColumn( j, column );
            float_t dot = 0.0;
            for( size_t i = 0; i < M; ++i )
            {
                dot += column.at( i, 0 ) * q_column.at( i, 0 );
            }
            Rn.at( k, j ) = dot;
            
            for( size_t i = 0; i < M; ++i )
            {
                column.at( i, 0 ) -= q_column.at( i, 0 ) * Rn.at( k, j );
            }
            
            A.setColumn( j, column );
        }
    }
}

} // namespace vmml

#endif

