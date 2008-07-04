#ifndef __VMML__MATRIX_MXN_FUNCTORS__HPP__
#define __VMML__MATRIX_MXN_FUNCTORS__HPP__

#include "matrix_mxn.hpp"

namespace vmml
{

template< typename float_t = double >
struct invert_2x2_matrix_functor
{
    void operator()( 
        const matrix_mxn< 2, 2, float_t >& M, 
        matrix_mxn< 2, 2, float_t >& Minverse
    )
    {
        const float_t& a = M[ 0 ][ 0 ];
        const float_t& b = M[ 0 ][ 1 ];
        const float_t& c = M[ 1 ][ 0 ];
        const float_t& d = M[ 1 ][ 1 ];
        
        float_t reciprocal_of_determinant = 1.0 / ( a * d - b * c );
        
        // set Minverse to the adjugate of M
        Minverse[ 0 ][ 0 ] = d;
        Minverse[ 0 ][ 1 ] = -b;
        Minverse[ 1 ][ 0 ] = -c;
        Minverse[ 1 ][ 1 ] = a;
    
        Minverse *= reciprocal_of_determinant;
    }

}; // struct x2_invert_functor

} // namespace vmml

#endif

