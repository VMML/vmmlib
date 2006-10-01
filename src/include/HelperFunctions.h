/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Jonas Boesch
* @author Stefan Eilemann
* @author Renato Pajarola
* @author David H. Eberly ( Wild Magic )
* @author Andrew Willmott ( VL )
*
* @license BSD license, check LICENSE
*
* parts of the source code of VMMLib were inspired by David Eberly's 
* Wild Magic and Andrew Willmott's VL.
* 
*/ 


/*
*
* a collection of small helper functions 
*
*/

#include "VMMLib.h"

namespace vmml
{

template< class Real >
inline Real min( const Real a, const Real b )
{
    return ( a < b ) ? a : b;
}

template< class Real >
inline Real max( const Real a, const Real b )
{
    return ( a > b ) ? a : b;
}

template< class Real >
inline Real squared( const Real a )
{
    return ( a == 0.0 ) ? 0.0 : a * a;
}

/*
 * Computes (a2 + b2)1/2 without destructive underflow or overflow.
 */
template< class Real >
inline Real pythag( Real a, Real b )
{
    a = fabs(a);
    b = fabs(b);
    if ( a > b )
        return a * sqrt( 1.0 + squared( b / a ) );
    else
        return ( b == 0.0 ) ? 0.0 : b * sqrt( 1.0 + squared( a / b ) );
}

template< class Real >
inline Real sign( Real a, Real b )
{
    return ( b >= 0.0 ) ? fabs( a ) : -fabs( a );
}


}; // namespace vmml