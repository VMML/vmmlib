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

#ifndef _JacobiSolver_H_
#define _JacobiSolver_H_

#include "VMMLib.h"
#include <math.h>
#include <assert.h>

namespace vmml
{

/*
 * This function computes the eigenvalues and eigenvectors of a 3x3 matrix.
 *
 * @param a matrix to be diagonalized.
 * @param d eigenvalues of A.
 * @param v matrix whose columns are the normalized eigenvectors of A.
 * @param rotationCount number of Jacobi rotations required.
 * @return true if the transformation has been done. False if not.
 *
 *
 * modified from numerical recipies for n=3 and float values
 *
 */

template < class Real > 
bool solveJacobi3x3( Matrix3< Real >& a, Vector3< Real >& d, 
                  Matrix3< Real >& v, size_t& rotationCount )
{
    v = Matrix3< Real >::IDENTITY;

    Vector3< Real > b, z;

    for ( size_t i = 0; i < 3; ++i ) 
    {
        b[i] = d[i] = a[i][i];
        z[i] = 0.0;
    }
    Real t, theta, s, c, tau; 
    size_t rot = 0;
    for ( size_t i = 1; i <= 150; ++i ) 
    {
        Real sm = 0.0;
        for ( size_t ip = 0; ip < 2; ++ip ) // < n-1 
        {
            for ( size_t iq = ip + 1; iq < 3; ++iq ) // < n
            {
                sm += fabs( a[ip][iq] );
            }
        }
        if ( sm == 0.0 ) 
        {
            rotationCount = rot;
            return true;
        }
        Real tresh = ( i < 4 ) ? 0.2 * sm / 9.0 : 0.0;
    
        for ( ssize_t ip = 0; ip < 2; ++ip ) // ip < n - 1  
        {
            for ( ssize_t iq = ip + 1; iq < 3; ++iq )
            {
                Real g = 100.0 * fabs( a[ip][iq] );
                // this has to be fabs( x ) + g == fabs( x ) and NOT
                // g == 0.0 because of floating point evilness
                // ( inaccuracies when comparing (anyfloat) to 0.0 )
                if ( i > 4 
                     && fabs( d[ip] ) + g == fabs( d[ip] ) 
                     && fabs( d[iq] ) + g == fabs( d[iq] ) 
                   )
                {
                    a[ip][iq] = 0.0;
                }
                else 
                {
                    if ( fabs( a[ip][iq] ) > tresh ) 
                    {
                        Real h = d[iq] - d[ip];
                        if ( fabs( h ) + g == fabs( h ) )
                        {
                            if ( h != 0.0 )
                                t = ( a[ip][iq] ) / h;
                            else
                                t = 0.0;
                        } 
                        else 
                        {
                            if( a[ip][iq] != 0.0 )
                                theta = 0.5 * h / ( a[ip][iq] );
                            else
                                theta = 0.0;
                            t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta ) );
                            if ( theta < 0.0 ) 
                                t = -t;
                        }
                        c = 1.0 / sqrt( 1 + t * t );
                        s = t * c;
                        tau = s / ( 1.0 + c );
                        h = t * a[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a[ip][iq] = 0.0;
                        
                        for ( ssize_t j = 0; j <= ip - 1; ++j ) 
                        {
                            g = a[j][ip];
                            h = a[j][iq];
                            a[j][ip] = g - s * ( h + g * tau );
                            a[j][iq] = h + s * ( g - h * tau );
                        }
                        for ( ssize_t j = ip + 1; j <= iq - 1; ++j ) 
                        {
                              g = a[ip][j];
                              h = a[j][iq];
                              a[ip][j] = g - s * ( h + g * tau );
                              a[j][iq] = h + s * ( g - h * tau );
                        }
                        for ( size_t j = iq + 1; j < 3; ++j ) 
                        {
                              g = a[ip][j];
                              h = a[iq][j];
                              a[ip][j] = g - s * ( h + g * tau );
                              a[iq][j] = h + s * ( g - h * tau );
                        }
                        for ( size_t j = 0; j < 3; ++j ) 
                        {
                            g = v[j][ip];
                            h = v[j][iq];
                            v[j][ip] = g - s * ( h + g * tau );
                            v[j][iq] = h + s * ( g - h * tau );
                        }
                        ++rot;
                    }
                }
            }
        }
                    
        for ( size_t ip = 0; ip < 3; ++ip ) 
        {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
  }
  return false;

};

}; // namespace vmml

#endif
