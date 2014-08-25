/*
 * Copyright (c) 2006-2014, Visualization and Multimedia Lab,
 *                          University of Zurich <http://vmml.ifi.uzh.ch>,
 *                          Eyescale Software GmbH,
 *                          Blue Brain Project, EPFL
 *
 * This file is part of VMMLib <https://github.com/VMML/vmmlib/>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution.  Neither the name of the Visualization and Multimedia
 * Lab, University of Zurich nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific prior
 * written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __VMML__JACOBI_SOLVER__HPP__
#define __VMML__JACOBI_SOLVER__HPP__

#include <vmmlib/vmmlib.hpp>

#include <cmath>
#include <cassert>

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

template < typename T > 
bool solve_jacobi_3x3(
    matrix< 3, 3, T >&    a,
    vector< 3, T >&       d, 
    matrix< 3, 3, T >&    v,
    size_t& rotationCount )
{
    identity( v );

    vector< 3, T > b, z;

    for ( size_t i = 0; i < 3; ++i ) 
    {
        b[i] = d[i] = a( i,i );
        z[i] = 0.0;
    }
    T t, theta, s, c, tau; 
    size_t rot = 0;
    for ( size_t i = 1; i <= 150; ++i ) 
    {
        T sm = 0.0;
        for ( size_t ip = 0; ip < 2; ++ip ) // < n-1 
        {
            for ( size_t iq = ip + 1; iq < 3; ++iq ) // < n
            {
                sm += fabs( a( iq, ip ) );
            }
        }
        if ( sm == 0.0 ) 
        {
            rotationCount = rot;
            return true;
        }
        T tresh = ( i < 4 ) ? 0.2 * sm / 9.0 : 0.0;
    
        for ( ssize_t ip = 0; ip < 2; ++ip ) // ip < n - 1  
        {
            for ( ssize_t iq = ip + 1; iq < 3; ++iq )
            {
                T g = 100.0 * fabs( a( iq,ip ) );
                // this has to be fabs( x ) + g == fabs( x ) and NOT
                // g == 0.0 because of floating point evilness
                // ( inaccuracies when comparing (anyfloat) to 0.0 )
                if ( i > 4 
                     && fabs( d[ip] ) + g == fabs( d[ip] ) 
                     && fabs( d[iq] ) + g == fabs( d[iq] ) 
                   )
                {
                    a( iq, ip ) = 0.0;
                }
                else 
                {
                    if ( fabs( a( iq, iq ) ) > tresh ) 
                    {
                        T h = d[iq] - d[ip];
                        if ( fabs( h ) + g == fabs( h ) )
                        {
                            if ( h != 0.0 )
                                t = ( a( iq, ip ) ) / h;
                            else
                                t = 0.0;
                        } 
                        else 
                        {
                            if( a( iq, ip ) != 0.0 )
                                theta = 0.5 * h / ( a( iq, ip ) );
                            else
                                theta = 0.0;
                            t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta ) );
                            if ( theta < 0.0 ) 
                                t = -t;
                        }
                        c = 1.0 / sqrt( 1 + t * t );
                        s = t * c;
                        tau = s / ( 1.0 + c );
                        h = t * a( iq, ip );
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a( iq, ip ) = 0.0;
                        
                        for ( ssize_t j = 0; j <= ip - 1; ++j ) 
                        {
                            g = a( ip, j );
                            h = a( iq, j );
                            a( ip, j ) = g - s * ( h + g * tau );
                            a( iq, j ) = h + s * ( g - h * tau );
                        }
                        for ( ssize_t j = ip + 1; j <= iq - 1; ++j ) 
                        {
                              g = a( j, ip );
                              h = a( iq, j );
                              a( j, ip ) = g - s * ( h + g * tau );
                              a( iq, j ) = h + s * ( g - h * tau );
                        }
                        for ( size_t j = iq + 1; j < 3; ++j ) 
                        {
                              g = a( j, ip );
                              h = a( j, iq );
                              a( j, ip ) = g - s * ( h + g * tau );
                              a( j, iq ) = h + s * ( g - h * tau );
                        }
                        for ( size_t j = 0; j < 3; ++j ) 
                        {
                            g = v( ip, j );
                            h = v( iq, j );
                            v( ip, j ) = g - s * ( h + g * tau );
                            v( iq, j ) = h + s * ( g - h * tau );
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

}

} // namespace vmml

#endif

