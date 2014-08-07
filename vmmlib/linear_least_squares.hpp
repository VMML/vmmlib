/*
 * Copyright (c) 2006-2014, Visualization and Multimedia Lab,
 *                          University of Zurich <http://vmml.ifi.uzh.ch>;
 *                          Eyescale Software GmbH;
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
#ifndef __VMML__LINEAR_LEAST_SQUARES__HPP__
#define __VMML__LINEAR_LEAST_SQUARES__HPP__

#include <vmmlib/vector3.h>
#include "matrix_mxn.hpp"

namespace vmml
{

// beta_hat = ( Xt * X )^-1 * Xt * y;
// beta_hat -> best-fit 'solutions' to a,b = [ a,b ]transposed.

// TODO FIXME what about the constant???

struct linear_least_squares
{

template< size_t number_of_data_points, size_t number_of_unknowns, typename float_t >
vmml::matrix_mxn< number_of_unknowns, 1, float_t >
solve( std::vector< Vector3< float_t > >& data_points, float_t tolerance = 1e-9 )
{
    
	assert( data_points.size() >= number_of_data_points );
	
    matrix_mxn< number_of_data_points, number_of_unknowns, float_t > X;
    matrix_mxn< number_of_data_points, 1, float_t > y;

    for( size_t nb_index = 0; nb_index < number_of_data_points; ++nb_index )
    {
		float_t x = data_points[ nb_index ].x;

		for( size_t index = 0; index < number_of_unknowns; ++index )
		{
			X[ nb_index ][ index ] = x;
			x *= x;			
		}

        y[ nb_index ][ 0 ] = data_points[ nb_index ].y;
    }

    std::cout << X << std::endl;
    
    vmml::matrix_mxn< number_of_unknowns, number_of_data_points > Xt;
    X.transposeTo( Xt );
    
    std::cout << Xt << std::endl;
    
    vmml::matrix_mxm< number_of_unknowns > XtX;
    XtX.mul( Xt, X );
    
    vmml::matrix_mxm< number_of_unknowns > XtX_inverse;
	XtX.computeInverse( XtX_inverse, tolerance );

    vmml::matrix_mxn< number_of_unknowns, number_of_data_points > XtX_inv_mul_Xt;
    
    XtX_inv_mul_Xt.mul( XtX_inverse, Xt );
    
    vmml::matrix_mxn< number_of_unknowns, 1 > beta_hat;
    beta_hat.mul( XtX_inv_mul_Xt, y );

    std::cout << beta_hat << std::endl;

	return beta_hat;
}



template< size_t number_of_data_points, size_t number_of_unknowns, typename float_t >
vmml::matrix_mxn< number_of_unknowns, 1, float_t >
solve_using_qr_decomposition( std::vector< Vector3< float_t > >& data_points, float_t tolerance = 1e-9 )
{

	assert( data_points.size() >= number_of_data_points );
	
    matrix_mxn< number_of_data_points, number_of_unknowns, float_t > X;
    matrix_mxn< number_of_data_points, 1, float_t > y;

	// r = y - X * Beta
	// X = Q*R
	// Qt * r = Qt * y - ( Qt * Q ) R * Beta = 
	//
	// | ( Qt * y )_n - R_n * Beta |   | U |
	// | ( Qt * y )_m-n            | = | L |
	// 
	// S = rt * Q * Qt * r = rt * r
	//
	// S = Ut * U + Lt * L
	//
	// => R_n * BetaHat = ( Qt * y )_n
	//


    vmml::matrix_mxn< number_of_unknowns, 1 > beta_hat;
	return beta_hat;
}


template< size_t number_of_data_points, size_t number_of_unknowns, typename float_t >
vmml::matrix_mxn< number_of_unknowns, 1, float_t >
solve_using_svd( std::vector< Vector3< float_t > >& data_points, float_t tolerance = 1e-9 )
{


}

};


} // namespace vmml

#endif

