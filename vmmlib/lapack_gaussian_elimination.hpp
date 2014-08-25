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
#ifndef __VMML__VMMLIB_LAPACK_GAUSSIAN_ELIMINATION__HPP__
#define __VMML__VMMLIB_LAPACK_GAUSSIAN_ELIMINATION__HPP__

#include <vmmlib/matrix.hpp>
#include <vmmlib/vector.hpp>
#include <vmmlib/exception.hpp>

#include <vmmlib/lapack_types.hpp>
#include <vmmlib/lapack_includes.hpp>

#include <string>

/**
*
* this is a wrapper for the following lapack routines:
*
* xGESV
*
*
*/


namespace vmml
{

// XYYZZZ
// X    = data type: S - float, D - double
// YY   = matrix type, GE - general, TR - triangular
// ZZZ  = function name

namespace lapack
{

//
//
//
// SGESV/DGESV
//
//

template< typename float_t >
struct xgesv_params
{
    lapack_int      n; // order of matrix A = M * N
    lapack_int      nrhs; // number of columns of B
    float_t*        a;   // input A, output P*L*U
    lapack_int      lda; // leading dimension of A (for us: number of rows)
    lapack_int*     ipiv; // pivot indices, integer array of size N
    float_t*        b;  // input b, output X
    lapack_int      ldb; // leading dimension of b
    lapack_int      info;

    friend std::ostream& operator << ( std::ostream& os,
        const xgesv_params< float_t >& p )
    {
        os
            << "n "         << p.n
            << " nrhs "     << p.nrhs
            << " lda "      << p.lda
            << " ldb "      << p.ldb
            << " info "     << p.info
            << std::endl;
        return os;
    }

};


#if 0
/* Subroutine */ int dgesv_(integer *n, integer *nrhs, doublereal *a, integer
	*lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
#endif


template< typename float_t >
inline void
xgesv_call( xgesv_params< float_t >& )
{
    VMMLIB_ERROR( "not implemented for this type.", VMMLIB_HERE );
}


template<>
inline void
xgesv_call( xgesv_params< float >& p )
{
    sgesv_(
        &p.n,
        &p.nrhs,
        p.a,
        &p.lda,
        p.ipiv,
        p.b,
        &p.ldb,
        &p.info
    );

}


template<>
inline void
xgesv_call( xgesv_params< double >& p )
{
    dgesv_(
        &p.n,
        &p.nrhs,
        p.a,
        &p.lda,
        p.ipiv,
        p.b,
        &p.ldb,
        &p.info
    );
}


template< size_t M, size_t N, typename float_t >
struct gaussian_elimination
{
    // computes x ( Ax = b ). x replaces b on output.
    void compute(
        matrix< N, N, float_t >& A,
        matrix< N, M, float_t >& b
        );

    void compute(
        matrix< N, N, float_t >& A,
        vector< N, float_t >& b
        );

    gaussian_elimination();
    ~gaussian_elimination();

    const lapack::xgesv_params< float_t >& get_params() { return p; }

    lapack::xgesv_params< float_t > p;

}; // struct lapack_linear_least_squares



template< size_t M, size_t N, typename float_t >
void
gaussian_elimination< M, N, float_t >::
compute(
        matrix< N, N, float_t >& A,
        matrix< N, M, float_t >& b
        )
{
    p.a = A.array;
    p.b = b.array;

    lapack::xgesv_call( p );

    if ( p.info != 0 )
    {
        if ( p.info < 0 )
            VMMLIB_ERROR( "invalid value in input matrix", VMMLIB_HERE );
        else
            VMMLIB_ERROR( "factor U is exactly singular, solution could not be computed.", VMMLIB_HERE );
    }
}



template< size_t M, size_t N, typename float_t >
void
gaussian_elimination< M, N, float_t >::
compute(
        matrix< N, N, float_t >& A,
        vector< N, float_t >& b
        )
{
    p.a = A.array;
    p.b = b.array;

    lapack::xgesv_call( p );

    if ( p.info != 0 )
    {
        if ( p.info < 0 )
            VMMLIB_ERROR( "invalid value in input matrix", VMMLIB_HERE );
        else
            VMMLIB_ERROR( "factor U is exactly singular, solution could not be computed.", VMMLIB_HERE );
    }
}



template< size_t M, size_t N, typename float_t >
gaussian_elimination< M, N, float_t >::
gaussian_elimination()
{
    p.n     = N;
    p.nrhs  = M;
    p.lda   = N;
    p.ldb   = N;
    p.ipiv = new lapack_int[ N ];

}



template< size_t M, size_t N, typename float_t >
gaussian_elimination< M, N, float_t >::
~gaussian_elimination()
{
    delete[] p.ipiv;
}


} // namespace lapack

} // namespace vmml

#endif
