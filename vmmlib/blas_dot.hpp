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
#ifndef __VMML__VMMLIB_BLAS_DOT__HPP__
#define __VMML__VMMLIB_BLAS_DOT__HPP__


#include <vmmlib/vector.hpp>
#include <vmmlib/exception.hpp>
#include <vmmlib/blas_includes.hpp>
#include <vmmlib/blas_types.hpp>

/**
 *
 *   a wrapper for blas's DOT routine.
 *   REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
 *   .. Scalar Arguments ..
 *   INTEGER INCX,INCY,N
 *
 *   .. Array Arguments ..
 *   REAL SX(*),SY(*)
 *
 *
 *  Purpose
 *  =======
 *
 *     SDOT forms the dot product of two vectors.
 *     uses unrolled loops for increments equal to one. *
 *
 *   more information in: http://netlib.org/blas/sdot.f
 **
 */


namespace vmml
{
namespace blas
{
    template< typename float_t >
    struct dot_params
    {
        blas_int         n;
        float_t*        x;
        blas_int        inc_x;
        float_t*        y;
        blas_int        inc_y;

        friend std::ostream& operator << ( std::ostream& os,
                                          const dot_params< float_t >& p )
        {
            os
            << " (1)\tn "     << p.n << std::endl
            << " (2)\tx "    << p.x << std::endl
            << " (3)\tincX "     << p.inc_x << std::endl
            << " (4)\ty "        << p.y << std::endl
            << " (5)\tincY "      << p.inc_y << std::endl
            << std::endl;
            return os;
        }

    };



    template< typename float_t >
    inline float_t
    dot_call( dot_params< float_t >& )
    {
        VMMLIB_ERROR( "not implemented for this type.", VMMLIB_HERE );
    }


    template<>
    inline float
    dot_call( dot_params< float >& p )
    {
        //std::cout << "calling blas sdot (single precision) " << std::endl;
        float vvi = cblas_sdot(
                         p.n,
                         p.x,
                         p.inc_x,
                         p.y,
                         p.inc_y
                         );
        return vvi;
    }

    template<>
    inline double
    dot_call( dot_params< double >& p )
    {
        //std::cout << "calling blas ddot (double precision) " << std::endl;
        double vvi = cblas_ddot(
                          p.n,
                          p.x,
                          p.inc_x,
                          p.y,
                          p.inc_y
                          );
        return vvi;
    }

} // namespace blas



template< size_t M, typename float_t >
struct blas_dot
{

    typedef vector< M, float_t > vector_t;

    blas_dot();
    ~blas_dot() {};

    bool compute( const vector_t& A_, const vector_t& B_, float_t& dot_prod_ );


    blas::dot_params< float_t > p;

    const blas::dot_params< float_t >& get_params(){ return p; };


}; // struct blas_dot


template< size_t M, typename float_t >
blas_dot< M, float_t >::blas_dot()
{
    p.n        = M;
    p.x        = 0;
    p.inc_x    = 1;
    p.y        = 0;
    p.inc_y    = 1;
}


template< size_t M, typename float_t >
bool
blas_dot< M, float_t >::compute( const vector_t& A_, const vector_t& B_, float_t& dot_prod_ )
{
    // blas needs non-const data
    vector_t* AA = new vector_t( A_ );
    vector_t* BB = new vector_t( B_ );

    p.x         = AA->array;
    p.y         = BB->array;

    dot_prod_ = blas::dot_call< float_t >( p );

    //std::cout << dot_prod_ << std::endl; //debug

    delete AA;
    delete BB;

    return true;
}

} // namespace vmml

#endif
