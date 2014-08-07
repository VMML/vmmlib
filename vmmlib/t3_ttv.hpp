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

/* @author Rafa Ballester
 *
 * Tensor times vector multiplication for tensor3 (t3)
 * Only mode 1 is implemented. This is used for fast inner product calculation. The basic idea is that the inner product between a tensor X and a 1-rank tensor (expressed as the outer product of three vectors u, v and w) is the same as ((X x u).v).w .
 *
 */

#ifndef __VMML__T3_TTV__HPP__
#define __VMML__T3_TTV__HPP__

#include <vmmlib/tensor3.hpp>

#ifdef VMMLIB_USE_OPENMP
#  include <omp.h>
#endif

namespace vmml
{
    class t3_ttv
    {
    public:

        template< size_t I1, size_t I2, size_t I3, typename T >
        static void multiply_first_mode(const tensor3< I1, I2, I3, T >& t3_in_, const vector< I1, T >& u, matrix< I2, I3, T >& m_res_);

    protected:
    };

    template< size_t I1, size_t I2, size_t I3, typename T >
    void t3_ttv::multiply_first_mode(const tensor3< I1, I2, I3, T >& t3_in_,
                                     const vector< I1, T >& u,
                                     matrix< I2, I3, T >& m_res_ )
    {
        for( unsigned j = 0; j < I2; ++j)
        {
            for( unsigned k = 0; k < I3; ++k)
            {
                T vtv = 0;
                for( unsigned i = 0; i < I1; ++i)
                    vtv += t3_in_.at(i, j, k) * u.at(i);

                m_res_.at(j, k) = vtv;
            }
        }

    }

} //end vmml namespace

#endif

