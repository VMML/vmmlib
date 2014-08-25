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

/* @author Susanne Suter
 *
 * Tool to validate a vector, a matrix or a tensor
 * 
 */

#ifndef __VMML__VALIDATOR__HPP__
#define __VMML__VALIDATOR__HPP__

#include <vmmlib/tensor3.hpp>


namespace vmml
{
	
	class validator
	{
	public:    
		
		
		template< size_t M, typename T >
		static bool is_valid( const vector< M, T >& data_ );

		template< size_t M, size_t N, typename T >
		static bool is_valid( const matrix< M, N, T >& data_ );

		template< size_t M, size_t N, size_t L, typename T >
		static bool is_valid( const tensor3< M, N, L, T >& data_ );

		
	}; //end validator class
	
#define VMML_TEMPLATE_CLASSNAME     validator
	


template< size_t M, size_t N, size_t L, typename T >
bool
VMML_TEMPLATE_CLASSNAME::is_valid( const tensor3< M, N, L, T >& data_ )
{
    typedef typename tensor3< M, N, L, T>::const_iterator	const_t3_iterator;
	bool valid = true;
    for( const_t3_iterator it = data_.begin(); valid && it != data_.end(); ++it )
    {
        if ( std::isnan( *it ) )
            valid = false;
        if ( std::isinf( *it ) )
            valid = false;
    }
	
#ifdef VMMLIB_THROW_EXCEPTIONS
    if ( ! valid )
        VMMLIB_ERROR( "vector contains nan or inf.", VMMLIB_HERE );
#endif
	
    return valid;
	
	
}
	
	
template< size_t M, size_t N, typename T >
bool
VMML_TEMPLATE_CLASSNAME::is_valid( const matrix< M, N, T >& data_ )
{
    typedef typename matrix< M, N, T>::const_iterator	const_m_iterator;
	
    bool valid = true;
    for( const_m_iterator it = data_.begin(); valid && it != data_.end(); ++it )
    {
        if ( std::isnan( *it ) )
            valid = false;
        if ( std::isinf( *it ) )
            valid = false;
    }
	
#ifdef VMMLIB_THROW_EXCEPTIONS
    if ( ! valid )
        VMMLIB_ERROR( "matrix contains nan or inf.", VMMLIB_HERE );
#endif
	
    return valid;
	
}

template< size_t M, typename T >
bool
VMML_TEMPLATE_CLASSNAME::is_valid( const vector< M, T >& data_ )
{
    typedef typename vector< M, T>::const_iterator	const_v_iterator;
	bool valid = true;
    for( const_v_iterator it = data_.begin(); valid && it != data_.end(); ++it )
    {
        if ( std::isnan( *it ) )
            valid = false;
        if ( std::isinf( *it ) )
            valid = false;
    }
	
#ifdef VMMLIB_THROW_EXCEPTIONS
    if ( ! valid )
        VMMLIB_ERROR( "vector contains nan or inf.", VMMLIB_HERE );
#endif
	
    return valid;
	
}

#undef VMML_TEMPLATE_CLASSNAME
	
	
} // namespace vmml



#endif
	
	
