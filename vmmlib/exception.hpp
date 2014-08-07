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
#ifndef __VMMLIB__EXCEPTION__HPP__
#define __VMMLIB__EXCEPTION__HPP__

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cassert>

#include <vmmlib/vmmlib_config.hpp>


#define VMMLIB_HERE ( except_here( __FILE__, __LINE__ ) )

#ifdef VMMLIB_THROW_EXCEPTIONS
#define VMMLIB_ERROR( desc, here ) throw( exception( desc, here ) )
#else
#define VMMLIB_ERROR( desc, here ) error_noexcept( desc, here )
#endif


namespace vmml
{

struct except_here
{
    except_here( const char* file_, int line_ ) : file( file_ ), line( line_ ) {}
    const char*     file;
    int             line;  
}; // struct except_here


inline void
error_noexcept( const std::string& desc, const except_here& here )
{
    std::cerr 
        << "vmmlib error at " << here.file << ":" << here.line << "\n"
        << desc 
        << std::endl;
    assert( 0 );
}



class exception : public std::exception
{
public:
    exception( const std::string& desc, const except_here& here )
    : _description( desc )
    , _here( here )
    {}
    
    virtual ~exception() throw() {}
    
    virtual const char* what() const throw()
    {
        std::stringstream ss;
        ss
            << _here.file << "(" << _here.line << "): - " 
            << _description 
            << std::endl;
        return ss.str().c_str();
    }
    
protected:
    std::string         _description;
    const except_here&  _here;

private:
    // disallow std ctor
    exception() : _here( *new except_here( "", 0 ) ){};
    // disallow assignment operator
    virtual const exception& operator=( const exception& ){ return *this; }


};


} // namespace stream_process

#endif

