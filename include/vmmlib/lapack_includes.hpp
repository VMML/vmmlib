#ifndef __VMML__LAPACK_INCLUDES__HPP__
#define __VMML__LAPACK_INCLUDES__HPP__


#ifdef __APPLE__

#include <Accelerate/Accelerate.h>

#else

// FIXME - include clapack headers
extern "C" {
#include <f2c.h>
#include <clapack/clapack.h>
}


#endif



#endif

