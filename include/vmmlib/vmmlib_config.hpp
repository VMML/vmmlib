#ifndef __VMML__VMMLIB_CONFIG__HPP__
#define __VMML__VMMLIB_CONFIG__HPP__

// #define VMMLIB_NO_SFINAE

#ifdef _MSC_VER
#  ifndef HAVE_SSIZE_T
#    include <basetsd.h>
typedef SSIZE_T ssize_t;
#    define HAVE_SSIZE_T
#  endif
#endif

// enabling this switch will have the following effect:
// operator T* will not be compiled, but for vectors,
// operator[] will instead be used. This means you can
// use vec[2] as before, but things like glVertex3fv( vec )
// will not work anymore.
//#define VMMLIB_NO_CONVERSION_OPERATORS

// enabling will align the data of vector, frustum and matrix to the given value
// Note: This won't work if other libraries in your application are also using
//       vmmlib with a different alignment. vmmlib won't check this.
//#define VMMLIB_FORCE_ALIGNMENT 16

#ifndef VMMLIB_CUSTOM_CONFIG
#  ifndef NDEBUG
#    define VMMLIB_SAFE_ACCESSORS
#  endif
#  define VMMLIB_THROW_EXCEPTIONS
#endif

#ifdef VMMLIB_FORCE_ALIGNMENT
#  ifdef __GNUC__
#    define VMMLIB_ALIGN( var ) var __attribute__((aligned(VMMLIB_FORCE_ALIGNMENT)))
#  elif defined WIN32
#    define VMMLIB_ALIGN( var ) __declspec (align (VMMLIB_FORCE_ALIGNMENT)) var
#  else
#    error "Alignment macro undefined"
#  endif
#else
#  define VMMLIB_ALIGN( var ) var
#endif

#endif
