#ifndef __VMML__VMMLIB_LAPACK_LINEAR_LEAST_SQUARES__HPP__
#define __VMML__VMMLIB_LAPACK_LINEAR_LEAST_SQUARES__HPP__

#include <vmmlib/matrix.hpp>
#include <vmmlib/vector.hpp>
#include <vmmlib/exception.hpp>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
// FIXME - include clapack headers
#endif

#include <string>

namespace vmml
{

// XYYZZZ 
// X    = data type: S - float, D - double
// YY   = matrix type, GE - general, TR - triangular
// ZZZ  = function name

namespace lapack
{

template< typename float_t >
struct llsq_params
{
    __CLPK_integer  n; // order of matrix A = M * N
    __CLPK_integer  nrhs; // number of columns of B
    float_t*        a;   // input A, output P*L*U
    __CLPK_integer  lda; // leading dimension of A (for us: number of rows)
    __CLPK_integer* ipiv; // pivot indices, integer array of size N
    float_t*        b;  // input b, output X
    __CLPK_integer  ldb; // leading dimension of b
    __CLPK_integer  info;
    
    friend std::ostream& operator << ( std::ostream& os, 
        const llsq_params< float_t >& p )
    {
        os 
            << "n "         << p.n 
            << " nrhs "     << p.nrhs
            << " lda "      << p.lda
            << " ldb "      << p.ldvt 
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
llsq_call( llsq_params< float_t >& p )
{
    VMMLIB_ERROR( "not implemented for this type.", VMMLIB_HERE );
}


template<>
inline void
llsq_call( llsq_params< float >& p )
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
llsq_call( llsq_params< double >& p )
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



} // namespace lapack

template< size_t M, size_t N, typename float_t >
struct lapack_linear_least_squares
{
    void compute(
        matrix< N, N, float_t >& A, 
        matrix< N, M, float_t >& b 
        )
    {
        p.a = A.array;
        p.b = b.array;
        
        lapack::llsq_call( p );

        if ( p.info != 0 )
        {
            if ( p.info < 0 )
                VMMLIB_ERROR( "invalid value in input matrix", VMMLIB_HERE );
            else
                VMMLIB_ERROR( "factor U is exactly singular, solution could not be computed.", VMMLIB_HERE );
        }
    }


    lapack_linear_least_squares()
    {
        p.n     = N;
        p.nrhs  = M;
        p.lda   = N;
        p.ldb   = N;
        p.ipiv = new __CLPK_integer[ N ];
    
    }
    
    ~lapack_linear_least_squares()
    {
        delete[] p.ipiv;
    }

    lapack::llsq_params< float_t > p;
    
}; // struct lapack_linear_least_squares



} // namespace vmml

#endif

