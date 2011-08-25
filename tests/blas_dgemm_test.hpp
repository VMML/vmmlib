#ifndef __VMML__BLAS_DGEMM_TEST__HPP__
#define __VMML__BLAS_DGEMM_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class blas_dgemm_test : public unit_test
	{
	public:
		blas_dgemm_test() : unit_test( "matrix matrix multiplication (dgemm) using blas" ) {}
		virtual bool run();
		
	protected:
		
	}; // class blas_dgemm_test
	
} // namespace vmml

#endif