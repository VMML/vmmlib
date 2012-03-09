#ifndef __VMML__BLAS_DAXPY_TEST__HPP__
#define __VMML__BLAS_DAXPY_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class blas_daxpy_test : public unit_test
	{
	public:
		blas_daxpy_test() : unit_test( "daxpy product using blas" ) {}
		virtual bool run();
		
	protected:
		
	}; // class blas_daxpy_test
	
} // namespace vmml

#endif

