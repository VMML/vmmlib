#ifndef __VMML__BLAS_DOT_TEST__HPP__
#define __VMML__BLAS_DOT_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class blas_dot_test : public unit_test
	{
	public:
		blas_dot_test() : unit_test( "dot product (dot) using blas" ) {}
		virtual bool run();
		
	protected:
		
	}; // class blas_dot_test
	
} // namespace vmml

#endif

