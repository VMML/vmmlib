#ifndef __VMML__CUBLAS_DGEMM_TEST__HPP__
#define __VMML__CUBLAS_DGEMM_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class cublas_dgemm_test : public unit_test
	{
	public:
		cublas_dgemm_test();
		
		virtual bool run();
		
	protected:
		
	}; // class cublas_dgemm_test
	
} // namespace vmml

#endif

