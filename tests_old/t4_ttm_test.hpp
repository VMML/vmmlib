#ifndef __VMML__T4_TTM_TEST__HPP__
#define __VMML__T4_TTM_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class t4_ttm_test : public unit_test
	{
	public:
		t4_ttm_test() : unit_test( "tensor4 TTM (tensor matrix multiplication)" ) {}
		virtual bool run();
		
	protected:
		
	}; // class t4_ttm_test
	
} // namespace vmml

#endif

