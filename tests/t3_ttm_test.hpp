#ifndef __VMML__T3_TTM_TEST__HPP__
#define __VMML__T3_TTM_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class t3_ttm_test : public unit_test
	{
	public:
		t3_ttm_test() : unit_test( "t3 TTM (tensor matrix multiplication)" ) {}
		virtual bool run();
		
	protected:
		
	}; // class t3_ttm_test
	
} // namespace vmml

#endif

