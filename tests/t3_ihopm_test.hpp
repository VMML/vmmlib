#ifndef __VMML__T3_IHOPM_TEST__HPP__
#define __VMML__T3_IHOPM_TEST__HPP__

#include "unit_test.hpp"
#include <string>

namespace vmml
{
	
	class t3_ihopm_test : public unit_test
	{
	public:
		t3_ihopm_test() : unit_test( "tensor3 incremental rank-r CP-ALS (or HOPM)" ) {}
		bool run();
		
	protected:
		
	}; // class t3_ihopmtest
	
} // namespace vmml

#endif


