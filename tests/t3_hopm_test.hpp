#ifndef __VMML__T3_HOPM_TEST__HPP__
#define __VMML__T3_HOPM_TEST__HPP__

#include "unit_test.hpp"
#include <string>

namespace vmml
{
	
	class t3_hopm_test : public unit_test
	{
	public:
		t3_hopm_test() : unit_test( "tensor3 HOPM" ) {}
		bool run();
		
	protected:
		
	}; // class t3_hopmtest
	
} // namespace vmml

#endif


