#ifndef __VMML__T3_HOOI_TEST__HPP__
#define __VMML__T3_HOOI_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class t3_hooi_test : public unit_test
	{
	public:
		t3_hooi_test() : unit_test( "tensor3 HOOI" ) {}
		bool run();
		
	protected:
		
	}; // class t3_hooi_test
	
} // namespace vmml

#endif


