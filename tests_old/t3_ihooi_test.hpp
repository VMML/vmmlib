#ifndef __VMML__T3_IHOOI_TEST__HPP__
#define __VMML__T3_IHOOI_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class t3_ihooi_test : public unit_test
	{
	public:
		t3_ihooi_test() : unit_test( "tensor3 IHOOI" ) {}
		bool run();
		
	protected:
		
	}; // class t3_hooi_test
	
} // namespace vmml

#endif


