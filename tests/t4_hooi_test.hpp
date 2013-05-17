#ifndef __VMML__T4_HOOI_TEST__HPP__
#define __VMML__T4_HOOI_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class t4_hooi_test : public unit_test
	{
	public:
		t4_hooi_test() : unit_test( "tensor4 HOOI" ) {}
		bool run();
		
	protected:
		
	}; // class t4_hooi_test
	
} // namespace vmml

#endif

