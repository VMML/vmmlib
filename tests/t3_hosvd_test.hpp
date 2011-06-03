#ifndef __VMML__T3_HOSVD_TEST__HPP__
#define __VMML__T3_HOSVD_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class t3_hosvd_test : public unit_test
	{
	public:
		t3_hosvd_test() : unit_test( "tensor3 HOSVD" ) {}
		bool run();
		
	protected:
		
	}; // class t3_hosvd_test
	
} // namespace vmml

#endif


