#ifndef __VMML__T3_VIRTUAL_PADDER_TEST__HPP__
#define __VMML__T3_VIRTUAL_PADDER_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class t3_virtual_padder_test : public unit_test
	{
	public:
		t3_virtual_padder_test() : unit_test( "t3_virtual_padder" ) {}
		virtual bool run();
		
	protected:
		
	}; // class t3_virtual_padder_test
	
} // namespace vmml

#endif

