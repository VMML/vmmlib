#ifndef __VMML__T4_CONV_TEST__HPP__
#define __VMML__T4_CONV_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class t4_converter_test : public unit_test
	{
	public:
		t4_converter_test() : unit_test( "t4_converter" ) {}
		virtual bool run();
		
	protected:
		
	}; // class t4_converter_test
	
} // namespace vmml

#endif

