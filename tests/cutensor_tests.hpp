#ifndef __VMML__T3_CUTENSOR_TESTS__HPP__
#define __VMML__T3_CUTENSOR_TESTS__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class cutensor_tests : public unit_test
	{
	public:
		cutensor_tests() : unit_test( "CUDA tensor3 classes tests" ) {}
		virtual bool run();
		
	protected:
		
	}; // class cutensor_tests
	
} // namespace vmml

#endif

