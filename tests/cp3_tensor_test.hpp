#ifndef __VMML__CP3_TENSOR_TEST__HPP__
#define __VMML__CP3_TENSOR_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class cp3_tensor_test : public unit_test
	{
	public:
		cp3_tensor_test() : unit_test( "cp3 tensor test" ) {}
		bool run();
		
	protected:
		
	}; // class cp3_tensor_test
	
} // namespace vmml

#endif


