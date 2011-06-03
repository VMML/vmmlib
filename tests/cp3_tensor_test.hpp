#ifndef __VMML__CP3_TENSOR_TEST__HPP__
#define __VMML__CP3_TENSOR_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class cp3_tensor_test : public unit_test
	{
	public:
		cp3_tensor_test() : unit_test( "cp3_tensor (R) x (I1xR) x (I2xR) x (I3xR)" ) {}
		bool run();
		
	protected:
		
	}; // class cp3_tensor_test
	
} // namespace vmml

#endif


