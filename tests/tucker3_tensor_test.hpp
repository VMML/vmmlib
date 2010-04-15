#ifndef __VMML__TUCKER3_TENSOR_TEST__HPP__
#define __VMML__TUCKER3_TENSOR_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class tucker3_tensor_test : public unit_test
		{
		public:
			tucker3_tensor_test() : unit_test( "tucker3_tensor (J1xJ2xJ3) x (I1xJ1) x (I2xJ2) x (I3xJ3)" ) {}
			bool run();
			
		protected:
			
		}; // class tucker3_tensor_test
	
} // namespace vmml

#endif

