#ifndef __VMML__TUCKER3_TENSOR_TEST__HPP__
#define __VMML__TUCKER3_TENSOR_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class tucker3_tensor_test : public unit_test
		{
		public:
			tucker3_tensor_test() : unit_test( "tucker3_tensor (R1xR2xR3) x (I1xR1) x (I2xR2) x (I3xR3)" ) {}
			bool run();
			
		protected:
			
		}; // class tucker3_tensor_test
	
} // namespace vmml

#endif

