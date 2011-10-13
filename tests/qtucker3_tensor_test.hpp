#ifndef __VMML__QTUCKER3_TENSOR_TEST__HPP__
#define __VMML__QTUCKER3_TENSOR_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class qtucker3_tensor_test : public unit_test
	{
	public:
		qtucker3_tensor_test() : unit_test( "quantized tucker3_tensor (R1xR2xR3) x (I1xR1) x (I2xR2) x (I3xR3)" ) {}
		bool run();
		
	protected:
		
	}; // class qtucker3_tensor_test
	
} // namespace vmml

#endif

