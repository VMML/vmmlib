#ifndef __VMML__TENSOR4_TEST__HPP__
#define __VMML__TENSOR4_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class tensor4_test : public unit_test
	{
	public:
		tensor4_test() : unit_test( "tensor4" ) {}
		virtual bool run();
		
	protected:
		
	}; // class tensor4_test
	
} // namespace vmml

#endif

