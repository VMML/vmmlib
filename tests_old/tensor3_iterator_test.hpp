#ifndef __VMML__TENSOR3_ITERATOR_TEST__HPP__
#define __VMML__TENSOR3_ITERATOR_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class tensor3_iterator_test : public unit_test
	{
	public:
		tensor3_iterator_test() : unit_test( "tensor3 iterator test" ) {}
		bool run();
		
	protected:
		
	}; // class tensor3_iterator_test
	
} // namespace vmml

#endif


