#ifndef __VMML__MATRIX_PSEUDOINVERSE_TEST__HPP__
#define __VMML__MATRIX_PSEUDOINVERSE_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{
	
	class matrix_pseudoinverse_test : public unit_test
	{
	public:
		matrix_pseudoinverse_test() : unit_test( "matrix pseudoinverse test" ) {}
		bool run();
		
	protected:
		
	}; // class matrix_pseudoinverse_test
	
} // namespace vmml

#endif


