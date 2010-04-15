#ifndef __VMML__TENSOR3_TEST__HPP__
#define __VMML__TENSOR3_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class tensor3_test : public unit_test
{
public:
	tensor3_test() : unit_test( "tensor3 (I1xI2xI3)" ) {}
    bool run();

protected:

}; // class tensor3_test

} // namespace vmml

#endif

