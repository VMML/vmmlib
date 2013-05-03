#ifndef __VMML__UTIL_TEST__HPP__
#define __VMML__UTIL_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class util_test : public unit_test
{
public:
	util_test() : unit_test( "util (m)" ) {}
    bool run();

protected:

}; // class util_test

} // namespace vmml

#endif

