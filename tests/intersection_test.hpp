#ifndef __VMML__INTERSECTION_TEST__HPP__
#define __VMML__INTERSECTION_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class intersection_test : public unit_test
{
public:
    intersection_test() : unit_test( "intersection" ) {}
    bool run();

protected:

}; // class intersection_test

} // namespace vmml

#endif
