#ifndef __VMML__FRUSTUM_TEST__HPP__
#define __VMML__FRUSTUM_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{

class frustum_test : unit_test
{
public:
    frustum_test() : unit_test( "frustum" ) {}
    bool run();
protected:

}; // class frustum_test

} // namespace vmml

#endif

