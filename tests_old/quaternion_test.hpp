#ifndef VMML_QUATERNION_TEST__HPP
#define VMML_QUATERNION_TEST__HPP

#include "unit_test.hpp"

namespace vmml
{
class quaternion_test : public unit_test
{
public: 
	quaternion_test() : unit_test( "quaternion" ) {}
    virtual bool run();
};

}
#endif

