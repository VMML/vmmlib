#ifndef __VMML__FRUSTUM_TEST__HPP__
#define __VMML__FRUSTUM_TEST__HPP__

#include "unit_test.hpp"

#include <vmmlib/frustum_culler.hpp>

namespace vmml
{

class frustum_test : public unit_test
{
public:
    frustum_test() : unit_test( "frustum" ) {}
    bool run();

private:
    bool _test( const frustum_culler< float > fc );
};

}

#endif

