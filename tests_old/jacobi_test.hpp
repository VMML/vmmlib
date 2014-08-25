#ifndef __VMML__JACOBI_TEST__HPP__
#define __VMML__JACOBI_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{

class jacobi_test : public unit_test
{
public:
    jacobi_test() : unit_test( "jacobi test" ) {}

    bool run();

protected:

}; // class jacobi_test

} // namespace vmml

#endif

