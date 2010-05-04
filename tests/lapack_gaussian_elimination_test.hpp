#ifndef __VMML__LAPACK_GAUSSIAN_ELIMINATION_TEST__HPP__
#define __VMML__LAPACK_GAUSSIAN_ELIMINATION_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{

class lapack_gaussian_elimination_test : public unit_test
{
public:
	lapack_gaussian_elimination_test() : unit_test( "gaussian elimination using lapack" ) {}
    virtual bool run();

protected:

}; // class lapack_linear_least_squares_test

} // namespace vmml

#endif

