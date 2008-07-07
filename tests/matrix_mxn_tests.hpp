#ifndef __VMML__MATRIX_MXN_TESTS__HPP__
#define __VMML__MATRIX_MXN_TESTS__HPP__

#include "matrix_mxn.hpp"
#include "matrix_mxm.hpp"
#include "matrix_mxn_functors.hpp"

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class matrix_mxn_tests : public unit_test
{
public:
	matrix_mxn_tests() : unit_test( "matrix_mxn" ) {}
    bool run();

protected:

}; // class matrix_mxn_tests

} // namespace vmml

#endif

