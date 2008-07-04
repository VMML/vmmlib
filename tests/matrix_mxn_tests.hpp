#ifndef __VMML__MATRIX_MXN_TESTS__HPP__
#define __VMML__MATRIX_MXN_TESTS__HPP__

#include "matrix_mxn.hpp"
#include "matrix_mxm.hpp"
#include "matrix_mxn_functors.hpp"

#include <string>

namespace vmml
{

class matrix_mxn_tests
{
public:
    bool run();

    std::string error_string;
protected:
    void _passed( const std::string& testname );
}; // class matrix_mxn_tests

} // namespace vmml

#endif

