#ifndef __VMML__LOPASS_FILTER_TEST__HPP__
#define __VMML__LOPASS_FILTER_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class lopass_filter_test : public unit_test
{
public:
    lopass_filter_test() : unit_test( "lopass_filter" ) {}
    bool run();

protected:

}; // class lopass_filter_test

} // namespace vmml

#endif

