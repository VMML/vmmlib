#ifndef __VMML__LOWPASS_FILTER_TEST__HPP__
#define __VMML__LOWPASS_FILTER_TEST__HPP__

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class lowpass_filter_test : public unit_test
{
public:
    lowpass_filter_test() : unit_test( "lowpass_filter" ) {}
    bool run();

protected:

}; // class lowpass_filter_test

} // namespace vmml

#endif
