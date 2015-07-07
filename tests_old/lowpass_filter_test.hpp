#ifndef VMMLIB__LOWPASS_FILTER_TEST__HPP
#define VMMLIB__LOWPASS_FILTER_TEST__HPP

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
