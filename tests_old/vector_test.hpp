#ifndef VMMLIB__VECTOR_TEST__HPP
#define VMMLIB__VECTOR_TEST__HPP

#include "unit_test.hpp"

#include <string>

namespace vmml
{

class vector_test : public unit_test
{
public:
	vector_test() : unit_test( "vector (m)" ) {}
    bool run();

protected:

}; // class vector_test

} // namespace vmml

#endif

