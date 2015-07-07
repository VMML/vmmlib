#ifndef VMMLIB__QR_DECOMPOSITION_TEST__HPP
#define VMMLIB__QR_DECOMPOSITION_TEST__HPP

#include "unit_test.hpp"

namespace vmml
{

class qr_decomposition_test : public unit_test
{
public:
	qr_decomposition_test() : unit_test( "QR decomposition" ) {}

    bool run();
protected:

}; // class qr_decomposition_test

} // namespace vmml

#endif

