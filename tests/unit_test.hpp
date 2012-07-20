#ifndef __VMML__UNIT_TEST__HPP__
#define __VMML__UNIT_TEST__HPP__

#include <string>
#include <iostream>

#include "unit_test_globals.hpp"

namespace vmml
{

class unit_test
{
public:
	unit_test( const std::string& test_name );
    virtual ~unit_test() {}

	virtual bool run() = 0;
	
    friend std::ostream& operator << ( std::ostream& os, 
        const unit_test& unit_test_ )
	{
		os << unit_test_._log;
		return os;
	}
    
protected:
	virtual void log( const std::string& msg, bool status_ok, bool warning_only = false );
	virtual void log_error( const std::string& msg, bool warning_only = false );
       
    unit_test_globals&  _globals;
	std::string         _log;
    double              _tolerance;

}; // class unit_test

} // namespace vmml

#endif

