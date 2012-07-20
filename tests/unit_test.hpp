#ifndef __VMML__UNIT_TEST__HPP__
#define __VMML__UNIT_TEST__HPP__

#include <string>
#include <iostream>

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
    
    void set_prefixes(  const std::string ok_       = "[  ok  ]", 
                        const std::string& warn_    = "[ WARN ]", 
                        const std::string& fail_    = "[ FAIL ]" );
        
    const std::string& get_fail_prefix();
    const std::string& get_warn_prefix();
    const std::string& get_ok_prefix();

protected:
	virtual void log( const std::string& msg, bool status_ok, bool warning_only = false );
	virtual void log_error( const std::string& msg, bool warning_only = false );
       
	std::string _log;
    
    std::string _fail_prefix;
    std::string _warn_prefix;
    std::string _ok_prefix;

    double  _tolerance;

}; // class unit_test

} // namespace vmml

#endif

