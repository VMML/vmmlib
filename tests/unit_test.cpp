#include "unit_test.hpp"

namespace vmml
{

unit_test::unit_test( const std::string& test_name )
    : _tolerance( 1e-9 )
    , _ok_prefix(   "[  ok  ]" )
    , _warn_prefix( "[ WARN ]" )
    , _fail_prefix( "[ FAIL ]" )
{
	_log = "UNIT TEST: ";
	_log += test_name;
	_log += "\n";
}



void
unit_test::log( const std::string& msg, bool status_ok, bool warning_only )
{
    if ( status_ok )
        _log += _ok_prefix;
    else 
        _log += warning_only ? _warn_prefix : _fail_prefix;
    
	_log += " ";
	_log += msg;
	_log += "\n";
}



void
unit_test::log_error( const std::string& msg, bool warning_only  )
{
    log(msg, false, warning_only);
}



void
unit_test::set_prefixes( const std::string ok_, const std::string& warn_, 
    const std::string& fail_ ) 
{
    _ok_prefix      = ok_;
    _warn_prefix    = warn_;
    _fail_prefix    = fail_;
}

    
const std::string&
unit_test::get_fail_prefix()
{
    return _fail_prefix;
}


const std::string&
unit_test::get_warn_prefix()
{
    return _warn_prefix;
}


const std::string&
unit_test::get_ok_prefix()
{
    return _ok_prefix;
}

} // namespace vmml

