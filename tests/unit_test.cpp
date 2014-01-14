#include "unit_test.hpp"

namespace vmml
{

unit_test::unit_test( const std::string& test_name )
    : _globals( unit_test_globals::get_instance() )
    , _log()
    , _tolerance( 1e-9 )
{
	_log = "UNIT TEST: ";
	_log += test_name;
	_log += "\n";
}

void unit_test::log( const std::string& msg, bool status_ok, bool warning_only )
{
    _globals.notify_test();

    if ( status_ok )
        _log += _globals.get_ok_prefix();
    else
    {
        if ( warning_only )
        {
            _log += _globals.get_warn_prefix();
            _globals.notify_warn();
        }
        else
        {
            _log += _globals.get_fail_prefix();
            _globals.notify_fail();
        }
    }

	_log += " ";
	_log += msg;
	_log += "\n";
}

void unit_test::log_error( const std::string& msg, bool warning_only  )
{
    log(msg, false, warning_only);
}


} // namespace vmml
