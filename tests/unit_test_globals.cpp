#include "unit_test_globals.hpp"

namespace vmml
{

unit_test_globals* unit_test_globals::_instance = 0;

unit_test_globals::unit_test_globals()
    : _fail_prefix( "[ FAIL ]" )
    , _warn_prefix( "[ WARN ]" )
    , _ok_prefix(   "[  ok  ]" )
    , _test_count( 0 )
    , _warn_count( 0 )
    , _fail_count( 0 )
{}

void
unit_test_globals::notify_test( size_t count )
{
    _test_count += count;
}


void
unit_test_globals::notify_warn( size_t count )
{
    _warn_count += count;
}


void
unit_test_globals::notify_fail( size_t count )
{
    _fail_count += count;
}



void
unit_test_globals::set_ok_prefix( const std::string& prefix_ )
{
    _ok_prefix = prefix_;
}


void
unit_test_globals::set_warn_prefix( const std::string& prefix_ )
{
    _warn_prefix = prefix_;
}


void
unit_test_globals::set_fail_prefix( const std::string& prefix_ )
{
    _fail_prefix = prefix_;
}


const std::string&
unit_test_globals::get_fail_prefix() const
{
    return _fail_prefix;
}


const std::string&
unit_test_globals::get_warn_prefix() const
{
    return _warn_prefix;
}


const std::string&
unit_test_globals::get_ok_prefix() const
{
    return _ok_prefix;
}


unit_test_globals&
unit_test_globals::get_instance()
{
    if ( ! _instance )
        _instance = new unit_test_globals();
    return *_instance;
}


size_t
unit_test_globals::get_test_count() const
{
    return _test_count;
}


size_t
unit_test_globals::get_warn_count() const
{
    return _warn_count;
}


size_t
unit_test_globals::get_fail_count() const
{
    return _fail_count;
}


void
unit_test_globals::destroy_instance()
{
    if ( _instance )
        delete _instance;
}

} // namespace vmml
