#ifndef VMMLIB__UNIT_TEST_GLOBALS__HPP
#define VMMLIB__UNIT_TEST_GLOBALS__HPP

#include <string>
#include <iostream>

namespace vmml
{

class unit_test_globals
{
public:
    // singleton getter
    static unit_test_globals& get_instance();
    
    static const size_t     OK      = 0;
    static const size_t     WARN    = 1;
    static const size_t     FAIL    = 2;
    
    void notify_test( size_t count = 1 );
    void notify_warn( size_t count = 1 );
    void notify_fail( size_t count = 1 );
    
    size_t get_test_count() const;
    size_t get_warn_count() const;
    size_t get_fail_count() const;

    void set_ok_prefix(     const std::string& prefix_ );
    void set_warn_prefix(   const std::string& prefix_ );
    void set_fail_prefix(   const std::string& prefix_ );
    
    const std::string& get_fail_prefix() const;
    const std::string& get_warn_prefix() const;
    const std::string& get_ok_prefix() const;

    static void destroy_instance(); 
    
    friend std::ostream& operator<<( std::ostream& os,
        const unit_test_globals& g ) 
    {
        os  << "====== TEST OVERVIEW ====== \n"
            << "tests:        "     << g.get_test_count() << "\n"
            << "warnings:     "     << g.get_warn_count() << "\n"
            << "failed tests: " << g.get_fail_count() << "\n"
            << "=========================== \n";
        return os;
    }
       
protected:
    static unit_test_globals* _instance;

    std::string _fail_prefix;
    std::string _warn_prefix;
    std::string _ok_prefix;

    size_t      _test_count;
    size_t      _warn_count;
    size_t      _fail_count;

private:
    unit_test_globals();


}; // class unit_test_globals

} // namespace vmml

#endif

