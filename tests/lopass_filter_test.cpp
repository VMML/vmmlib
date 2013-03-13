#include "lopass_filter_test.hpp"

#include <vmmlib/lopass_filter.hpp>
#include <vmmlib/vector.hpp>
#include <sstream>
#include <cmath>

namespace vmml
{


bool
lopass_filter_test::run()
{
    vector< 4, double > v;
    double data[] = { 2, 4, 8, 16 };
    float tmp = 9;

    v.iter_set( data, data+4 );

    lopass_filter< 4, double > filter (v, .5f);
    double filtered = filter.filter();

    bool ok = filtered == tmp;

    log( "low pass filter, filter ( data, smooth_factor )", ok );
    if ( ! ok )
    {
        std::stringstream error;
        error << "Filter " << v << "\n"
              << "result should be: " << tmp << ",\n"
              << "result is: " << filtered << std::endl;
        log_error( error.str() );
    }

    return ok;
}

} // namespace vmml

