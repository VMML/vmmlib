#include "lopass_filter_test.hpp"

#include <vmmlib/lopass_filter.hpp>
#include <vmmlib/vector.hpp>
#include <deque>
#include <sstream>
#include <cmath>

namespace vmml
{


bool
lopass_filter_test::run()
{
    double data[] = { 0, 2, 4, 8, 16 };

    lopass_filter< 4, double > filter (.5f);
    filter.add_value(data[0]);
    filter.add_value(data[1]);
    filter.add_value(data[2]);
    filter.add_value(data[3]);
    filter.add_value(data[4]);

    float tmp = 9;
    double filtered = filter.filter();

    bool ok = filtered == tmp;

    log( "low pass filter, filter ( data, smooth_factor )", ok );
    if ( ! ok )
    {
        std::stringstream error;
        error << "Filter " << data << "\n"
              << "result should be: " << tmp << ",\n"
              << "result is: " << filtered << std::endl;
        log_error( error.str() );
    }

    return ok;
}

} // namespace vmml

