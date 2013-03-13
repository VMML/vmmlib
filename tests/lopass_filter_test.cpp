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
    vector< 5, double > v;
    v.iter_set(data, data+5);

    lopass_filter< 4, double > filter (.1f);
    filter.add(data[0]);
    filter.add(data[1]);
    filter.add(data[2]);
    filter.add(data[3]);
    filter.add(data[4]);

    filter.set_smooth_factor(.5f);

    float tmp = 9;
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
