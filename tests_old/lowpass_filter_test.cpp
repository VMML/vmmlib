#include "lowpass_filter_test.hpp"

#include <vmmlib/lowpass_filter.hpp>
#include <vmmlib/vector.hpp>
#include <deque>
#include <sstream>
#include <cmath>

namespace vmml
{
bool lowpass_filter_test::run()
{
    bool ok = true;
    bool global_ok = true;
    double data[] = { 0, 2, 4, 8, 16 };
    Vector< 5, double > v;
    v.iter_set(data, data+5);

    LowpassFilter< 4, double > filter (.5f);
    filter.add(data[0]);
    filter.add(data[1]);
    filter.add(data[2]);
    filter.add(data[3]);
    filter.add(data[4]);

    float tmp = 9;
    const double filtered = filter.get();

    TEST(filtered == tmp);

    log( "low pass filter, filter ( data, smooth_factor )", ok );
    if ( ! ok )
    {
        std::stringstream error;
        error << "Filter " << v << "\n"
              << "result should be: " << tmp << ",\n"
              << "result is: " << filtered << std::endl;
        log_error( error.str() );
    }

    return global_ok;
}

} // namespace vmml
