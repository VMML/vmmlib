#include "t3_ihooi_test.hpp"
#include "vmmlib/t3_ihooi.hpp"
#include <sstream>

namespace vmml {

    bool
    t3_ihooi_test::run() {
        bool ok = false;

        double precision = 0.001;

#define I 4
#define R 4
#define NBLOCKS 2

        tensor3< I, I, I, double > t3_input;
        double data_input[] = {
            0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858,
            0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
            0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
            0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071
        };
        t3_input.set(data_input, data_input + I * I * I);

        typedef matrix< I, R, double > t3_u_type;
        typedef tensor3< R, R, R, double> t3_core_type;

        t3_core_type t3_core;
        t3_u_type u1;
        t3_u_type u2;
        t3_u_type u3;

        typedef t3_hooi< R, R, R, I, I, I > hooi_type;
        t3_ihooi< R, R, R, NBLOCKS, I, I, I, double, double >::i_als(t3_input, u1, u2, u3, t3_core, hooi_type::init_hosvd(), 20, -1);

        //std::cout << "u1:\n" << u1 << std::endl << "u2:\n" << u2 << std::endl << "u3:\n" << u3 << std::endl << "lambda\n" << lambda << std::endl;

        //check test data
        t3_core_type t3_core_check;
        t3_u_type u1_check;
        t3_u_type u2_check;
        t3_u_type u3_check;

        double data_check[] = {
            3.0755, -0.1154, 0, 0,
            -0.0723, 0.2529, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0.0132, 0.6130, 0, 0,
            0.9410, 0.3878, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0.5961, 0.1045,
            0, 0, -0.0986, 0.0287,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, -0.0955, 0.1952,
            0, 0, -0.3350, 0.1218,
        };
        t3_core_check.set(data_check, data_check + I * I * I);

        double data_u1_check[] = {
            0.4409, -0.4897, -0.2280, -0.3190,
            0.3960, -0.4795, -0.2584, -0.8490,
            0.3230, -0.2798, -0.1235, -0.2404,
            0.7379, 0.6723, 0.9306, -0.3459,
        };
        u1_check.set(data_u1_check, data_u1_check + I * R);

        double data_u2_check[] = {
            -0.4196, 0.6487, 0.8321, 0.3869,
            -0.5298, -0.2388, 0.2436, 0.3261,
            -0.6362, -0.5215, -0.2148, 0.0691,
            -0.3720, 0.5003, -0.4496, 0.8598,
        };
        u2_check.set(data_u2_check, data_u2_check + I * R);

        double data_u3_check[] = {
            -0.6458, -0.7248, -0.1541, 0.2925,
            -0.4492, 0.1096, 0.9545, 0.1567,
            -0.4903, 0.5689, 0.1873, -0.8779,
            -0.3753, 0.3727, -0.1736, -0.3452,
        };

        u3_check.set(data_u3_check, data_u3_check + I * R);

        ok = t3_core.equals(t3_core_check, precision) && ok;
        ok = u1.equals(u1_check, precision);
        ok = u2.equals(u2_check, precision) && ok;
        ok = u3.equals(u3_check, precision) && ok;

        if (ok) {
            log("incremental Tucker approximation", ok);
        } else {
            std::stringstream error;
            error
                    << "incremental Tucker approximation" << std::setprecision(16) << std::endl
                    << "core should be: " << t3_core_check << std::endl
                    << "core is: " << t3_core << std::endl
                    << "u1 should be: " << std::endl << u1_check << std::endl
                    << "u1 is: " << std::endl << u1 << std::endl
                    << "u2 should be: " << std::endl << u2_check << std::endl
                    << "u2 is: " << std::endl << u2 << std::endl
                    << "u3 should be: " << std::endl << u3_check << std::endl
                    << "u3 is: " << std::endl << u3 << std::endl;

            log_error(error.str());
        }

        return ok;
    }

} //end vmml namespace
