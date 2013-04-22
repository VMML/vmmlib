#include "t4_converter_test.hpp"

#include <vmmlib/tensor4.hpp>
#include <vmmlib/t4_converter.hpp>
#include <sstream>

#define TEST( x ) \
{ \
    ok = x; \
    global_ok &= ok; \
}

namespace vmml
{

	bool t4_converter_test::run()
	{
        bool global_ok = true;
        bool ok = false;
        // indicates if failing the test produces only a warning
        bool fail_test = false;

        // define Tensor
		const size_t a = 3;
        const size_t b = 2;
        const size_t c = 4;
        const size_t d = 4;
        typedef int T;

        tensor4< a, b, c, d, T >  t4;
        tensor4< a, b, c, d, T >  t4r;
        t4r.zero();
        t4_converter<a, b, c, d, T> conv;

        t4.fill_random_signed();

        // current working directory '.' is the same in different OS
        remove(".~tmpRaw.raw"); // remove just for security

        conv.write_to_raw(t4, "./", ".~tmpRaw");
        conv.read_from_raw(t4r, "./", ".~tmpRaw");

        // remove file

        int done = remove(".~tmpRaw.raw");
        if (done != 0)
        {
            std::cerr << "*** Warning ***" << std::endl << "The temporary file, which was created, could not be deleted, it is in the working Directory under the name .~tmpRaw.raw and can safely be deleted manually." << std::endl << std::endl;
        }


        // Test writing and reading raw
        TEST( t4 == t4r );
        if ( ok )
		{
            log( "tensor4 IO write/read .raw file", true);

		} else
		{

			std::stringstream error;
			error << "tensor4 IO write and read for raw files. t4 read from file is: "
            << std::endl
			<< t4r
			<< std::endl
            << "original t4 is:" << std::endl
            << t4 << std::endl;
			log_error( error.str(),  fail_test);
		}



		return global_ok;
	}



} // namespace vmml

