#include "t3_padder_test.hpp"

#include <vmmlib/t3_padder.hpp>
#include <vmmlib/t3_converter.hpp>
#include <sstream>

namespace vmml
{
	
	bool t3_padder_test::run()
	{
        bool ok = false;
		
        // define Tensor
		const size_t i = 16;
        const size_t i1 = 3;
        const size_t i2 = 4;
        const size_t i3 = 4;
        const size_t b = 2;
        typedef int T;
		
        typedef tensor3< i1, i2, i3, T >  t3_t;
        typedef tensor3< b, b, b, T >  t3_block_t;
		typedef t3_converter< i1, i2, i3, T >  t3_conv_t;
		typedef t3_padder< i, i1, i2, i3, T >  t3_padder_t;
		
		std::string dir = "./";
		std::string filename = "tmp_sus.raw";
		std::string path = dir; path.append( filename );
		
		
		//create tmp file
		t3_t t3; 
        t3.fill_increasing_values();

		// current working directory '.' is the same in different OS
        remove( path.c_str() ); // remove just for security
		
		t3_conv_t conv;
        conv.write_to_raw( t3, dir, filename );
		
		//actual padder access test
		t3_block_t block; block.zero();
		t3_padder_t t3_p( dir, filename );
		
		t3_p.get_tensor3_block( block, 0, 0, 0 );
		std::cout << "full data block\n" << block << std::endl;
		
		t3_p.get_tensor3_block( block, 4, 4, 4 );
		std::cout << "empty data block\n" << block << std::endl;

		t3_p.get_tensor3_block( block, 2, 2, 2 );
		std::cout << "half full data block\n" << block << std::endl;

 		
        // remove file
        int done = remove( filename.c_str() );
        if (done != 0)
        {
            std::cerr << "*** Warning ***" << std::endl << "The temporary file, which was created, could not be deleted, it is in the working Directory under the name .~tmpRaw.raw and can safely be deleted manually." << std::endl << std::endl;
        }
		
        // Test writing and reading raw
        if ( false )
		{
            log( "t3 padder test", true);
			
		} else
		{
			
			std::stringstream error;
			error << "t3 padder test: "
            << std::endl
			<< t3
			<< std::endl;
			log_error( error.str(),  false);
		}
		
		
		
		return ok;
	}
	
	
	
} // namespace vmml

