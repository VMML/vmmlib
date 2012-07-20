#include "tensor4_test.hpp"

#include <vmmlib/tensor4.hpp>
#include <sstream>

namespace vmml
{
	
	bool
	tensor4_test::run()
	{
		bool ok = false;
		
		tensor4< 2, 3, 4, 2, int >  t4;
		//t4.zero();
		
		if ( ok )
		{	
			log( "tensor4 method x", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "tensor4 method x. t4 should be: " 
			<< std::endl 
			//<< t4
			<< std::endl;
			log_error( error.str() );
		}
		
		
		return ok;
	}
	
	
	
} // namespace vmml

