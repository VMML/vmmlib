
#include "cp3_tensor_test.hpp"

#include <vmmlib/cp3_tensor.hpp>

#include <sstream>

namespace vmml
{
	
	bool
	cp3_tensor_test::run()
	{
		bool ok = false;
		
		
		
		if ( ok)
		{	
			log( "cp3 tensor test ", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "cp3 tensor test: " << std::endl;
			
			log_error( error.str() );
		}
		
		
		ok = true;
		return ok;
	}
	
	
} // namespace vmml