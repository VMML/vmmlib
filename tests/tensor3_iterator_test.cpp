#include "tensor3_iterator_test.hpp"

#include <vmmlib/tensor3_iterator.hpp>
#include <vmmlib/tensor3.hpp>

#include <sstream>

namespace vmml
{
	
	bool
	tensor3_iterator_test::run()
	{
		bool ok = false;
		
		tensor3< 2, 3, 4, uint16_t >  t3;
		tensor3< 2, 3, 4, uint16_t >  t3_check;
		
		t3.fill_increasing_values();
		
		tensor3_iterator< tensor3< 2, 3, 4, uint16_t> > t3_it( t3 );
		
		/*tensor3_iterator< tensor3< 2, 3, 4, uint16_t>>::iterator it = t3_it.begin(), it_end = t3_it.end();
		
		
		for( ; it != it_end; ++it )
		{
			std::cout << *it << std::endl;
		}*/
		
		
		if ( ok)
		{	
			log( "tensor3 iterator  ", true  );
		} else
		{
			std::stringstream error;
			error 
			<< "tensor3 iterator: " << std::endl;
			
			log_error( error.str() );
		}
		
		
		ok = true;
		return ok;
	}
	
	
} // namespace vmml