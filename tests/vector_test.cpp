#include "vector_test.hpp"

#include <vmmlib/vector.hpp>
#include <sstream>

namespace vmml
{
bool
vector_test::run()
{
    bool ok = true;
    
    vector< 4 > v;
    double data[] = { 1, 2, 3, 4 };
    
    v.copyFrom1DimCArray( data );
    
    // tests copyFrom1DimCArray function
	ok = true;
	{
		size_t tmp = 1;
		for( size_t index = 0; ok && index < 4; ++index, ++tmp )
		{
            ok = v.at( index ) == tmp;
		}
        
        tmp = 4;
        float dataf[] = { 4, 3, 2, 1 };
        v.copyFrom1DimCArray( dataf );
		for( size_t index = 0; ok && index < 4; ++index, --tmp )
		{
            ok = v.at( index ) == tmp;
		}

		log( "copyFrom1DimCArray( ..  )", ok  );
		if ( ! ok )
		{
			std::stringstream error;
			error << v << std::endl;
			log_error( error.str() );
		}
	
	}



    return ok;
}

} // namespace vmml

