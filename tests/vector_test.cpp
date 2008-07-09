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


    // tests operator+ function
	ok = true;
	{
        vector< 4 > v_other;
        vector< 4 > v_result;
        
        v = data;
        
        double datad[] = { 4, 3, 2, 1 };
        v_other = datad;

        v_result = v + v_other;

		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 5;
		}

        v_result = v;
        v_result += v_other;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 5;
		}

		log( "operator+ / operator+=", ok  );
		if ( ! ok )
		{
			std::stringstream error;
			error 
                << "\n"
                << "v        " << v 
                << "v_other  " << v_other
                << "v_result " << v_result
                << std::endl;
			log_error( error.str() );
		}
	}


    // tests operator- function
	ok = true;
	{
        vector< 4 > v_other;
        vector< 4 > v_result;
        
        v = data;
        
        double datad[] = { 1, 2, 3, 4 };
        v_other = datad;

        v_result = v - v_other;

		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 0;
		}

        v_result = v;
        v_result -= v_other;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 0;
		}

		log( "operator- / operator-=", ok  );
		if ( ! ok )
		{
			std::stringstream error;
			error 
                << "\n"
                << "v        " << v 
                << "v_other  " << v_other
                << "v_result " << v_result
                << std::endl;
			log_error( error.str() );
		}
	}


    // tests operator* function
	ok = true;
	{
        vector< 4 > v_other;
        vector< 4 > v_result;
        
        v = data;
        
        double datad[] = { 24, 12, 8, 6 };
        v_other = datad;

        v_result = v * v_other;

		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 24;
		}

        v_result = v;
        v_result *= v_other;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 24;
		}

		log( "operator* / operator*=", ok  );
		if ( ! ok )
		{
			std::stringstream error;
			error 
                << "\n"
                << "v        " << v 
                << "v_other  " << v_other
                << "v_result " << v_result
                << std::endl;
			log_error( error.str() );
		}
	}


    // tests operator/ function
	ok = true;
	{
        vector< 4 > v_other;
        vector< 4 > v_result;
        
        v = data;
        
        double datad[] = { 2, 4, 6, 8 };
        v_other = datad;

        v_result = v / v_other;

		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 0.5;
		}

        v_result = v;
        v_result /= v_other;
        
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == 0.5;
		}

		log( "operator/ / operator/=", ok  );
		if ( ! ok )
		{
			std::stringstream error;
			error 
                << "\n"
                << "v        " << v 
                << "v_other  " << v_other
                << "v_result " << v_result
                << std::endl;
			log_error( error.str() );
		}
	}


    return ok;
}

} // namespace vmml

