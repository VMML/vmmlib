#include "vector_test.hpp"

#include <vmmlib/vector.hpp>
#include <sstream>
#include <cmath>

namespace vmml
{
bool
vector_test::run()
{
    bool ok = true;
    
    vector< 4 > v;
    double data[] = { 1, 2, 3, 4 };
    
    double tolerance = 1e-15;
    
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

        v = data;
        v_result = v + 2;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == index + 3;
		}
        
        v_result = v;
        v_result += 2;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == index + 3;
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


        v_result = v - 1.0;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == index;
		}

        v_result = v;
        v_result -= 1.0;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == index;
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

        v_result = v * 2.0;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == v.at( index ) * 2.0;
		}

        v_result = v;
        v_result *= 2.0;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == v.at( index ) * 2.0;
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


        v_result = v / 1.5;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == v.at( index ) / 1.5;
		}

        v_result = v;
        v_result /= 1.5;
		for( size_t index = 0; ok && index < 4; ++index )
		{
            ok = v_result.at( index ) == v.at( index ) / 1.5;
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

    // tests norm / normSquared (length/lengthSquared) computation
	ok = true;
	{
        vector< 4 > vec;
        vec = data;
        
        double normSquared = vec.normSquared();
        ok = normSquared == 1 * 1 + 2 * 2 + 3 * 3 + 4 * 4;

        double norm = vec.norm();
        if ( ok ) 
            ok = sqrt( normSquared ) == norm;

		log( "norm / normSquared", ok  );

    }


    // tests normalize
	ok = true;
	{
        vector< 4 > vec;
        vec = data;
        vec.normalize();
        ok = vec.norm() == 1.0;

		log( "normalize() maximum precision", ok  );
        if ( ! ok )
        {
            ok = vec.norm() - 1.0 < 1e-15;
            log( "normalize() with tolerance 1e-15", ok  );
        }
        if ( ! ok )
        {
            std::stringstream ss;
            ss << "norm after normalize() " << vec.norm() << std::endl;
            log_error( ss.str() );
        }

    }


    // tests 
	ok = true;
	{
        vector< 4 > vec;

    }


    return ok;
}

} // namespace vmml

