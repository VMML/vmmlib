#include "util_test.hpp"

#include <vmmlib/vector.hpp>
#include <vmmlib/matrix.hpp>
#include <vmmlib/util.hpp>
#include <sstream>
#include <cmath>

#define TEST( x ) \
{ \
    ok = (x); \
    global_ok &= ok; \
}

namespace vmml
{


bool
util_test::run()
{
    bool global_ok = true;
    bool ok = true;
    
    // tests scalar equals
    ok = true;
    {
        double d1 = 0.12345;
        double d2 = 0.23456;
        double d3 = 0.35801;
        
        float  f1 = 0.12345f;
        float  f2 = 0.23456f;
        float  f3 = 0.35801f;
        
        TEST( equals( d1 + d2, d3 ) )
        
        log( "scalar equals (double)", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << d1 << " + " << d2 << " != " << d3 << std::endl;
            log_error( error.str() );
        }
        
        TEST( equals( f1 + f2, f3 ) )
        
        log( "scalar equals (float)", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << f1 << " + " << f2 << " != " << f3 << std::endl;
            log_error( error.str() );
        }
        
        TEST( equals( f1 + d2, d3 ) )
        
        log( "scalar equals (mixed)", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << f1 << " + " << d2 << " != " << f3 << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests create_translation
    ok = true;
    {
        vector< 3, double > v(1, 2, 3);
        matrix< 4, 4, double > tmp = matrix< 4, 4, double >::IDENTITY;
        tmp.set_column(3, vector< 4, double >(v, 1));

        matrix<4, 4, double> m = create_translation(v);
        
        TEST( equals( m, tmp ) )
        
        log( "create_translation", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << m << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests apply_translation
    ok = true;
    {
        vector< 3, double > v(1.0, 2.0, 3.0);
        matrix< 4, 4, double > tmp = matrix< 4, 4, double >::IDENTITY;
        tmp.set_column(3, vector< 4, double >(v, 1));

        matrix<4, 4, double> m = matrix<4, 4, double>::IDENTITY;
        apply_translation(m, 1.0, 2.0, 3.0);
        
        TEST( equals( m, tmp ) )
        
        log( "apply_translation", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << m << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests create_rotation
    ok = true;
    {
        vec4f v( -0.5f, 0.5f, 0.0f, 1.0f );
        v = create_rotation( M_PI_F, vec3f( 1.0f, 1.0f, 0.0f) ) * v;
        vec4f tmp( 0.5f, -0.5f, 0.0, 1.0f );
        
        TEST( equals( v, tmp ) )
        
        log( "create_rotation", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << v << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests apply_rotation
    ok = true;
    {
        vec4f v( -0.5f, 0.0f, 0.5f, 1.0f );
        mat4f m = mat4f::IDENTITY;
        apply_rotation( m, M_PI_F, 1.0f, 1.0f, 1.0f );
        v = m * v;
        vec4f tmp( 0.5f, 0.0f, -0.5f, 1.0f );
        
        TEST( equals( v, tmp ) )
        
        log( "apply_rotation", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << v << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests create_scaling
    ok = true;
    {
        vec4f v(1.0f, 2.0f, 3.0f, 1.0f);
        v = create_scaling(2.0f) * v;
        vec4f tmp(2.0f, 4.0f, 6.0f, 1.0f);
        
        TEST( equals( v, tmp ) )
        
        log( "create_scaling", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << v << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests manhattan
    ok = true;
    {
        vec3f v1(1.0, 2.0, 3.0);
        vec3f v2(2.0, 1.0, 4.0);
        float d = manhattan(v2 - v1);
        float tmp = 3;

        TEST( d == tmp )
        
        log( "manhattan", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << d << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests zero
    ok = true;
    {
        vec4f v(1.0, 2.0, 3.0, 1.0);
        zero( v );
        vec4f tmp = vec4f::ZERO;

        TEST( equals( v, tmp ) )
        
        log( "zero", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << v << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests min
    ok = true;
    {
        vec3f v1(1.0, 2.0, 3.0);
        vec3f v2(2.0, 1.0, 4.0);
        vec3f v = min(v1, v2);
        vec3f tmp(1.0f, 1.0f, 3.0f);

        TEST( equals( v, tmp ) )
        
        log( "min", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << v << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    // tests max
    ok = true;
    {
        vec3f v1(1.0, 2.0, 3.0);
        vec3f v2(2.0, 1.0, 4.0);
        vec3f v = max(v1, v2);
        vec3f tmp(2.0f, 2.0f, 4.0f);

        TEST( equals( v, tmp ) )
        
        log( "max", ok  );
        if ( !ok )
        {
            std::stringstream error;
            error << v << " != " << tmp << std::endl;
            log_error( error.str() );
        }
    }
    
    return global_ok;
}

} // namespace vmml

