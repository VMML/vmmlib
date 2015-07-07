#include "util_test.hpp"

#include <vmmlib/vector.hpp>
#include <vmmlib/matrix.hpp>
#include <vmmlib/util.hpp>
#include <sstream>
#include <cmath>

namespace vmml
{
bool util_test::run()
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
        Vector< 3, double > v(1, 2, 3);
        Matrix< 4, 4, double > tmp = Matrix< 4, 4, double >::IDENTITY;
        tmp.set_column(3, Vector< 4, double >(v, 1));

        Matrix<4, 4, double> m = create_translation(v);

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
        Vector< 3, double > v(1.0, 2.0, 3.0);
        Matrix< 4, 4, double > tmp = Matrix< 4, 4, double >::IDENTITY;
        tmp.set_column(3, Vector< 4, double >(v, 1));

        Matrix<4, 4, double> m = Matrix<4, 4, double>::IDENTITY;
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
        Vector< 4, float > v( -0.5f, 0.5f, 0.0f, 1.0f );
        v = create_rotation( M_PI_F, Vector< 3, float >( 1.0f, 1.0f, 0.0f) ) * v;
        Vector< 4, float > tmp( 0.5f, -0.5f, 0.0, 1.0f );

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
        Vector< 4, float > v( -0.5f, 0.0f, 0.5f, 1.0f );
        Matrix4f m = Matrix4f::IDENTITY;
        apply_rotation( m, M_PI_F, 1.0f, 1.0f, 1.0f );
        v = m * v;
        Vector< 4, float > tmp( 0.5f, 0.0f, -0.5f, 1.0f );

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
        Vector< 4, float > v(1.0f, 2.0f, 3.0f, 1.0f);
        v = create_scaling(2.0f) * v;
        Vector< 4, float > tmp(2.0f, 4.0f, 6.0f, 1.0f);

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
        Vector< 3, float > v1(1.0, 2.0, 3.0);
        Vector< 3, float > v2(2.0, 1.0, 4.0);
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
        Vector< 4, float > v(1.0, 2.0, 3.0, 1.0);
        zero( v );
        Vector< 4, float > tmp = Vector< 4, float >::ZERO;

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
        Vector< 3, float > v1(1.0, 2.0, 3.0);
        Vector< 3, float > v2(2.0, 1.0, 4.0);
        Vector< 3, float > v = min(v1, v2);
        Vector< 3, float > tmp(1.0f, 1.0f, 3.0f);

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
        Vector< 3, float > v1(1.0, 2.0, 3.0);
        Vector< 3, float > v2(2.0, 1.0, 4.0);
        Vector< 3, float > v = max(v1, v2);
        Vector< 3, float > tmp(2.0f, 2.0f, 4.0f);

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
