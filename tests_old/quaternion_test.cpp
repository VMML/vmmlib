#include "quaternion_test.hpp"

#include <vmmlib/math.hpp>
#include <vmmlib/quaternion.hpp>

#include <iostream>

using namespace std;
namespace vmml
{
typedef Quaternion< float > Quaternionf;
typedef Quaternion< double > Quaterniond;


bool quaternion_test::run()
{
    bool global_ok = true;
    Quaternion< double > q;
    double QData[] = { 1., 6., 3., 8.  };
    for( size_t index = 0; index < 4; ++index )
    {
        q.array[ index ] = QData[ index ];
    }

    // operator==/!= tests
    {
        bool ok = true;
        Quaterniond qq, qqq;
        qq.array[ 0 ] = qqq.array[ 0 ] = 1.0;
        qq.array[ 1 ] = qqq.array[ 1 ] = 6.0;
        qq.array[ 2 ] = qqq.array[ 2 ] = 3.0;
        qq.array[ 3 ] = qqq.array[ 3 ] = 8.0;
        TEST(qq == qqq);
        log( "operator==, operator!=", ok );
    }


    // operator= tests
    {
        bool ok = true;
        Quaterniond tquaternion_test = q;
        TEST(tquaternion_test == q );

        tquaternion_test.iter_set< double* >( QData, QData + 4 );
        if ( ok )
            TEST(tquaternion_test == q);



        log( "operator=", ok );
    }

    // ctor tests
    {
        bool ok = true;
        Quaterniond qq( q );
        TEST(q == qq);

        Quaterniond t( 1., 6., 3., 8 );
        if ( ok )
            TEST(q == t);

        Vector< 3, double > xyz;
        double xyzData[] = { 1., 6., 3. };
        xyz = xyzData;

        Quaterniond s( xyz, 8 );
        if ( ok )
            TEST(q == s);

        Matrix< 3, 3, double > mat;
        double matData[] = { 1., 0., 0., 0., 0., 1., 0., -1., 0. };
        mat = matData;
        const Quaterniond u( mat );

        TEST( u.w() == sqrt( 2. ) / 2. &&
              u.x() == - ( 1 / sqrt( 2. ) ) &&
              (u.y() == 0 && u.z() == 0 ) );
        log( "constructors", ok );

    }

    // set test
    {
        bool ok = true;
        Quaterniond qqq;
        qqq.set ( 1., 6., 3., 8. );
        TEST( qqq == q );
        log( "set( x,y,z,w )", ok );
    }

    // abs
    {
        bool ok = true;
        TEST(q.abs() == sqrt( 110.0 ) && q.squared_abs() == 110);
        log( "abs(), squared_abs()", ok );
    }


    //conjugate() test
    {
        bool ok = true;
        Quaterniond conj(  -1., -6., -3., 8. );
        TEST( q.get_conjugate() == conj );

        conj.conjugate();
        if ( ok )
            TEST(q == conj);

        log( "conjugate()", ok );
    }

    // quat / scalar operations
    {
        bool ok = true;
        Quaterniond t;
        t.set( 1, 2, 3, 4 );

        Quaterniond t3;
        t3.set( 3, 6, 9, 12 );

        double f = 3.0;
        double rf = 1./ f;

        t *= f;
        TEST(t == t3);

        t.set( 1, 2, 3, 4 );
        t /= rf;
        if (ok)
            TEST(t == t3);

        t.set( 1, 2, 3, 4 );
        if (ok)
            TEST(t * f == t3);

        t.set( 1, 2, 3, 4 );
        if (ok)
            TEST(t / rf == t3);

        log( "quaternion / scalar operations: operator*, /, *=, /=", ok );
    }

    {
        bool ok = true;

        Quaterniond qq;
        qq.set( 8, 3, 6, 1 );

        Quaterniond qpqq;
        qpqq.set(  9., 9., 9., 9.   );

        // +, +=
        TEST(q + qq == qpqq );

        qq += q;
        if (ok)
            TEST(qq == qpqq);

        // -, -=
        qq.set( 8, 3, 6, 1 );
        if (ok)
            TEST(qpqq - qq == q);

        qpqq -= qq;
        if (ok)
            TEST(qpqq == q);

        // *, *=
        qq.set( 2, 3, 4, 1 );
        Quaterniond q2( 3, 2, 1, 4 );
        Quaterniond p = qq * q2;
        Quaterniond pCorrect( 6, 24, 12, -12 );
        if (ok)
            TEST(p == pCorrect);

        p = qq;
        p *= q2;
        if (ok)
            TEST(p == pCorrect);

        log( "quaternion / quaternion operations: operator+, -, *, +=, -=, *=", ok );


    }

    {
        bool ok = true;

        Quaterniond qq( 1, 2, 3, 4 );
        Quaterniond q2( -6, 5, -4, 2 );

        Vector< 3, double > v = qq.cross( q2 );

        Vector< 3, double > v0( 1, 2, 3 );
        Vector< 3, double > v1( -6, 5, -4 );

        TEST( v == v0.cross( v1 ) );

        log( "cross product ( vec3 = quat x quat ).", ok );
    }

    // TODO
//    {
//        // FIXME: find correct values -> actually test if you get the
//        // correct values
//        bool ok = true;
//
//        quaterniond qq(1., 2., 3., 4.);
//        quaterniond pp(1., 2., 3., 4.);
//        double a = 0.;
//
//        quaterniond x(1., 2., 3., 4.);
//
//        quaterniond::slerp(a, qq, pp);
//
//        quaterniond correct(0.0, 1.0, 1.0, 0.0);
//
//        TEST(correct == x);
//
//        log( "FIXME: todo slerp(a, p, q ).", ok );
//    }


#if 0
	Quaterniond rquaternion_test( 8. / 110., -1. / 110., -6. / 110., -3. / 110. );

	if ( _quaternion.invert() != rquaternion_test )
	{
		cout << "test: Quaternion::invert() failed!" << endl;
        failed();
        assert( 0 );
	}

	if ( tquaternion_test.dot( _quaternion ) != 40. )
	{
		cout << "test: Quaternion::dot( quaternion ) failed!" << endl;
        failed();
        assert( 0 );
	}

    Vector3d vector3_test( -30., -4., 18. );
    if ( tquaternion_test.cross( _quaternion ) != vector3_test )
	{
		cout << "test: Quaternion::cross( quaternion ) failed!" << endl;
        failed();
        assert( 0 );
	}

	Vector3d yaxis( 1., 0., 0. );
	Quaterniond svector3_test( 0., 18., -30., -4. );
	Quaterniond result_test( _quaternion.rotate( M_PI / 2.f, yaxis, vector3_test ) );
	if ( abs( result_test.abs() - svector3_test.abs() ) > 1e-13 )
	{
		cout << "test: Quaternion::rotate( T, Vector3, Vector3 ) failed!" << endl;
		failed();
		assert( 0 );
	}

    if ( ok )
        cout << "Quaternion: all tests passed!" << endl;

#endif
	return global_ok;
}

}
