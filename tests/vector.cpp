
/* Copyright (c) 2014, Stefan.Eilemann@epfl.ch
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - Neither the name of Eyescale Software GmbH nor the names of its
 *   contributors may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <vmmlib/vector.hpp>
#include <vmmlib/math.hpp>

#define BOOST_TEST_MODULE vector
#include <boost/test/unit_test.hpp>

using namespace vmml;

BOOST_AUTO_TEST_CASE(vector_base)
{
    Vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    v.iter_set( data, data+4 );

    // tests copyFrom1DimCArray function
    size_t tmp = 1;
    for( size_t index = 0; index < 4; ++index, ++tmp )
    {
        BOOST_CHECK(v.at( index ) == tmp);
    }

    tmp = 4;
    float dataf[] = { 4, 3, 2, 1 };
    v.iter_set( dataf, dataf + 4 );
    for( size_t index = 0; index < 4; ++index, --tmp )
    {
        BOOST_CHECK(v.at( index ) == tmp);
    }
}

BOOST_AUTO_TEST_CASE(vector_plus)
{
    Vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    // tests operator+ function
    Vector< 4, double > v_other;
    Vector< 4, double > v_result;

    v = data;

    double datad[] = { 4, 3, 2, 1 };
    v_other = datad;

    v_result = v + v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == 5);
    }

    v_result = v;
    v_result += v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == 5);
    }

    v = data;
    v_result = v + 2.;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == index + 3);
    }

    v_result = v;
    v_result += 2;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == index + 3);
    }
}

BOOST_AUTO_TEST_CASE(vector_minus)
{
    Vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    // tests operator- function
    Vector< 4, double > v_other;
    Vector< 4, double > v_result;
    v = data;

    double datad[] = { 1, 2, 3, 4 };
    v_other = datad;

    v_result = v - v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == 0);
    }

    v_result = v;
    v_result -= v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == 0);
    }


    v_result = v - 1.0;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == index);
    }

    v_result = v;
    v_result -= 1.0;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == index);
    }
}

BOOST_AUTO_TEST_CASE(vector_times)
{
    Vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    // tests operator* function
    Vector< 4, double > v_other;
    Vector< 4, double > v_result;

    v = data;

    double datad[] = { 24, 12, 8, 6 };
    v_other = datad;

    v_result = v * v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == 24);
    }

    v_result = v;
    v_result *= v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == 24);
    }

    v_result = v * 2.0;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == v.at( index ) * 2.0);
    }

    v_result = v;
    v_result *= 2.0;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(v_result.at( index ) == v.at( index ) * 2.0);
    }
}

BOOST_AUTO_TEST_CASE(vector_div)
{
    Vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    // tests operator/ function
    Vector< 4, double > v_other;
    Vector< 4, double > v_result;

    v = data;

    double datad[] = { 2, 4, 6, 8 };
    v_other = datad;

    v_result = v / v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(( v_result.at( index ) - 0.5 ) < 1e-12);
    }

    v_result = v;
    v_result /= v_other;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(( v_result.at( index ) - 0.5 ) < 1e-12);
    }


    v_result = v / 1.5;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(( v_result.at( index ) - ( v.at( index ) / 1.5 ) ) < 1e-12);
    }

    v_result = v;
    v_result /= 1.5;
    for( size_t index = 0; index < 4; ++index )
    {
        BOOST_CHECK(( v_result.at( index ) - ( v.at( index ) / 1.5 ) ) < 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(vector_norm)
{
    Vector< 4, double > v;
    double data[] = { 1, 2, 3, 4 };

    // tests norm / normSquared (length/lengthSquared) computation
    Vector< 4, double > vec;
    vec = data;

    double normSquared = vec.squared_length();
    BOOST_CHECK(normSquared == 1 * 1 + 2 * 2 + 3 * 3 + 4 * 4);

    double norm = vec.length();
    BOOST_CHECK(sqrt( normSquared ) == norm);

    // tests normalize
    vec = data;
    vec.normalize();
    BOOST_CHECK_CLOSE( vec.length(), 1.0, 0.0000001 );


    // constructor tests
    double vData[] = { 1, 2, 3, 4 };
    Vector< 4, double > v4( 1, 2, 3, 4 );

    Vector< 2, double > v2C;
    v2C = vData;
    Vector< 2, double > v2( 1, 2 );

    BOOST_CHECK(v2 == v2C );

    Vector< 3, double > v3C;
    v3C = vData;
    Vector< 3, double > v3( 1, 2, 3 );

    BOOST_CHECK(v3 == v3C );

    Vector< 4, double > v4C;
    v4C = vData;

    BOOST_CHECK(v4 == v4C);

    double vData2[] = { 23, 23, 23, 23 };
    v4C = vData2;

    Vector< 4, double > v4_( 23 );
    BOOST_CHECK(v4_ == v4C);

    v3 = vData;
    v4C = vData;
    Vector< 4, double > v4from3_1( v3, vData[ 3 ] );
    BOOST_CHECK(v4from3_1 == v4C);

    double hvData[] = { 1., 2., 3., 0.25 };
    double xvData[] = { 4.0, 8.0, 12.0 };

    Vector< 4, double > homogenous;
    homogenous.iter_set( hvData, hvData + 4 );
    Vector< 3, double > nonh;
    nonh.iter_set( xvData, xvData + 3 );

    Vector< 4, double > htest( nonh );

    // to-homogenous-coordinates ctor
    BOOST_CHECK((htest == Vector< 4, double >( 4, 8., 12., 1. ) ));

    Vector< 3, double > nhtest( homogenous );

    // from homogenous-coordiates ctor
    BOOST_CHECK(nhtest == nonh );

    // set tests
    vec.set( 2, 3, 4, 5 );
    Vector< 4, double > vecCorrect;
    double vCData[] = { 2, 3, 4, 5 };
    vecCorrect = vCData;
    BOOST_CHECK(vec == vecCorrect);

    vec.set( 2 );

    double vCData2[] = { 2, 2, 2, 2 };
    vecCorrect = vCData2;
    BOOST_CHECK( vec == vecCorrect );

    Vector< 3, double > v1( 2, 3, 4 );

    // component accessors
    Vector< 4, double > vd( 1, 2, 3, 4 );
    BOOST_CHECK( vd.x() == 1 && vd.y() == 2 && vd.z() == 3 && vd.w() == 4 );
}

BOOST_AUTO_TEST_CASE(vector_dot)
{
    // dot product
    Vector< 3, float > v0( 1, 2, 3 );
    Vector< 3, float > v1( -6, 5, -4 );
    BOOST_CHECK( v0.dot( v1 ) == -8 );
}

BOOST_AUTO_TEST_CASE(vector_cross)
{
    // cross product
    Vector< 3, float > v0( 1, 2, 3 );
    Vector< 3, float > v1( -6, 5, -4 );
    Vector< 3, float > vcorrect( -23, -14, 17 );
    BOOST_CHECK(v0.cross( v1 ) == vcorrect);

    // ???
    Vector< 4, float > vf( -1.0f, 3.0f, -99.0f, -0.9f );
    Vector< 4, size_t > vui( 0, 5, 2, 4 );

    size_t index = vf.find_min_index();
    float f = vf.find_min();

    BOOST_CHECK( index == 2 && f == -99.0f );

    index = vf.find_max_index();
    f = vf.find_max();
    BOOST_CHECK( index == 1 && f == 3.0f );

    index = vui.find_min_index();
    size_t ui = vui.find_min();
    BOOST_CHECK( index == 0 && ui == 0 );

    index = vui.find_max_index();
    ui = vui.find_max();
    BOOST_CHECK( index == 1 && ui == 5 );
}

BOOST_AUTO_TEST_CASE(vector_product)
{
    const Vector< 3, float > v0( 1, 2, 3 );
    BOOST_CHECK_EQUAL( v0.product(), 6 );

    const Vector< 3, int > v1( -6, 5, -4 );
    BOOST_CHECK_EQUAL( v1.product(), 120 );
}

BOOST_AUTO_TEST_CASE(vector_tbd1)
{
    Vector< 4, float > v1( -1.0f, 3.0f, -99.0f, -0.9f );
    float f = 4.0f;
    Vector< 4, float > v_scaled = f * v1;

    BOOST_CHECK(v_scaled == (Vector< 4, float >( -4.0f, 12.0f, -396.0f, -3.6f ) ));


    // ???
    Vector< 3, float > vf( 3.0, 2.0, 1.0 );
    Vector< 3, double > vd( vf );
    Vector< 3, double >::const_iterator it = vd.begin(), it_end = vd.end();
    Vector< 3, float >::const_iterator fit = vf.begin();
    for( ; it != it_end; ++it, ++fit )
    {
        BOOST_CHECK(*it == *fit);
    }
    vd = vf;
    for( ; it != it_end; ++it, ++fit )
    {
        BOOST_CHECK(*it == *fit);
    }
}

BOOST_AUTO_TEST_CASE(vector_tbd2)
{
    // ???
    Vector< 4, float > vf( 3.0, 2.0, 1.0, 1.0 );
    Vector< 3, float >& v3 = vf.get_sub_vector< 3 >();
    BOOST_CHECK(v3.x() == vf.x() && v3.y() == vf.y());
    v3.normalize();

    BOOST_CHECK(v3.x() == vf.x() && v3.y() == vf.y());

    //elementwise sqrt
    Vector< 4, double > vsq(9.0, 4.0, 1.0, 2.0);
    Vector< 4, double > vsq_check( 3.0, 2.0, 1.0, std::sqrt( 2.0 ));
    vsq.sqrt_elementwise();
    BOOST_CHECK_EQUAL( vsq, vsq_check );

    //elementwise sqrt
    Vector< 4, float > vr( 9.0, 4.0, 1.0, 2.0 );
    Vector< 4, float > vr_check( 0.1111111119389534, 0.25, 1, 0.5 );
    vr.reciprocal();
    BOOST_CHECK(vr == vr_check);
}

BOOST_AUTO_TEST_CASE(vector_l2norm)
{
    Vector< 4, float > vr( 9.0, 4.0, 1.0, 2.0 );
    double v_norm_check = 10.09950493836208;
    double v_norm = vr.norm();

    BOOST_CHECK((v_norm - v_norm_check) < 0.0001);
}
