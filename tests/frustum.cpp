
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

#include <vmmlib/frustum.hpp>
#include <vmmlib/frustum_culler.hpp>

#define BOOST_TEST_MODULE frustum
#include <boost/test/unit_test.hpp>

static void _testCull( const vmml::FrustumCuller< float >& fc )
{
    const vmml::Vector< 4, float > sphereIn( 0.f, 0.f, -10.f, 1.f );
    const vmml::Vector< 4, float > sphereOut( 0.f, 0.f, 0.f, .5f );
    const vmml::Vector< 4, float > sphereBorder( 0.f, 0.f, -1.f, 1.f );

    BOOST_CHECK_EQUAL( fc.test_sphere( sphereIn ), vmml::VISIBILITY_FULL );
    BOOST_CHECK_EQUAL( fc.test_sphere( sphereOut ), vmml::VISIBILITY_NONE );
    BOOST_CHECK_EQUAL( fc.test_sphere( sphereBorder ), vmml::VISIBILITY_PARTIAL );

    const vmml::Vector< 2, float > xy( -1.f, 1.f );
    const vmml::Vector< 2, float > zIn( -2.f, -4.f );
    const vmml::Vector< 2, float > zOut( 0.f, -.5f );
    const vmml::Vector< 2, float > zBorder( -.5f, -1.5f );

    BOOST_CHECK_EQUAL( fc.test_aabb( xy, xy, zIn ), vmml::VISIBILITY_FULL );
    BOOST_CHECK_EQUAL( fc.test_aabb( xy, xy, zOut ), vmml::VISIBILITY_NONE );
    BOOST_CHECK_EQUAL( fc.test_aabb( xy, xy, zBorder ), vmml::VISIBILITY_PARTIAL );
}


BOOST_AUTO_TEST_CASE(frustum_base)
{
    const vmml::Frustum< float > Frustum( -1.f, 1., -1.f, 1., 1.f, 100.f );
    const vmml::Matrix< 4, 4, float > mvp = Frustum.compute_matrix();

    vmml::FrustumCuller< float > fc;
    fc.setup( mvp );
    _testCull( fc );

    //   e_____f
    //  /     /|
    // | a b | |
    // | c d |/h
    //  -----
    const vmml::Vector< 3, float > a( -1.f,  1.f, -1.f );
    const vmml::Vector< 3, float > b(  1.f,  1.f, -1.f );
    const vmml::Vector< 3, float > c( -1.f, -1.f, -1.f );
    const vmml::Vector< 3, float > d(  1.f, -1.f, -1.f );
    const vmml::Vector< 3, float > e( -100.f,  100.f, -100.f );
    const vmml::Vector< 3, float > f(  100.f,  100.f, -100.f );
    const vmml::Vector< 3, float > g( -100.f, -100.f, -100.f );
    const vmml::Vector< 3, float > h(  100.f, -100.f, -100.f );

    fc.setup( a, b, c, d, e, f, g, h );
    _testCull( fc );
}
