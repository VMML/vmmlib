
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

#include <vmmlib/aabb.hpp>

#define BOOST_TEST_MODULE axisAlignedBoundingBox
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(axisAlignedBoundingBox_base)
{
    vmml::AABBf box1;
    BOOST_CHECK_EQUAL( box1.isEmpty(), true );

    const vmml::Vector3f p1( 0.f, 0.f, 0.f );
    box1.merge( p1 );
    BOOST_CHECK_EQUAL( box1.isEmpty(), true );
    BOOST_CHECK_EQUAL( box1.getMin(), p1 );
    BOOST_CHECK_EQUAL( box1.getMax(), p1 );

    const vmml::Vector3f p2( 1.f, 1.f, 1.f );
    box1.merge( p2 );
    BOOST_CHECK_EQUAL( box1.isEmpty(), false );
    BOOST_CHECK_EQUAL( box1.getMin(), p1 );
    BOOST_CHECK_EQUAL( box1.getMax(), p2 );
    BOOST_CHECK_EQUAL( box1.getDimension(), p2 );

    const vmml::AABBf box2( -p2, p2 );
    BOOST_CHECK_EQUAL( box2.isEmpty(), false );
    BOOST_CHECK_EQUAL( box2.getMin(), -p2 );
    BOOST_CHECK_EQUAL( box2.getMax(), p2 );
    BOOST_CHECK_EQUAL( box2.getDimension(), vmml::Vector3f( 2.f, 2.f, 2.f ));

    box1.merge( box2 );
    BOOST_CHECK_EQUAL( box1, box2 );
}
