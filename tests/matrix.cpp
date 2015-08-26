/*
 * Copyright (c) 2006-2015, Visualization and Multimedia Lab,
 *                          University of Zurich <http://vmml.ifi.uzh.ch>,
 *                          Eyescale Software GmbH,
 *                          Blue Brain Project, EPFL
 *
 * This file is part of VMMLib <https://github.com/VMML/vmmlib/>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution.  Neither the name of the Visualization and Multimedia
 * Lab, University of Zurich nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific prior
 * written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <vmmlib/vector.hpp>
#include <vmmlib/matrix.hpp>
#include <vmmlib/math.hpp>

#define BOOST_TEST_MODULE vector
#include <boost/test/unit_test.hpp>

using namespace vmml;

double test_tolerance = 1e-12;

BOOST_AUTO_TEST_CASE(matrix_plus)
{
    Matrix<3, 3, double> A;
    Matrix<3, 3, double> B;
    Matrix<3, 3, double> C;
    Matrix<3, 3, double> groundTruth;

    double A_data[] =
    {
        0.4856,   0.1044,   0.8898,
        0.9691,   0.0663,   0.0371,
        0.2318,   0.5179,   0.0420
    };
    A = A_data;

    double B_data[] =
    {
        0.8542,   0.9748,   0.9285,
        0.7433,   0.7517,   0.7792,
        0.2684,   0.0834,   0.1592
    };
    B = B_data;

    C = A+B;

    double groundTruth_data[] =
    {
        1.339800000000000,  1.079200000000000,  1.818300000000000,
        1.712400000000000,  0.818000000000000,  0.816300000000000,
        0.500200000000000,  0.601300000000000,  0.201200000000000
    };
    groundTruth = groundTruth_data;

    BOOST_CHECK(C.equals(groundTruth,test_tolerance));
}

BOOST_AUTO_TEST_CASE(matrix_product)
{
    Matrix<3, 4, double> A;
    Matrix<4, 3, double> B;
    Matrix<3, 3, double> C;
    Matrix<3, 3, double> groundTruth;

    double A_data[] =
    {
        0.4856,   0.1044,   0.8898,   0.4373,
        0.9691,   0.0663,   0.0371,   0.5836,
        0.2318,   0.5179,   0.0420,   0.5793
    };
    A = A_data;

    double B_data[] =
    {
        0.8542,   0.9748,   0.9285,
        0.7433,   0.7517,   0.7792,
        0.2684,   0.0834,   0.1592,
        0.7066,   0.7254,   0.3720
    };
    B = B_data;

    C = A*B;

    double groundTruth_data[] =
    {
        1.040218540000000,  0.943267100000000,  0.836559840000000,
        1.299415410000000,  1.420953970000000,  1.174475830000000,
        1.003564810000000,  1.038991090000000,  0.840959980000000
    };
    groundTruth = groundTruth_data;

    BOOST_CHECK(C.equals(groundTruth,test_tolerance));
}


BOOST_AUTO_TEST_CASE(matrix_determinant)
{
    Matrix< 2, 2, double > m22;
    Matrix< 3, 3, double > m33;
    Matrix< 4, 4, double > m44;

    double m_data[] =
    {
        1,   2,   0,   0,
        2,   1,   2,   0,
        0,   2,   1,   2,
        0,   0,   2,   1
    };

    m44 = m_data;
    m44.get_sub_matrix( m33 );
    m33.get_sub_matrix( m22 );

    double det44 = m44.det();
    double det33 = m33.det();
    double det22 = m22.det();

    BOOST_CHECK(det44 == 5.0 && det33 == -7.0 && det22 == -3.0);
}

BOOST_AUTO_TEST_CASE(matrix_convolution)
{
    Matrix< 4, 4, int > m1;
    int test_data[] = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2 };
    Matrix< 3, 3, int > kernel;
    int test_kernel[] = {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9 };

    m1 = test_data;
    kernel = test_kernel;
    m1.convolve(kernel);

    Matrix< 4, 4, int > correct_result;
    int corr_res_data[] = {
        159, 192, 237, 264,
        297, 296, 293, 290,
        273, 250, 217, 196,
        231, 198, 153, 126 };
    correct_result = corr_res_data;
    BOOST_CHECK(m1 == correct_result);
}

BOOST_AUTO_TEST_CASE(matrix_inverse)
{
    Matrix< 3, 3 > M, M_inverse, M_inverse_correct;
    double Mdata[] = { 8, 1, 6, 3, 5, 7, 4, 9, 2 };
    M.set( Mdata, Mdata + 9 );

    double M_inverse_correct_data[] =
        {   .14722222222222222222222222222222, -.14444444444444444444444444444444, .63888888888888888888888888888889e-1,
            -.61111111111111111111111111111111e-1, .22222222222222222222222222222222e-1, .10555555555555555555555555555556,
            -.19444444444444444444444444444444e-1, .18888888888888888888888888888889, -.10277777777777777777777777777778 };

    M_inverse_correct.set( M_inverse_correct_data, M_inverse_correct_data + 9 );
    M.inverse( M_inverse );

    BOOST_CHECK(M_inverse.equals(M_inverse_correct,double(test_tolerance)));
}

BOOST_AUTO_TEST_CASE(matrix_frobenius_norm)
{
    Matrix< 4, 4, int > data_2;
    int test_data[] = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2 };
    data_2 = test_data;
    float_t norm_check = 22.0907220343745223090082;
    float_t norm = data_2.frobenius_norm();
    BOOST_CHECK(abs(norm - norm_check) < test_tolerance);
}
