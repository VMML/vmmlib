/***************************************************************************
*   Copyright (C) 2004 by jonas boesch									  *
*   boesch@ifi.unizh.ch													*
*                                                                         *
***************************************************************************/

#include "VMMLib.h"

#include "JacobiSolver.h"
#include "SingularValueDecomposition.h"

#include <iostream>

typedef vmml::Vector3< float >  Vector3f;
typedef vmml::Vector4< float >  Vector4f;
typedef vmml::Vector3< double > Vector3d;
typedef vmml::Vector4< double > Vector4d;

typedef vmml::Matrix3< float >  Matrix3f;
typedef vmml::Matrix4< float >  Matrix4f;
typedef vmml::Matrix3< double > Matrix3d;
typedef vmml::Matrix4< double > Matrix4d;

int main()
{
    Vector3f a;
    a = 0;
    a.x = 1.0;
    a[2] = 0.3;

    std::cout << a << std::endl << " sizeof " << sizeof(a) << std::endl;

    Vector3f aCopy( a );
    a.normalise();
    Vector3f::normalise( (float*) &aCopy );
    std::cout << a << " " << aCopy << std::endl;
    assert ( a == aCopy );
    
    
    Vector4d b;
    b = 0.6;
    b[3] = 0.1;

    std::cout << b << std::endl << " sizeof " << sizeof(b) << std::endl;
    
    Vector4f c( a, 0 );
    std::cout << c << std::endl << " sizeof " << sizeof(c) << std::endl;

    a.scale( 2 );
    std::cout << a << std::endl << " sizeof " << sizeof(a) << std::endl;

    Vector3f aa0( 0.2 );
    aa0.normalise();
    Vector3f aa1( 0.1, 0.2, 0.3 );
    aa1.normalise();
    Vector3f aa2( 0.9, 0.5, 0.1 );
    aa2.normalise();
    
    a = aa0.normal( aa1, aa2 ); 
    std::cout << a << std::endl << " sizeof " << sizeof(a) << std::endl;
    
    Matrix3d m3_0( vmml::Matrix3< double >::IDENTITY );

    
    
    Matrix4d m4_0( vmml::Matrix4< double >::IDENTITY );

    m3_0.m[2][1] = 3.0;
    m4_0[2][0] = 3.0;

    std::cout << m3_0 << m4_0 << std::endl;

    m3_0 = m3_0.transpose();
    m4_0 = m4_0.transpose();
  
    std::cout << m3_0 << m4_0 << std::endl;
    
    m3_0 = - m3_0;
    Matrix4d m4_1;
    m4_0[2][3] = 5.0;
    m4_0[3][0] = -1.0;
    m4_0[1][0] = 3.0;
    std::cout << "det m4_0: " << m4_0.determinant() << std::endl;  
    if ( m4_0.inverse( m4_1 ) )
        std::cout << "inv ok" << std::endl;
    else
        std::cout << "inv not ok" << std::endl;      
    std::cout << m3_0 << m4_0 << m4_1 << std::endl;

    double aa,bb;
    aa = m4_0.cofactor( 0, 0 );
    bb = m4_0.cofactor( 1, 2, 3, 1, 2, 3 );
    std::cout << "cof test: " << aa << ", " << bb << std::endl;
    
    Matrix3d aaa( aa0.x, aa0.y, aa0.z, aa1.x, aa1.y, aa1.z, aa2.x, aa2.y, aa2.z );
    Matrix3d v = Matrix3d::IDENTITY;
    Vector3d d = 0.0;
    size_t rotations;
    vmml::solveJacobi3x3< double > ( aaa, d, v, rotations );
    
}