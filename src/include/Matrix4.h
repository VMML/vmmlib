/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Jonas Boesch
* @author Stefan Eilemann
* @author Renato Pajarola
* @author David H. Eberly ( Wild Magic )
* @author Andrew Willmott ( VL )
*
* @license BSD license, check LICENSE
*
* parts of the source code of VMMLib were inspired by David Eberly's 
* Wild Magic and Andrew Willmott's VL.
* 
*/ 


#ifndef _Matrix4_H_
#define _Matrix4_H_

/* 
*   4x4 Matrix Class
*
*
*/ 

#include "Matrix3.h"
#include "Vector3.h"
#include "Vector4.h"

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

// * * * * * * * * * *
//
// - declaration -
//
// * * * * * * * * * *

namespace vmml
{

template< class Real > 
class Matrix4
{
public:
    union
    {
        struct
        {
            Real m00, m01, m02, m03, 
                 m10, m11, m12, m13, 
                 m20, m21, m22, m23, 
                 m30, m31, m32, m33;
        };
        Real m[4][4]; 
        Real ml[16]; // linear
    };
    
    Matrix4();
    Matrix4( const Matrix4& mm );
    Matrix4( Real v00, Real v01, Real v02, Real v03, 
             Real v10, Real v11, Real v12, Real v13,
             Real v20, Real v21, Real v22, Real v23,
             Real v30, Real v31, Real v32, Real v33 );
    Matrix4( const Vector4<Real>& v0, const Vector4<Real>& v1, 
             const Vector4<Real>& v2, const Vector4<Real>& v3, 
             bool column_vectors = false );

    // dangerous, but implemented to allow easy conversion between 
    // Matrix< float > and Matrix< double >
    //the pointer 'values must be a valid 16 component c array of the resp. type
    Matrix4( float* values );
    Matrix4( double* values );
 
    inline Matrix4& operator= ( const Matrix4& mm );
    inline bool operator== ( const Matrix4& mm ) const;
    inline bool operator!= ( const Matrix4& mm ) const;

    void set( const Matrix4& mm );
    // dangerous, but implemented to allow easy conversion between 
    // Matrix< float > and Matrix< double >
    //the pointer 'values must be a valid 16 component c array of the resp. type
    void set( const float* mm );
    void set( const double* mm );

    inline Real* operator[] ( size_t row ) const;
    Vector4< Real > getRow( size_t row ) const;
    Vector4< Real > getColumn( size_t col ) const; 
    
    // vec3
    void setRow( size_t row, Vector3< Real > rowvec );
    void setColumn( size_t col, Vector3< Real > colvec );
    // vec4
    void setRow( size_t row, Vector4< Real > rowvec );
    void setColumn( size_t col, Vector4< Real > colvec );

    // arithmetic operations
    Matrix4 operator+ ( const Matrix4& mm ) const;
    Matrix4 operator- ( const Matrix4& mm ) const;
    Matrix4 operator* ( const Matrix4& mm ) const;
    Matrix4 operator* ( Real scalar ) const; // matrix = matrix * scalar 
    inline Matrix4 operator/ ( Real scalar ) const
        { scalar = 1.0 / scalar; return operator*(scalar); }; 
    
    // vector = matrix * vector
    Vector3< Real > operator* ( const Vector3< Real >& vv ) const;
    Vector4< Real > operator* ( const Vector4< Real >& vv ) const;

    Matrix4& operator+= ( const Matrix4& mm );
    Matrix4& operator-= ( const Matrix4& mm );
    Matrix4& operator*= ( const Matrix4& mm );
    Matrix4& operator*= ( Real scalar ); // matrix = matrix * scalar 
    inline Matrix4& operator/= ( Real scalar )
        { scalar = 1.0 / scalar; return operator*=( scalar ); };

    Matrix4 negate() const;
    Matrix4 operator-() const;

    Matrix4 transpose() const;
    Real determinant() const;
    Matrix4 adjoint() const;
    
    bool inverse( Matrix4& result, Real limit = 0.0000000001 );

    void rotateX( const Real angle );
    void rotateY( const Real angle );
    void rotateZ( const Real angle );
    void scale( const Real scale[3] );
    void setTranslation( const Real x, const Real y, const Real z );

    void tensor( const Vector3< Real >& u, const Vector3< Real >& v );
    void tensor( const Vector4< Real >& u, const Vector4< Real >& v );

    // computes the cofactor/minor of a submatrix n-1xn-1 = 3x3
    // specify the index of the row/column to be ignored
    Real cofactor( const size_t ignoreRow, const size_t ignoreCol ) const;
    // specify the indices of the rows/columns to be used 
    Real cofactor( const size_t row0, const size_t row1, const size_t row2,
                   const size_t col0, const size_t col1, const size_t col2 )
        const;
                   

    friend std::ostream& operator<<( std::ostream& o, const Matrix4& mm )
    {
        o << "Matrix4( " << mm.m[0][0] << ", " << mm.m[0][1] << ", " 
          << mm.m[0][2] << ", " << mm.m[0][3] << std::endl;
        o << "         " << mm.m[1][0] << ", " << mm.m[1][1] << ", " 
          << mm.m[1][2] << ", " << mm.m[1][3] << std::endl;
        o << "         " << mm.m[2][0] << ", " << mm.m[2][1] << ", " 
          << mm.m[2][2] << ", " << mm.m[2][3] << std::endl;
        o << "         " << mm.m[3][0] << ", " << mm.m[3][1] << ", " 
          << mm.m[3][2] << ", " << mm.m[3][3] << " )" << std::endl;
        return o;
    };  
    
    static const Matrix4 IDENTITY;
    static const Matrix4 ZERO;

    typedef Matrix4<float>  Matrix4f;
    typedef Matrix4<double> Matrix4d;
};

// * * * * * * * * * *
//
// - implementation -
//
// * * * * * * * * * *

template< class Real > 
const Matrix4< Real > Matrix4< Real >::IDENTITY( 1, 0, 0, 0, 0, 1, 0, 0,
                                                 0, 0, 1, 0, 0, 0, 0, 1 );

template< class Real > 
const Matrix4< Real > Matrix4< Real >::ZERO( 0, 0, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0 );


template< class Real > 
Matrix4< Real >::Matrix4()
{}

template< class Real > 
Matrix4< Real >::Matrix4( const Matrix4< Real >& mm )
{
    memcpy(m,mm.m, 16 * sizeof( Real ) );
}

template< class Real > 
Matrix4< Real >::Matrix4( Real v00, Real v01, Real v02, Real v03, 
                          Real v10, Real v11, Real v12, Real v13,
                          Real v20, Real v21, Real v22, Real v23,
                          Real v30, Real v31, Real v32, Real v33 )
    : m00( v00 )
    , m01( v01 )
    , m02( v02 )
    , m03( v03 )
    , m10( v10 )
    , m11( v11 )
    , m12( v12 )
    , m13( v13 )
    , m20( v20 )
    , m21( v21 )
    , m22( v22 )
    , m23( v23 )
    , m30( v30 )
    , m31( v31 )
    , m32( v32 )
    , m33( v33 )
{}
         
template< class Real > 
Matrix4< Real >::Matrix4( const Vector4< Real >& v0, const Vector4< Real >& v1, 
                          const Vector4< Real >& v2, const Vector4< Real >& v3,
                          bool column_vectors )
{
    if ( column_vectors )
        for ( size_t i = 0; i < 4; ++i )
        {
            m[i][0] = v0[i];
            m[i][1] = v1[i];
            m[i][2] = v2[i];
            m[i][3] = v3[i];
        }
    else
        for ( size_t i = 0; i < 4; ++i )
        {
            m[0][i] = v0[i];
            m[1][i] = v1[i];
            m[2][i] = v2[i];
            m[3][i] = v3[i];
        }
}

template< class Real > 
Matrix4< Real >::Matrix4( float* values )
{
    assert( values && "Matrix4: Initialisation of a Matrix from a Nullpointer was requested." );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< Real > ( values[i] );
}

template< class Real > 
Matrix4< Real >::Matrix4( double* values )
{
    assert( values && "Matrix4: Initialisation of a Matrix from a Nullpointer was requested." );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< Real > ( values[i] );
}


template< class Real > 
Matrix4< Real >& Matrix4< Real >::operator= ( const Matrix4< Real >& mm )
{
    memcpy( ml, mm.ml, 16 * sizeof(Real) );
    return *this;
}

template< class Real > 
bool Matrix4< Real >::operator== (const Matrix4< Real >& mm) const
{
    bool equal = true;
    for ( size_t i = 0; i < 16 && equal; ++i )
    {
        equal = ( ml[i] == mm.ml[i] );
    }
    return equal;
}

template< class Real > 
inline bool Matrix4< Real >::operator!= (const Matrix4< Real >& mm) const
{
    return !operator==(mm);
}

template< class Real > 
void Matrix4< Real >::set( const Matrix4& mm )
{
    memcpy( ml, mm.ml, 16 * sizeof( Real ) );
}

template< class Real > 
void Matrix4< Real >::set( const float* values )
{
    assert( values && "Matrix4: Nullpointer argument as source for initialisation!" );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< Real > ( values[i] );
}

template< class Real > 
void Matrix4< Real >::set( const double* values )
{
    assert( values && "Matrix4: Nullpointer argument as source for initialisation!" );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< Real > ( values[i] );
}


// returns a row vector
template< class Real > 
Real* Matrix4< Real >::operator[] (size_t row) const
{
    assert( row < 4 && "Matrix4: Requested Row ( operator[] ) with invalid index!" );
    return const_cast< Real* > ( m[row] );
}

template< class Real > 
Vector4< Real > Matrix4< Real >::getRow( size_t row ) const
{
    assert( row < 4 && "Matrix4: Requested Row ( getRow ) with invalid index!" );
    return Vector4< Real >( m[row] );
}

template< class Real > 
Vector4< Real > Matrix4< Real >::getColumn( size_t col ) const
{
    assert( col < 4 && "Matrix4: Requested Column ( getColumn ) with invalid index!" );
    return Vector4< Real > ( m[0+col], m[4+col], m[8+col], m[12+col] );
}

template< class Real > 
void Matrix4< Real >::setRow( size_t row, Vector4< Real > rowvec )
{
    m[row][0] = rowvec[0];
    m[row][1] = rowvec[1];
    m[row][2] = rowvec[2];
    m[row][3] = rowvec[3];
}

template< class Real > 
void Matrix4< Real >::setColumn( size_t col, Vector4< Real > colvec )
{
    m[0][col] = colvec[0];
    m[1][col] = colvec[1];
    m[2][col] = colvec[2];
    m[3][col] = colvec[3];
}

template< class Real > 
void Matrix4< Real >::setRow( size_t row, Vector3< Real > rowvec )
{
    m[row][0] = rowvec[0];
    m[row][1] = rowvec[1];
    m[row][2] = rowvec[2];
}

template< class Real > 
void Matrix4< Real >::setColumn( size_t col, Vector3< Real > colvec )
{
    m[0][col] = colvec[0];
    m[1][col] = colvec[1];
    m[2][col] = colvec[2];
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::operator+ (const Matrix4< Real >& mm) const
{
    Matrix4< Real > result;
    for ( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] + mm.ml[i];
    return result;
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::operator- (const Matrix4< Real >& mm) const
{
    Matrix4< Real > result;
    for ( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] - mm.ml[i];
    return result;
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::operator* (const Matrix4< Real >& mm) const
{
    Matrix4< Real > result;
    size_t i, j, k;
    Real tmp;

    for (j = 0; j < 4; j++)
        for (i = 0; i < 4; i++) 
        {
            tmp = 0.0;
            for (k = 0; k < 4; k++)
                tmp += m[j][k] * mm.m[k][i];
            result[j][i] = tmp;
        }
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::operator* ( Real scalar ) const
{
    Matrix4< Real > result;
    for ( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] * scalar;
    return result;
}

template< class Real > 
Matrix4< Real >& Matrix4< Real >::operator+= (const Matrix4< Real >& mm) 
{
    for ( size_t i = 0; i < 16; ++i )
        ml[i]  += mm.ml[i];
    return *this;
}

template< class Real > 
Matrix4< Real >& Matrix4< Real >::operator-= ( const Matrix4& mm )
{
    for ( size_t i = 0; i < 16; ++i )
        ml[i]  -= mm.ml[i];
    return *this;
}

template< class Real > 
Matrix4< Real >& Matrix4< Real >::operator*= ( const Matrix4& mm ) 
{
    for( int j = 0; j < 4; ++j )
        for( int i = 0; i < 4; ++i )
        {
            Real tmp = 0.0;
            for( int k = 0; k < 4; ++k )
                tmp += m.m[j][k] * mm.m[k][i];
            m[j][i] = tmp;
        }
    return *this;
}

template< class Real > 
Matrix4< Real >& Matrix4< Real >::operator*= ( Real scalar )
{
    // matrix = matrix * scalar 
    for ( size_t i = 0; i < 16; ++i )
        ml[i]  *= scalar;
    return *this;
}

template< class Real > 
Vector4< Real > Matrix4< Real >::operator* (const Vector4< Real >& vv) const
{
	Vector4< Real > result;
	result[0] = vv[0] * m[0][0] + vv[1] * m[1][0] + vv[2] * m[2][0] + vv[3] * m[3][0];
	result[1] = vv[0] * m[0][1] + vv[1] * m[1][1] + vv[2] * m[2][1] + vv[3] * m[3][1];
	result[2] = vv[0] * m[0][2] + vv[1] * m[1][2] + vv[2] * m[2][2] + vv[3] * m[3][2];
	result[3] = vv[0] * m[0][3] + vv[1] * m[1][3] + vv[2] * m[2][3] + vv[3] * m[3][3];
	return result;
}

template< class Real > 
Vector3< Real > Matrix4< Real >::operator* (const Vector3< Real >& vv) const
{
	Vector3< Real > result;
	result[0] = vv[0] * m[0][0] + vv[1] * m[1][0] + vv[2] * m[2][0];
	result[1] = vv[0] * m[0][1] + vv[1] * m[1][1] + vv[2] * m[2][1];
	result[2] = vv[0] * m[0][2] + vv[1] * m[1][2] + vv[2] * m[2][2];
	return result;
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::transpose() const
{
    Matrix4< Real > result;
    for ( size_t i = 0; i < 4; ++i )
        for ( size_t j = 0; j < 4; ++j )
            result.m[i][j] = m[j][i]; 
    return result;
}

template< class Real > 
Real Matrix4< Real >::determinant() const
{
    return m[0][0] * cofactor( 1, 2, 3, 1, 2, 3 ) 
         - m[0][1] * cofactor( 1, 2, 3, 0, 2, 3 ) 
         + m[0][2] * cofactor( 1, 2, 3, 0, 1, 3 ) 
         - m[0][3] * cofactor( 1, 2, 3, 0, 1, 2 ); 
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::adjoint() const
{
    return Matrix4( 
             cofactor( 1, 2, 3, 1, 2, 3 ),
            -cofactor( 0, 2, 3, 1, 2, 3 ),
             cofactor( 0, 1, 3, 1, 2, 3 ),
            -cofactor( 0, 1, 2, 1, 2, 3 ),
            -cofactor( 1, 2, 3, 0, 2, 3 ),
             cofactor( 0, 2, 3, 0, 2, 3 ),
            -cofactor( 0, 1, 3, 0, 2, 3 ), 
             cofactor( 0, 1, 2, 0, 2, 3 ), 
             cofactor( 1, 2, 3, 0, 1, 3 ),
            -cofactor( 0, 2, 3, 0, 1, 3 ),
             cofactor( 0, 1, 3, 0, 1, 3 ),
            -cofactor( 0, 1, 2, 0, 1, 3 ), 
            -cofactor( 1, 2, 3, 0, 1, 2 ),
             cofactor( 0, 2, 3, 0, 1, 2 ),
            -cofactor( 0, 1, 3, 0, 1, 2 ),
             cofactor( 0, 1, 2, 0, 1, 2 ) );
}

template< class Real > 
bool Matrix4< Real >::inverse( Matrix4< Real >& result, Real limit )
{
	Matrix4< Real > adj( adjoint() );

    Vector4< Real > adv0 ( adj.m[0] );
    Vector4< Real > v0 ( m[0] );
    
	Real det = adv0.dot( v0 );
	
    if ( fabs(det) <= limit )
        return false;
        
	result = adj.transpose();
	result *= ( 1.0/det );

	return true;		
}

template< typename Real >
void Matrix4<Real>::rotateX( const Real angle )
{
    //matrix multiplication: ml = ml * rotation x axis
    const Real sinus = sin(angle);
    const Real cosin = cos(angle);

    Real temp = m[0][0];
    m[0][0] = m[0][0] * cosin - m[0][2] * sinus;
    m[0][2] = temp    * sinus + m[0][2] * cosin;

    temp = m[1][0];
    m[1][0] = m[1][0] * cosin - m[1][2] * sinus;
    m[1][2] = temp    * sinus + m[1][2] * cosin;

    temp = m[2][0];
    m[2][0] = m[2][0] * cosin - m[2][2] * sinus;
    m[2][2] = temp    * sinus + m[2][2] * cosin;

    temp = m[3][0];
    m[3][0] = m[3][0] * cosin - m[3][2] * sinus;
    m[3][2] = temp    * sinus + m[3][2] * cosin;
}
#if 0
template<>
inline void Matrix4<float>::rotateX( const float angle )
{
    //matrix multiplication: ml = ml * rotation x axis
    const float sinus = sinf(angle);
    const float cosin = cosf(angle);

    float temp = m[0][0];
    m[0][0] = m[0][0] * cosin - m[0][2] * sinus;
    m[0][2] = temp    * sinus + m[0][2] * cosin;

    temp = m[1][0];
    m[1][0] = m[1][0] * cosin - m[1][2] * sinus;
    m[1][2] = temp    * sinus + m[1][2] * cosin;

    temp = m[2][0];
    m[2][0] = m[2][0] * cosin - m[2][2] * sinus;
    m[2][2] = temp    * sinus + m[2][2] * cosin;

    temp = m[3][0];
    m[3][0] = m[3][0] * cosin - m[3][2] * sinus;
    m[3][2] = temp    * sinus + m[3][2] * cosin;
}

#endif

template< typename Real >
void Matrix4<Real>::rotateY( const Real angle )
{
    //matrix multiplication: ml = ml * rotation y axis
    const Real sinus = sin(angle);
    const Real cosin = cos(angle);

    Real temp = m[0][1];
    m[0][1] = m[0][1] *  cosin + m[0][2] * sinus;
    m[0][2] = temp    * -sinus + m[0][2] * cosin;

    temp = m[1][1];
    m[1][1] = m[1][1] *  cosin + m[1][2] * sinus;
    m[1][2] = temp    * -sinus + m[1][2] * cosin;

    temp = m[2][1];
    m[2][1] = m[2][1] *  cosin + m[2][2] * sinus;
    m[2][2] = temp    * -sinus + m[2][2] * cosin;

    temp = m[3][1];
    m[3][1] = m[3][1] *  cosin + m[3][2] * sinus;
    m[3][2] = temp    * -sinus + m[3][2] * cosin;
}

#if 0
template<>
inline void Matrix4<float>::rotateY( const float angle )
{
    //matrix multiplication: ml = ml * rotation y axis
    const float sinus = sinf(angle);
    const float cosin = cosf(angle);

    float temp = m[0][1];
    m[0][1] = m[0][1] *  cosin + m[0][2] * sinus;
    m[0][2] = temp    * -sinus + m[0][2] * cosin;

    temp = m[1][1];
    m[1][1] = m[1][1] *  cosin + m[1][2] * sinus;
    m[1][2] = temp    * -sinus + m[1][2] * cosin;

    temp = m[2][1];
    m[2][1] = m[2][1] *  cosin + m[2][2] * sinus;
    m[2][2] = temp    * -sinus + m[2][2] * cosin;

    temp = m[3][1];
    m[3][1] = m[3][1] *  cosin + m[3][2] * sinus;
    m[3][2] = temp    * -sinus + m[3][2] * cosin;
}
#endif

template< typename Real >
void Matrix4<Real>::rotateZ( const Real angle )
{
    //matrix multiplication: ml = ml * rotation z axis
    const Real sinus = sin(angle);
    const Real cosin = cos(angle);

    Real temp = m[0][0];
    m[0][0] = m[0][0] *  cosin + m[0][1] * sinus;
    m[0][1] = temp    * -sinus + m[0][1] * cosin;

    temp = m[1][0];
    m[1][0] = m[1][0] *  cosin + m[1][1] * sinus;
    m[1][1] = temp    * -sinus + m[1][1] * cosin;

    temp = m[2][0];
    m[2][0] = m[2][0] *  cosin + m[2][1] * sinus;
    m[2][1] = temp    * -sinus + m[2][1] * cosin;

    temp = m[3][0];
    m[3][0] = m[3][0] *  cosin + m[3][1] * sinus;
    m[3][1] = temp    * -sinus + m[3][1] * cosin;
}

#if 0
template<>
inline void Matrix4<float>::rotateZ( const float angle )
{
    //matrix multiplication: ml = ml * rotation z axis
    const float sinus = sinf(angle);
    const float cosin = cosf(angle);

    float temp = m[0][0];
    m[0][0] = m[0][0] *  cosin + m[0][1] * sinus;
    m[0][1] = temp    * -sinus + m[0][1] * cosin;

    temp = m[1][0];
    m[1][0] = m[1][0] *  cosin + m[1][1] * sinus;
    m[1][1] = temp    * -sinus + m[1][1] * cosin;

    temp = m[2][0];
    m[2][0] = m[2][0] *  cosin + m[2][1] * sinus;
    m[2][1] = temp    * -sinus + m[2][1] * cosin;

    temp = m[3][0];
    m[3][0] = m[3][0] *  cosin + m[3][1] * sinus;
    m[3][1] = temp    * -sinus + m[3][1] * cosin;
}
#endif

template< typename Real >
void Matrix4<Real>::scale( const Real scale[3] )
{
    ml[0]  *= scale[0];
    ml[4]  *= scale[0];
    ml[8]  *= scale[0];
    ml[12] *= scale[0];
    ml[1]  *= scale[1];
    ml[5]  *= scale[1];
    ml[9]  *= scale[1];
    ml[13] *= scale[1];
    ml[2]  *= scale[2];
    ml[6]  *= scale[2];
    ml[10] *= scale[2];
    ml[14] *= scale[2];
}

template< typename Real >
void Matrix4<Real>::setTranslation( const Real x, const Real y, const Real z )
{
    ml[3] = x;
    ml[7] = y;
    ml[11] = z;
}

template< class Real > 
void Matrix4< Real >::tensor( const Vector3< Real >& u, 
                              const Vector3< Real >& v )
{
    int i, j;
    for (j = 0; j < 3; j++) 
    {
        for (i = 0; i < 3; i++)
            m[j][i] = u[j] * v[i];
        m[j][3] = u[j];
    }
    for (i = 0; i < 3; i++)
        m[3][i] = v[i];
    m[3][3] = 1.0;    
}

template< class Real > 
void Matrix4< Real >::tensor( const Vector4< Real >& u,
                              const Vector4< Real >& v )
{
    int i, j;
    for (j = 0; j < 4; j++)
        for (i = 0; i < 4; i++)
            m[j][i] = u[j] * v[i];

}

template< class Real > 
Real Matrix4< Real >::cofactor( const size_t ignoreRow, 
                                const size_t ignoreCol ) const
{
    Matrix3< Real > cof;
    size_t srcrow, row, srccol, col;
    srcrow = row = srccol = col = 0;
    
    if ( ignoreRow == 0 )
        srcrow++;
    
    for ( ; row < 3; ++row )
    {
        srccol = ( ignoreCol == 0 ) ? 1 : 0;
        for ( col = 0; col < 3; ++col )
        {
            std::cout << "row,col: " << row << "," << col
                      << " srcrow, srccol: " << srcrow << ", " << srccol 
                      << std::endl; 
            cof[row][col] = m[srcrow][srccol++];
            if ( srccol == ignoreCol ) 
                srccol++;
        }
        srcrow++;
        if ( srcrow == ignoreRow )
            srcrow++;
    }
    return cof.determinant();
}

template< class Real > 
Real Matrix4< Real >::cofactor( const size_t row0, const size_t row1,
                                const size_t row2, const size_t col0,
                                const size_t col1, const size_t col2 ) const
{
    Matrix3< Real > cof( m[row0][col0], m[row0][col1], m[row0][col2], 
                         m[row1][col0], m[row1][col1], m[row1][col2],
                         m[row2][col0], m[row2][col1], m[row2][col2] );

    return cof.determinant();
}


template< class Real > 
Matrix4< Real > Matrix4< Real >::operator-() const
{
    Matrix4< Real > result( *this );
    result *= -1.0;
    return result;
}

template< class Real > 
Matrix4< Real > Matrix4< Real >::negate() const
{
    Matrix4< Real > result( *this );
    result *= -1.0;
    return result;
}
}
#endif
