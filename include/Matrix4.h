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
#include <iomanip>
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

template< typename T > 
class Matrix4
{
public:
    union
    {
        struct
        {
            T   m00, m10, m20, m30, 
                m01, m11, m21, m31, 
                m02, m12, m22, m32, 
                m03, m13, m23, m33;
        };
        T m[4][4]; 
        T ml[16]; // linear
    };
    
    Matrix4();
    Matrix4( const Matrix4& mm );
    Matrix4( T v00, T v01, T v02, T v03, 
             T v10, T v11, T v12, T v13,
             T v20, T v21, T v22, T v23,
             T v30, T v31, T v32, T v33 );
    Matrix4( const Vector4<T>& v0, const Vector4<T>& v1, 
             const Vector4<T>& v2, const Vector4<T>& v3, 
             bool columnVectors = false );

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

    Vector4< T > getColumn( size_t column ) const;
    Vector4< T > getRow( size_t row ) const; 
    
    // vec3
    void setColumn( size_t column, const Vector3< T >& columnvec );
    void setRow( size_t row, const Vector3< T >& rowvec );
    // vec4
    void setColumn( size_t column, const Vector4< T >& columnvec );
    void setRow( size_t row, const Vector4< T >& rowvec );

    // arithmetic operations
    Matrix4 operator+ ( const Matrix4& mm ) const;
    Matrix4 operator- ( const Matrix4& mm ) const;
    Matrix4 operator* ( const Matrix4& mm ) const;
    Matrix4 operator* ( T scalar ) const; // matrix = matrix * scalar 
    inline Matrix4 operator/ ( T scalar ) const
        { scalar = 1.0 / scalar; return operator*(scalar); }; 
    
    // vector = matrix * vector
    Vector3< T > operator* ( const Vector3< T >& vv ) const;
    Vector4< T > operator* ( const Vector4< T >& vv ) const;

    Matrix4& operator+= ( const Matrix4& mm );
    Matrix4& operator-= ( const Matrix4& mm );
    Matrix4& operator*= ( const Matrix4& mm );
    Matrix4& operator*= ( T scalar ); // matrix = matrix * scalar 
    inline Matrix4& operator/= ( T scalar )
        { scalar = 1.0 / scalar; return operator*=( scalar ); };

    Matrix4 negate() const;
    Matrix4 operator-() const;

    Matrix4 transpose() const;
    T determinant() const;
    Matrix4 adjoint() const;
    
    bool inverse( Matrix4& result, T limit = 0.0000000001 );

    void rotateX( const T angle );
    void rotateY( const T angle );
    void rotateZ( const T angle );
    void scale( const T scale[3] );
    void scaleTranslation( const T scale[3] );
    void setTranslation( const T x, const T y, const T z );
    void setTranslation( const T trans[3] );

    void tensor( const Vector3< T >& u, const Vector3< T >& v );
    void tensor( const Vector4< T >& u, const Vector4< T >& v );

    // computes the cofactor/minor of a submatrix n-1xn-1 = 3x3
    // specify the index of the row/column to be ignored
    T cofactor( const size_t ignoreRow, const size_t ignoreCol ) const;
    // specify the indices of the rows/columns to be used 
    T cofactor( const size_t row0, const size_t row1, const size_t row2,
                   const size_t col0, const size_t col1, const size_t col2 )
        const;
                   

    friend std::ostream& operator << ( std::ostream& os, const Matrix4& m )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << std::endl << "|" << std::setw(7) << m.m[0][0] << " " 
           << std::setw(7) << m.m[0][1] << " " 
           << std::setw(7) << m.m[0][2] << " " 
           << std::setw(7) << m.m[0][3] << "|" << std::endl
           << "|" << std::setw(7) << m.m[1][0] << " " 
           << std::setw(7) << m.m[1][1] << " " 
           << std::setw(7) << m.m[1][2] << " " 
           << std::setw(7) << m.m[1][3] << "|" << std::endl
           << "|" << std::setw(7) << m.m[2][0] << " " 
           << std::setw(7) << m.m[2][1] << " " 
           << std::setw(7) << m.m[2][2] << " " 
           << std::setw(7) << m.m[2][3] << "|" << std::endl
           << "|" << std::setw(7) << m.m[3][0] << " " 
           << std::setw(7) << m.m[3][1] << " " 
           << std::setw(7) << m.m[3][2] << " " 
           << std::setw(7) << m.m[3][3] << "|" << std::endl;
        os.precision( prec );
        os.setf( flags );
        return os;
    };  
    
    static const Matrix4 IDENTITY;
    static const Matrix4 ZERO;

};

#ifndef VMMLIB_DISABLE_TYPEDEFS
    typedef Matrix4<float>  Matrix4f;
    typedef Matrix4<double> Matrix4d;
#endif

// * * * * * * * * * *
//
// - implementation -
//
// * * * * * * * * * *

template< typename T > 
const Matrix4< T > Matrix4< T >::IDENTITY( 1, 0, 0, 0, 0, 1, 0, 0,
                                           0, 0, 1, 0, 0, 0, 0, 1 );

template< typename T > 
const Matrix4< T > Matrix4< T >::ZERO( 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0 );


template< typename T > 
Matrix4< T >::Matrix4()
{}

template< typename T > 
Matrix4< T >::Matrix4( const Matrix4< T >& mm )
{
    memcpy(m,mm.m, 16 * sizeof( T ) );
}

template< typename T > 
Matrix4< T >::Matrix4( T v00, T v01, T v02, T v03, 
                       T v10, T v11, T v12, T v13,
                       T v20, T v21, T v22, T v23,
                       T v30, T v31, T v32, T v33 )
    : m00( v00 )
    , m10( v10 )
    , m20( v20 )
    , m30( v30 )
    , m01( v01 )
    , m11( v11 )
    , m21( v21 )
    , m31( v31 )
    , m02( v02 )
    , m12( v12 )
    , m22( v22 )
    , m32( v32 )
    , m03( v03 )
    , m13( v13 )
    , m23( v23 )
    , m33( v33 )
{}
         
template< typename T > 
Matrix4< T >::Matrix4( const Vector4< T >& v0, const Vector4< T >& v1, 
                          const Vector4< T >& v2, const Vector4< T >& v3,
                          bool columnVectors )
{
    if ( columnVectors )
        for ( size_t i = 0; i < 4; ++i )
        {
            m[0][i] = v0[i];
            m[1][i] = v1[i];
            m[2][i] = v2[i];
            m[3][i] = v3[i];
        }
    else
        for ( size_t i = 0; i < 4; ++i )
        {
            m[i][0] = v0[i];
            m[i][1] = v1[i];
            m[i][2] = v2[i];
            m[i][3] = v3[i];
        }
}

template< typename T > 
Matrix4< T >::Matrix4( float* values )
{
    assert( values && "Matrix4: Initialisation of a Matrix from a Nullpointer was requested." );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< T > ( values[i] );
}

template< typename T > 
Matrix4< T >::Matrix4( double* values )
{
    assert( values && "Matrix4: Initialisation of a Matrix from a Nullpointer was requested." );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< T > ( values[i] );
}


template< typename T > 
Matrix4< T >& Matrix4< T >::operator= ( const Matrix4< T >& mm )
{
    memcpy( ml, mm.ml, 16 * sizeof(T) );
    return *this;
}

template< typename T > 
bool Matrix4< T >::operator== (const Matrix4< T >& mm) const
{
    for( size_t i = 0; i < 16; ++i )
    {
        if( ml[i] != mm.ml[i] )
            return false;
    }
    return true;
}

template< typename T > 
inline bool Matrix4< T >::operator!= (const Matrix4< T >& mm) const
{
    return !operator==(mm);
}

template< typename T > 
void Matrix4< T >::set( const Matrix4& mm )
{
    memcpy( ml, mm.ml, 16 * sizeof( T ) );
}

template< typename T > 
void Matrix4< T >::set( const float* values )
{
    assert( values && "Matrix4: Nullpointer argument as source for initialisation!" );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< T > ( values[i] );
}

template< typename T > 
void Matrix4< T >::set( const double* values )
{
    assert( values && "Matrix4: Nullpointer argument as source for initialisation!" );
    for ( size_t i = 0; i < 16; ++i )
        ml[i] = static_cast< T > ( values[i] );
}

template< typename T > 
Vector4< T > Matrix4< T >::getColumn( size_t column ) const
{
    assert( column < 4 && "Matrix4: Requested Column ( getColumn ) with invalid index!" );
    return Vector4< T >( m[column] );
}

template< typename T > 
Vector4< T > Matrix4< T >::getRow( size_t row ) const
{
    assert( row < 4 && "Matrix4: Requested Row ( getRow ) with invalid index!");
    return Vector4< T > ( m[0+row], m[4+row], m[8+row], m[12+row] );
}

template< typename T > 
void Matrix4< T >::setColumn( size_t column, const Vector4< T >& columnvec )
{
    m[column][0] = columnvec[0];
    m[column][1] = columnvec[1];
    m[column][2] = columnvec[2];
    m[column][3] = columnvec[3];
}

template< typename T > 
void Matrix4< T >::setRow( size_t row, const Vector4< T >& rowvec )
{
    m[0][row] = rowvec[0];
    m[1][row] = rowvec[1];
    m[2][row] = rowvec[2];
    m[3][row] = rowvec[3];
}

template< typename T > 
void Matrix4< T >::setColumn( size_t column, const Vector3< T >& columnvec )
{
    m[column][0] = columnvec[0];
    m[column][1] = columnvec[1];
    m[column][2] = columnvec[2];
}

template< typename T > 
void Matrix4< T >::setRow( size_t row, const Vector3< T >& rowvec )
{
    m[0][row] = rowvec[0];
    m[1][row] = rowvec[1];
    m[2][row] = rowvec[2];
}

template< typename T > 
Matrix4< T > Matrix4< T >::operator+ (const Matrix4< T >& mm) const
{
    Matrix4< T > result;
    for ( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] + mm.ml[i];
    return result;
}

template< typename T > 
Matrix4< T > Matrix4< T >::operator- (const Matrix4< T >& mm) const
{
    Matrix4< T > result;
    for( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] - mm.ml[i];
    return result;
}

template< typename T > 
Matrix4< T > Matrix4< T >::operator* (const Matrix4< T >& o) const
{
    Matrix4< T > r;

    r.m00 = m00*o.m00 + m10*o.m01 + m20*o.m02 + m30*o.m03;
    r.m10 = m00*o.m10 + m10*o.m11 + m20*o.m12 + m30*o.m13;
    r.m20 = m00*o.m20 + m10*o.m21 + m20*o.m22 + m30*o.m23;
    r.m30 = m00*o.m30 + m10*o.m31 + m20*o.m32 + m30*o.m33;

    r.m01 = m01*o.m00 + m11*o.m01 + m21*o.m02 + m31*o.m03;
    r.m11 = m01*o.m10 + m11*o.m11 + m21*o.m12 + m31*o.m13;
    r.m21 = m01*o.m20 + m11*o.m21 + m21*o.m22 + m31*o.m23;
    r.m31 = m01*o.m30 + m11*o.m31 + m21*o.m32 + m31*o.m33;

    r.m00 = m02*o.m00 + m12*o.m01 + m22*o.m02 + m32*o.m03;
    r.m10 = m02*o.m10 + m12*o.m11 + m22*o.m12 + m32*o.m13;
    r.m20 = m02*o.m20 + m12*o.m21 + m22*o.m22 + m32*o.m23;
    r.m30 = m02*o.m30 + m12*o.m31 + m22*o.m32 + m32*o.m33;

    r.m01 = m03*o.m00 + m13*o.m01 + m23*o.m02 + m33*o.m03;
    r.m11 = m03*o.m10 + m13*o.m11 + m23*o.m12 + m33*o.m13;
    r.m21 = m03*o.m20 + m13*o.m21 + m23*o.m22 + m33*o.m23;
    r.m31 = m03*o.m30 + m13*o.m31 + m23*o.m32 + m33*o.m33;

    return r;
}

template< typename T > 
Matrix4< T > Matrix4< T >::operator* ( T scalar ) const
{
    Matrix4< T > result;
    for ( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] * scalar;
    return result;
}

template< typename T > 
Matrix4< T >& Matrix4< T >::operator+= (const Matrix4< T >& mm) 
{
    for ( size_t i = 0; i < 16; ++i )
        ml[i]  += mm.ml[i];
    return *this;
}

template< typename T > 
Matrix4< T >& Matrix4< T >::operator-= ( const Matrix4& mm )
{
    for ( size_t i = 0; i < 16; ++i )
        ml[i]  -= mm.ml[i];
    return *this;
}

template< typename T > 
Matrix4< T >& Matrix4< T >::operator*= ( const Matrix4& o ) 
{
    T t[4];
    t[0] = m00*o.m00 + m10*o.m01 + m20*o.m02 + m30*o.m03;
    t[1] = m00*o.m10 + m10*o.m11 + m20*o.m12 + m30*o.m13;
    t[2] = m00*o.m20 + m10*o.m21 + m20*o.m22 + m30*o.m23;
    t[3] = m00*o.m30 + m10*o.m31 + m20*o.m32 + m30*o.m33;
    m00 = t[0]; m10 = t[1]; m20 = t[2]; m30 = t[3];

    t[0] = m01*o.m00 + m11*o.m01 + m21*o.m02 + m31*o.m03;
    t[1] = m01*o.m10 + m11*o.m11 + m21*o.m12 + m31*o.m13;
    t[2] = m01*o.m20 + m11*o.m21 + m21*o.m22 + m31*o.m23;
    t[3] = m01*o.m30 + m11*o.m31 + m21*o.m32 + m31*o.m33;
    m01 = t[0]; m11 = t[1]; m21 = t[2]; m31 = t[3];

    t[0] = m02*o.m00 + m12*o.m01 + m22*o.m02 + m32*o.m03;
    t[1] = m02*o.m10 + m12*o.m11 + m22*o.m12 + m32*o.m13;
    t[2] = m02*o.m20 + m12*o.m21 + m22*o.m22 + m32*o.m23;
    t[3] = m02*o.m30 + m12*o.m31 + m22*o.m32 + m32*o.m33;
    m02 = t[0]; m12 = t[1]; m22 = t[2]; m32 = t[3];

    t[0] = m03*o.m00 + m13*o.m01 + m23*o.m02 + m33*o.m03;
    t[1] = m03*o.m10 + m13*o.m11 + m23*o.m12 + m33*o.m13;
    t[2] = m03*o.m20 + m13*o.m21 + m23*o.m22 + m33*o.m23;
    t[3] = m03*o.m30 + m13*o.m31 + m23*o.m32 + m33*o.m33;
    m03 = t[0]; m13 = t[1]; m23 = t[2]; m33 = t[3];

    return *this;
}

template< typename T > 
Matrix4< T >& Matrix4< T >::operator*= ( T scalar )
{
    // matrix = matrix * scalar 
    for ( size_t i = 0; i < 16; ++i )
        ml[i]  *= scalar;
    return *this;
}

template< typename T > 
Vector4< T > Matrix4< T >::operator* (const Vector4< T >& vv) const
{
	return Vector4< T >( vv[0] * m00 + vv[1] * m10 + vv[2] * m20 + vv[3] * m30,
                         vv[0] * m01 + vv[1] * m11 + vv[2] * m21 + vv[3] * m31,
                         vv[0] * m02 + vv[1] * m12 + vv[2] * m22 + vv[3] * m32,
                         vv[0] * m03 + vv[1] * m13 + vv[2] * m23 + vv[3] * m33);
}

template< typename T > 
Vector3< T > Matrix4< T >::operator* (const Vector3< T >& vv) const
{
	Vector4< T > result;
    
	result[0] = vv[0] * m00 + vv[1] * m10 + vv[2] * m20 + m30;
	result[1] = vv[0] * m01 + vv[1] * m11 + vv[2] * m21 + m31;
	result[2] = vv[0] * m02 + vv[1] * m12 + vv[2] * m22 + m32;
	result[3] = vv[0] * m03 + vv[1] * m13 + vv[2] * m23 + m33;
	return Vector3<T>( result.x/result.w, result.y/result.w, 
                       result.z/result.w );
}

template< typename T > 
Matrix4< T > Matrix4< T >::transpose() const
{
    Matrix4< T > result;
    for ( size_t i = 0; i < 4; ++i )
        for ( size_t j = 0; j < 4; ++j )
            result.m[i][j] = m[j][i]; 
    return result;
}

template< typename T > 
T Matrix4< T >::determinant() const
{
    return m00 * cofactor( 1, 2, 3, 1, 2, 3 ) 
         - m01 * cofactor( 1, 2, 3, 0, 2, 3 ) 
         + m02 * cofactor( 1, 2, 3, 0, 1, 3 ) 
         - m03 * cofactor( 1, 2, 3, 0, 1, 2 ); 
}

template< typename T > 
Matrix4< T > Matrix4< T >::adjoint() const
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

template< typename T > 
bool Matrix4< T >::inverse( Matrix4< T >& result, T limit )
{
	Matrix4< T > adj( adjoint() );

    Vector4< T > adv0 ( adj.m[0] );
    Vector4< T > v0 ( m[0] );
    
	T det = adv0.dot( v0 );
	
    if ( fabs(det) <= limit )
        return false;
        
	result = adj.transpose();
	result *= ( 1.0/det );

	return true;		
}

template< typename T >
void Matrix4<T>::rotateX( const T angle )
{
    //matrix multiplication: ml = ml * rotation x axis
    const T sinus = sin(angle);
    const T cosin = cos(angle);

    T temp = m00;
    m00 = m00  * cosin - m02 * sinus;
    m02 = temp * sinus + m02 * cosin;

    temp = m10;
    m10 = m10  * cosin - m12 * sinus;
    m12 = temp * sinus + m12 * cosin;

    temp = m20;
    m20 = m20  * cosin - m22 * sinus;
    m22 = temp * sinus + m22 * cosin;

    temp = m30;
    m30 = m30  * cosin - m32 * sinus;
    m32 = temp * sinus + m32 * cosin;
}

template<>
inline void Matrix4<float>::rotateX( const float angle )
{
    //matrix multiplication: ml = ml * rotation x axis
    const float sinus = sinf(angle);
    const float cosin = cosf(angle);

    float temp = m00;
    m00 = m00  * cosin - m02 * sinus;
    m02 = temp * sinus + m02 * cosin;

    temp = m10;
    m10 = m10  * cosin - m12 * sinus;
    m12 = temp * sinus + m12 * cosin;

    temp = m20;
    m20 = m20  * cosin - m22 * sinus;
    m22 = temp * sinus + m22 * cosin;

    temp = m30;
    m30 = m30  * cosin - m32 * sinus;
    m32 = temp * sinus + m32 * cosin;
}

template< typename T >
void Matrix4<T>::rotateY( const T angle )
{
    //matrix multiplication: ml = ml * rotation y axis
    const T sinus = sin(angle);
    const T cosin = cos(angle);

    T temp = m01;
    m01 = m01  *  cosin + m02 * sinus;
    m02 = temp * -sinus + m02 * cosin;

    temp = m11;
    m11 = m11  *  cosin + m12 * sinus;
    m12 = temp * -sinus + m12 * cosin;

    temp = m21;
    m21 = m21  *  cosin + m22 * sinus;
    m22 = temp * -sinus + m22 * cosin;

    temp = m31;
    m31 = m31  *  cosin + m32 * sinus;
    m32 = temp * -sinus + m32 * cosin;
}

template<>
inline void Matrix4<float>::rotateY( const float angle )
{
    //matrix multiplication: ml = ml * rotation y axis
    const float sinus = sinf(angle);
    const float cosin = cosf(angle);

    float temp = m01;
    m01 = m01  *  cosin + m02 * sinus;
    m02 = temp * -sinus + m02 * cosin;

    temp = m11;
    m11 = m11  *  cosin + m12 * sinus;
    m12 = temp * -sinus + m12 * cosin;

    temp = m21;
    m21 = m21  *  cosin + m22 * sinus;
    m22 = temp * -sinus + m22 * cosin;

    temp = m31;
    m31 = m31  *  cosin + m32 * sinus;
    m32 = temp * -sinus + m32 * cosin;
}

template< typename T >
void Matrix4<T>::rotateZ( const T angle )
{
    //matrix multiplication: ml = ml * rotation z axis
    const T sinus = sin(angle);
    const T cosin = cos(angle);

    T temp = m00;
    m00 = m00  *  cosin + m01 * sinus;
    m01 = temp * -sinus + m01 * cosin;

    temp = m10;
    m10 = m10  *  cosin + m11 * sinus;
    m11 = temp * -sinus + m11 * cosin;

    temp = m20;
    m20 = m20  *  cosin + m21 * sinus;
    m21 = temp * -sinus + m21 * cosin;

    temp = m30;
    m30 = m30  *  cosin + m31 * sinus;
    m31 = temp * -sinus + m31 * cosin;
}

template<>
inline void Matrix4<float>::rotateZ( const float angle )
{
    //matrix multiplication: ml = ml * rotation z axis
    const float sinus = sinf(angle);
    const float cosin = cosf(angle);

    float temp = m00;
    m00 = m00  *  cosin + m01 * sinus;
    m01 = temp * -sinus + m01 * cosin;

    temp = m10;
    m10 = m10  *  cosin + m11 * sinus;
    m11 = temp * -sinus + m11 * cosin;

    temp = m20;
    m20 = m20  *  cosin + m21 * sinus;
    m21 = temp * -sinus + m21 * cosin;

    temp = m30;
    m30 = m30  *  cosin + m31 * sinus;
    m31 = temp * -sinus + m31 * cosin;
}

template< typename T >
void Matrix4<T>::scale( const T scale[3] )
{
    ml[0]  *= scale[0];
    ml[1]  *= scale[0];
    ml[2]  *= scale[0];
    ml[3]  *= scale[0];
    ml[4]  *= scale[1];
    ml[5]  *= scale[1];
    ml[6]  *= scale[1];
    ml[7]  *= scale[1];
    ml[8]  *= scale[2];
    ml[9]  *= scale[2];
    ml[10] *= scale[2];
    ml[11] *= scale[2];
}

template< typename T >
void Matrix4<T>::scaleTranslation( const T scale[3] )
{
    ml[12] *= scale[0];
    ml[13] *= scale[1];
    ml[14] *= scale[2];
}

template< typename T >
void Matrix4<T>::setTranslation( const T x, const T y, const T z )
{
    ml[12] = x;
    ml[13] = y;
    ml[14] = z;
}
template< typename T >
void Matrix4<T>::setTranslation( const T trans[3] )
{
    ml[12]  = trans[0];
    ml[13]  = trans[1];
    ml[14] = trans[2];
}

template< typename T > 
void Matrix4< T >::tensor( const Vector3< T >& u, 
                              const Vector3< T >& v )
{
    int i, j;
    for (j = 0; j < 3; j++) 
    {
        for (i = 0; i < 3; i++)
            m[i][j] = u[j] * v[i];
        m[3][j] = u[j];
    }
    for (i = 0; i < 3; i++)
        m[i][3] = v[i];
    m[3][3] = 1.0;    
}

template< typename T > 
void Matrix4< T >::tensor( const Vector4< T >& u,
                              const Vector4< T >& v )
{
    int i, j;
    for (j = 0; j < 4; j++)
        for (i = 0; i < 4; i++)
            m[i][j] = u[j] * v[i];

}

template< typename T > 
T Matrix4< T >::cofactor( const size_t ignoreRow, const size_t ignoreCol ) const
{
    Matrix3< T > cof;
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
            cof[col][row] = m[srcrow][srccol++];
            if ( srccol == ignoreCol ) 
                srccol++;
        }
        srcrow++;
        if ( srcrow == ignoreRow )
            srcrow++;
    }
    return cof.determinant();
}

template< typename T > 
T Matrix4< T >::cofactor( const size_t row0, const size_t row1,
                                const size_t row2, const size_t col0,
                                const size_t col1, const size_t col2 ) const
{
    Matrix3< T > cof( m[row0][col0], m[row0][col1], m[row0][col2], 
                      m[row1][col0], m[row1][col1], m[row1][col2],
                      m[row2][col0], m[row2][col1], m[row2][col2] );

    return cof.determinant();
}


template< typename T > 
Matrix4< T > Matrix4< T >::operator-() const
{
    Matrix4< T > result( *this );
    result *= -1.0;
    return result;
}

template< typename T > 
Matrix4< T > Matrix4< T >::negate() const
{
    Matrix4< T > result( *this );
    result *= -1.0;
    return result;
}
}
#endif
