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
            T m00, m01, m02, m03, 
                 m10, m11, m12, m13, 
                 m20, m21, m22, m23, 
                 m30, m31, m32, m33;
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

    inline T* operator[] ( size_t row ) const;
    Vector4< T > getRow( size_t row ) const;
    Vector4< T > getColumn( size_t col ) const; 
    
    // vec3
    void setRow( size_t row, Vector3< T > rowvec );
    void setColumn( size_t col, Vector3< T > colvec );
    // vec4
    void setRow( size_t row, Vector4< T > rowvec );
    void setColumn( size_t col, Vector4< T > colvec );

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
         
template< typename T > 
Matrix4< T >::Matrix4( const Vector4< T >& v0, const Vector4< T >& v1, 
                          const Vector4< T >& v2, const Vector4< T >& v3,
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
    bool equal = true;
    for ( size_t i = 0; i < 16 && equal; ++i )
    {
        equal = ( ml[i] == mm.ml[i] );
    }
    return equal;
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


// returns a row vector
template< typename T > 
T* Matrix4< T >::operator[] (size_t row) const
{
    assert( row < 4 && "Matrix4: Requested Row ( operator[] ) with invalid index!" );
    return const_cast< T* > ( m[row] );
}

template< typename T > 
Vector4< T > Matrix4< T >::getRow( size_t row ) const
{
    assert( row < 4 && "Matrix4: Requested Row ( getRow ) with invalid index!" );
    return Vector4< T >( m[row] );
}

template< typename T > 
Vector4< T > Matrix4< T >::getColumn( size_t col ) const
{
    assert( col < 4 && "Matrix4: Requested Column ( getColumn ) with invalid index!" );
    return Vector4< T > ( m[0+col], m[4+col], m[8+col], m[12+col] );
}

template< typename T > 
void Matrix4< T >::setRow( size_t row, Vector4< T > rowvec )
{
    m[row][0] = rowvec[0];
    m[row][1] = rowvec[1];
    m[row][2] = rowvec[2];
    m[row][3] = rowvec[3];
}

template< typename T > 
void Matrix4< T >::setColumn( size_t col, Vector4< T > colvec )
{
    m[0][col] = colvec[0];
    m[1][col] = colvec[1];
    m[2][col] = colvec[2];
    m[3][col] = colvec[3];
}

template< typename T > 
void Matrix4< T >::setRow( size_t row, Vector3< T > rowvec )
{
    m[row][0] = rowvec[0];
    m[row][1] = rowvec[1];
    m[row][2] = rowvec[2];
}

template< typename T > 
void Matrix4< T >::setColumn( size_t col, Vector3< T > colvec )
{
    m[0][col] = colvec[0];
    m[1][col] = colvec[1];
    m[2][col] = colvec[2];
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
    for ( size_t i = 0; i < 16; ++i )
        result.ml[i] = ml[i] - mm.ml[i];
    return result;
}

template< typename T > 
Matrix4< T > Matrix4< T >::operator* (const Matrix4< T >& mm) const
{
    Matrix4< T > result;
    size_t i, j, k;
    T tmp;

    for (j = 0; j < 4; j++)
        for (i = 0; i < 4; i++) 
        {
            tmp = 0.0;
            for (k = 0; k < 4; k++)
                tmp += m[j][k] * mm.m[k][i];
            result[j][i] = tmp;
        }
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
Matrix4< T >& Matrix4< T >::operator*= ( const Matrix4& mm ) 
{
    for( int j = 0; j < 4; ++j )
        for( int i = 0; i < 4; ++i )
        {
            T tmp = 0.0;
            for( int k = 0; k < 4; ++k )
                tmp += m[j][k] * mm.m[k][i];
            m[j][i] = tmp;
        }
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
	Vector4< T > result;
	result[0] = vv[0] * m[0][0] + vv[1] * m[1][0] + vv[2] * m[2][0] + vv[3] * m[3][0];
	result[1] = vv[0] * m[0][1] + vv[1] * m[1][1] + vv[2] * m[2][1] + vv[3] * m[3][1];
	result[2] = vv[0] * m[0][2] + vv[1] * m[1][2] + vv[2] * m[2][2] + vv[3] * m[3][2];
	result[3] = vv[0] * m[0][3] + vv[1] * m[1][3] + vv[2] * m[2][3] + vv[3] * m[3][3];
	return result;
}

template< typename T > 
Vector3< T > Matrix4< T >::operator* (const Vector3< T >& vv) const
{
	Vector3< T > result;
	result[0] = vv[0] * m[0][0] + vv[1] * m[1][0] + vv[2] * m[2][0];
	result[1] = vv[0] * m[0][1] + vv[1] * m[1][1] + vv[2] * m[2][1];
	result[2] = vv[0] * m[0][2] + vv[1] * m[1][2] + vv[2] * m[2][2];
	return result;
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
    return m[0][0] * cofactor( 1, 2, 3, 1, 2, 3 ) 
         - m[0][1] * cofactor( 1, 2, 3, 0, 2, 3 ) 
         + m[0][2] * cofactor( 1, 2, 3, 0, 1, 3 ) 
         - m[0][3] * cofactor( 1, 2, 3, 0, 1, 2 ); 
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

    T temp = m[0][0];
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

template< typename T >
void Matrix4<T>::rotateY( const T angle )
{
    //matrix multiplication: ml = ml * rotation y axis
    const T sinus = sin(angle);
    const T cosin = cos(angle);

    T temp = m[0][1];
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

template< typename T >
void Matrix4<T>::rotateZ( const T angle )
{
    //matrix multiplication: ml = ml * rotation z axis
    const T sinus = sin(angle);
    const T cosin = cos(angle);

    T temp = m[0][0];
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

template< typename T >
void Matrix4<T>::scale( const T scale[3] )
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

template< typename T >
void Matrix4<T>::scaleTranslation( const T scale[3] )
{
    ml[3] *= scale[0];
    ml[7] *= scale[1];
    ml[11] *= scale[2];
}

template< typename T >
void Matrix4<T>::setTranslation( const T x, const T y, const T z )
{
    ml[3] = x;
    ml[7] = y;
    ml[11] = z;
}

template< typename T > 
void Matrix4< T >::tensor( const Vector3< T >& u, 
                              const Vector3< T >& v )
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

template< typename T > 
void Matrix4< T >::tensor( const Vector4< T >& u,
                              const Vector4< T >& v )
{
    int i, j;
    for (j = 0; j < 4; j++)
        for (i = 0; i < 4; i++)
            m[j][i] = u[j] * v[i];

}

template< typename T > 
T Matrix4< T >::cofactor( const size_t ignoreRow, 
                                const size_t ignoreCol ) const
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
