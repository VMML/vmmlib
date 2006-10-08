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
* typename T
*/ 

#ifndef _Matrix3_H_
#define _Matrix3_H_

/* 
 *   3x3 Matrix Class
 */ 

#include "Vector3.h"
#include "Vector4.h"

#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <assert.h>

// - declaration -

namespace vmml
{

template< typename T > 
class Matrix3
{
public:
    union
    {
        struct
        {
            T m00, m01, m02, m10, m11, m12, m20, m21, m22;
        };
        T m[3][3]; // rows 
        T ml[9]; // linear
    };
        
    Matrix3();
    Matrix3( const Matrix3& mm );
    Matrix3( T a, T b, T c, 
             T d, T e, T f, 
             T g, T h, T i 
             );
    Matrix3( const Vector3<T>& v0, const Vector3<T>& v1, 
             const Vector3<T>& v2, bool column_vectors = false );

    // dangerous, but implemented to allow easy conversion between 
    // Matrix< float > and Matrix< double >
    //the pointer 'values must be a valid 9 component c array of the resp. type
    Matrix3( float* values );
    Matrix3( double* values );
 
    inline Matrix3& operator= ( const Matrix3& mm );
    inline Matrix3& operator= ( const T r );
    
    inline bool operator== (const Matrix3& mm) const;
    inline bool operator!= (const Matrix3& mm) const;

    void set( const Matrix3& mm );
    // dangerous, but implemented to allow easy conversion between 
    // Matrix< float > and Matrix< double >
    //the pointer 'values must be a valid 9 component c array of the resp. type
    void set( const float* mm );
    void set( const double* mm );

    inline T* operator[] ( size_t row ) const;
    Vector3< T > getRow( size_t row ) const;
    Vector3< T > getColumn( size_t col ) const; 
    
    void setRow( size_t row, Vector3< T > rowvec );
    void setColumn( size_t col, Vector3< T > colvec );

    // arithmetic operations
    Matrix3 operator+ ( const Matrix3& mm ) const;
    Matrix3 operator- ( const Matrix3& mm ) const;
    Matrix3 operator* ( const Matrix3& mm ) const;
    Matrix3 operator* ( T scalar ) const; // matrix = matrix * scalar 

    Matrix3& operator+= ( const Matrix3& mm );
    Matrix3& operator-= ( const Matrix3& mm );
    Matrix3& operator*= ( const Matrix3& mm );
    Matrix3& operator*= ( T scalar ); // matrix = matrix * scalar 

    // vector = matrix * vector
    Vector3< T > operator* ( const Vector3< T >& vv ) const;

    Matrix3 transpose() const;
    T determinant() const;
    bool isPositiveDefinite( T limit = -0.0000000001 );
    
    bool inverse( Matrix3& result, T limit = 0.0000000001 );

    void tensor( const Vector3< T >& u, const Vector3< T >& v ); 

    Matrix3 operator- () const;
    Matrix3 negate () const;

    friend std::ostream& operator << ( std::ostream& os, const Matrix3& m )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << std::endl << "|" << std::setw(7) << m.ml[0] << " " 
           << std::setw(7) << m.ml[1] << " " 
           << std::setw(7) << m.ml[2] << "|" << std::endl
           << "|" << std::setw(7) << m.ml[3] << " "
           << std::setw(7) << m.ml[4] << " " 
           << std::setw(7) << m.ml[5] << "|" << std::endl
           << "|" << std::setw(7) << m.ml[6] << " "
           << std::setw(7) << m.ml[7] << " " 
           << std::setw(7) << m.ml[8] << "|" << std::endl; 
        os.precision( prec );
        os.setf( flags );
        return os;
    };  
    
    static const Matrix3 IDENTITY;
    static const Matrix3 ZERO;
};


#ifndef VMMLIB_DISABLE_TYPEDEFS
    typedef Matrix3<float>  Matrix3f;
    typedef Matrix3<double> Matrix3d;
#endif

// - implementation -

template< typename T > 
const Matrix3< T > Matrix3< T >::IDENTITY( 1, 0, 0, 0, 1, 0, 0, 0, 1 );

template< typename T > 
const Matrix3< T > Matrix3< T >::ZERO( 0, 0, 0, 0, 0, 0, 0, 0, 0 );


template< typename T > 
Matrix3< T >::Matrix3()
{}

template< typename T > 
Matrix3< T >::Matrix3( const Matrix3< T >& mm )
{
    memcpy( m, mm.m, 9 * sizeof( T ));
}

template< typename T > 
Matrix3< T >::Matrix3( 
            T a, T b, T c, 
            T d, T e, T f, 
            T g, T h, T i )
{
    ml[0] = a;
    ml[1] = b;
    ml[2] = c;
    ml[3] = d;
    ml[4] = e;
    ml[5] = f;
    ml[6] = g;
    ml[7] = h;
    ml[8] = i;
}
         
template< typename T > 
Matrix3< T >::Matrix3( const Vector3< T >& v0, const Vector3< T >& v1,
                          const Vector3< T >& v2, bool column_vectors )
{
    if ( column_vectors )
    {
        ml[0] = v0.x;
        ml[1] = v0.y;
        ml[2] = v0.z;
        ml[3] = v1.x;
        ml[4] = v1.y;
        ml[5] = v1.z;
        ml[6] = v2.x;
        ml[7] = v2.y;
        ml[8] = v2.z;
    } 
    else
    {
        ml[0] = v0.x;
        ml[3] = v0.y;
        ml[6] = v0.z;
        ml[1] = v1.x;
        ml[4] = v1.y;
        ml[7] = v1.z;
        ml[2] = v2.x;
        ml[5] = v2.y;
        ml[8] = v2.z;
    }
}

template< typename T > 
Matrix3< T >::Matrix3( float* values )
{
    assert( values && "Matrix3: Initialisation of a Matrix from a Nullpointer was requested." );
    for ( size_t i = 0; i < 9; ++i )
        ml[i] = static_cast< T > ( values[i] );
}

template< typename T > 
Matrix3< T >::Matrix3( double* values )
{
    assert( values && "Matrix3: Initialisation of a Matrix from a Nullpointer was requested." );
    for ( size_t i = 0; i < 9; ++i )
        ml[i] = static_cast< T > ( values[i] );
}

template< typename T > 
Matrix3< T >& Matrix3< T >::operator= ( const T r )
{
    for ( size_t i = 0; i < 9; ++i )
    {
        ml[i] = r;
    }
    return *this;
}


template< typename T > 
Matrix3< T >& Matrix3< T >::operator= ( const Matrix3< T >& mm )
{
    memcpy(ml,mm.ml,9*sizeof(T));
    return *this;
}

template< typename T > 
bool Matrix3< T >::operator== (const Matrix3< T >& mm) const
{
    bool equal = true;
    for ( size_t i = 0; i < 9 && equal; ++i )
    {
        equal = ( ml[i] == mm.ml[i] );
    }
    return equal;

}

template< typename T > 
inline bool Matrix3< T >::operator!= (const Matrix3< T >& mm) const
{
    return !operator==(mm);
}

template< typename T > 
void Matrix3< T >::set( const Matrix3& mm )
{
    memcpy( ml, mm.ml, 9 * sizeof( T ) );
}

template< typename T > 
void Matrix3< T >::set( const float* mm )
{
    assert( mm && "Matrix3: Nullpointer argument as source for initialisation!" );
    for ( size_t i = 0; i < 9; ++i )
        ml[i] = static_cast< T > ( mm[i] );
}

template< typename T > 
void Matrix3< T >::set( const double* mm )
{
    assert( mm && "Matrix3: Nullpointer argument as source for initialisation!" );
    for ( size_t i = 0; i < 9; ++i )
        ml[i] = static_cast< T > ( mm[i] );
}


// returns a row vector
template< typename T > 
T* Matrix3< T >::operator[] ( size_t row ) const
{
    if ( row > 2 ) 
        std::cerr << "Matrix3::op[] - invalid row index " << row << "." << std::endl;
    assert( row < 3 && "Matrix3: Requested Row ( operator[] ) with invalid index!" );
    return const_cast< T* > ( m[row] );
}

template< typename T > 
Vector3< T > Matrix3< T >::getRow( size_t row ) const
{
    if ( row > 2 ) 
        std::cerr << "Matrix3::getRow - invalid row index " << row << "." << std::endl;
    assert( row < 3 && "Matrix3: Requested Row ( getRow ) with invalid index!" );
    return Vector3< T >( m[row] );
}

template< typename T > 
Vector3< T > Matrix3< T >::getColumn( size_t col ) const
{
    assert( col < 3 && "Matrix3: Requested Column ( getColumn ) with invalid index!" );
    return Vector3< T > ( m[0+col], m[3+col], m[6+col] );
}

template< typename T > 
void Matrix3< T >::setRow( size_t row, Vector3< T > rowvec )
{
    m[row][0] = rowvec[0];
    m[row][1] = rowvec[1];
    m[row][2] = rowvec[2];
}

template< typename T > 
void Matrix3< T >::setColumn( size_t col, Vector3< T > colvec )
{
    m[0][col] = colvec[0];
    m[1][col] = colvec[1];
    m[2][col] = colvec[2];
}

template< typename T > 
Matrix3< T > Matrix3< T >::operator+ ( const Matrix3< T >& mm ) const
{
    Matrix3< T > result;
    for ( size_t i = 0; i < 9; ++i ) 
        result.ml[i] = ml[i] + mm.ml[i];
    return result;
}

template< typename T > 
Matrix3< T > Matrix3< T >::operator- ( const Matrix3< T >& mm ) const
{
    Matrix3< T > result;
    for ( size_t i = 0; i < 9; ++i ) 
        result.ml[i] = ml[i] - mm.ml[i];
    return result;
}

template< typename T > 
Matrix3< T > Matrix3< T >::operator* ( const Matrix3< T >& mm ) const
{
    Matrix3< T > result;
    size_t i, j, k;
    T tmp;

    for (j = 0; j < 3; j++)
        for (i = 0; i < 3; i++) 
        {
            tmp = 0.0;
            for (k = 0; k < 3; k++)
                tmp += m.m[j][k] * mm.m[k][i];
            result[j][i] = tmp;
        }
}

template< typename T > 
Matrix3< T > Matrix3< T >::operator* ( T scalar ) const
{
    Matrix3< T > result;
    for ( size_t i = 0; i < 9; ++i )
        result.ml[i] = ml[i] * scalar;
    return result;
}

template< typename T > 
Matrix3< T >& Matrix3< T >::operator+= ( const Matrix3& mm )
{
    for ( size_t i = 0; i < 9; ++i )
        ml[i] += mm.ml[i];
    return *this;
}

template< typename T > 
Matrix3< T >& Matrix3< T >::operator-= ( const Matrix3& mm )
{
    for ( size_t i = 0; i < 9; ++i )
        ml[i] -= mm.ml[i];
    return *this;
}

template< typename T > 
Matrix3< T >& Matrix3< T >::operator*= ( const Matrix3& mm )
{
    size_t i, j, k;
    T tmp;

    for (j = 0; j < 3; j++)
        for (i = 0; i < 3; i++) 
        {
            tmp = 0.0;
            for (k = 0; k < 3; k++)
                tmp += m[j][k] * mm.m[k][i];
            m[j][i] = tmp;
        }
    return *this;
}

template< typename T > 
Matrix3< T >& Matrix3< T >::operator*= ( T scalar ) // matrix = matrix * scalar 
{
    for ( size_t i = 0; i < 9; ++i )
        ml[i] *= scalar;
    return *this;
}

template< typename T > 
Vector3< T > Matrix3< T >::operator* ( const Vector3< T >& vv ) const
{  
    Vector3< T > result;
    for (size_t i = 0; i < 3; ++i)
        result[i] = m[i][0]* vv[0] + m[i][1] * vv[1] + m[i][2] * vv[2];
    return result;
}

template< typename T > 
Matrix3< T > Matrix3< T >::transpose() const
{
    Matrix3< T > result;
    result.m[0][0] = m[0][0];
    result.m[0][1] = m[1][0];
    result.m[0][2] = m[2][0];
    result.m[1][0] = m[0][1];
    result.m[1][1] = m[1][1];
    result.m[1][2] = m[2][1];
    result.m[2][0] = m[0][2];
    result.m[2][1] = m[1][2];
    result.m[2][2] = m[2][2];
    return result;
}

template< typename T > 
T Matrix3< T >::determinant() const
{
    Vector3< T > cof;
    cof[0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    cof[1] = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    cof[2] = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    return m[0][0] * cof[0] + m[0][1] * cof[1] + m[0][2] * cof[2];
}

template< typename T > 
bool Matrix3< T >::isPositiveDefinite( T limit )
{
    Vector3< T > d;
    d[0] = m[0][0];
    d[1] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
    d[2] = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] + 
         m[0][1]*m[1][2]*m[2][0] - m[0][1]*m[1][0]*m[2][2] +
         m[0][2]*m[1][0]*m[2][1] - m[0][2]*m[1][1]*m[2][0];

    if (d[0] < limit || d[1] < limit || d[2] < limit) 
        return false;
    return true;
}


template< typename T > 
bool Matrix3< T >::inverse( Matrix3< T >& result, T limit )
{
    // Invert a 3x3 using cofactors.  This is about 8 times faster than
    // the Numerical Recipes code which uses Gaussian elimination.

    result[0][0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    result[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];
    result[0][2] = m[0][1] * m[1][2] - m[0][2] * m[1][1];
    result[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    result[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];
    result[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];
    result[2][0] = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    result[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];
    result[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];

    T det = m[0][0] * result[0][0] + m[0][1] * result[1][0] + m[0][2] * result[2][0];

    if ( fabs( det ) <= limit )
        return false; // matrix is not invertible

    T detinv = 1.0 / det;
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
            result[i][j] *= detinv;
    }
    return true;

}

template< typename T > 
void Matrix3< T >::tensor( const Vector3< T >& u, const Vector3< T >& v)
{
    int i, j;
    for (j = 0; j < 3; j++)
        for (i = 0; i < 3; i++)
            m[j][i] = u[j] * v[i];
}

template< typename T > 
Matrix3< T > Matrix3< T >::operator-() const
{
    Matrix3< T > result( *this );
    result *= -1.0;
    return result;
}

template< typename T > 
Matrix3< T > Matrix3< T >::negate() const
{
    Matrix3< T > result( *this );
    result *= -1.0;
    return result;
}
}
#endif
