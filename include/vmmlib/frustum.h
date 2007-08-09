/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Stefan Eilemann
*
* @license revised BSD license, check LICENSE
*/ 

#ifndef __VMML__FRUSTUM__H__
#define __VMML__FRUSTUM__H__


// - declaration -

namespace vmml
{

template< class T > 
class Frustum
{
public:
    union
    {
        struct
        {
            T left;
            T right;
            T bottom;
            T top;
            T nearPlane;
            T farPlane;
        };
        T data[6];
    };
            
    // contructors
    Frustum(); // warning: components NOT initialised ( for performance )
    Frustum( const T left, const T right, const T bottom, const T top,
             const T nearPlane, const T farPlane );
        
    //the pointer 'values' must be a valid 6 component c array of the resp. type
    Frustum( const float* values );
    Frustum( const double* values );

    ~Frustum();

    void set( const T _left, const T _right, const T _bottom, 
        const T _top, const T _near, const T _far );
             
    Matrix4< T > computeMatrix() const;
    Matrix4< T > computeOrthoMatrix() const;
    void       adjustNear( const T nearPlane );

    friend std::ostream& operator << ( std::ostream& os, const Frustum& frustum)
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << "[" << std::setw(10) << frustum.left << " " 
           << std::setw(10) << frustum.right  << " " 
           << std::setw(10) << frustum.bottom << " " 
           << std::setw(10) << frustum.top    << " " 
           << std::setw(10) << frustum.nearPlane   << " " 
           << std::setw(10) << frustum.farPlane    << "]";
        os.precision( prec );
        os.setf( flags );
        return os;
    };  

    static const Frustum DEFAULT;
};

    
typedef Frustum< float >  Frustumf;
typedef Frustum< double > Frustumd;

typedef Frustum< float >  frustumf;
typedef Frustum< double > frustumd;

} // namespace vmml

// - implementation - //
#include <vmmlib/matrix4.h>

namespace vmml
{

template< typename T > 
const Frustum< T > Frustum< T >::DEFAULT( -1.0, 1.0, -1.0, 1.0, 0.1, 100.0 );



template < class T > 
Frustum< T >::Frustum() 
{} 



template < class T > 
Frustum<T>::Frustum( const T _left, const T _right, const T _bottom, 
                     const T _top, const T _near, const T _far )
    : left( _left ),
      right( _right ),
      bottom( _bottom ),
      top( _top ),
      nearPlane( _near ),
      farPlane( _far )
{} 



template < class T > 
Frustum< T >::Frustum( const float* values )
{
    assert( values && 
            "Frustum: Nullpointer argument as source for initialisation!" );
    left   = static_cast< T > ( values[0] );
    right  = static_cast< T > ( values[1] );
    bottom = static_cast< T > ( values[2] );
    top    = static_cast< T > ( values[3] );
    nearPlane   = static_cast< T > ( values[4] );
    farPlane    = static_cast< T > ( values[5] );
}



template < class T > 
Frustum< T >::Frustum( const double* values )
{
    assert( values &&
            "Frustum: Nullpointer argument as source for initialisation!" );
    left   = static_cast< T > ( values[0] );
    right  = static_cast< T > ( values[1] );
    bottom = static_cast< T > ( values[2] );
    top    = static_cast< T > ( values[3] );
    nearPlane   = static_cast< T > ( values[4] );
    farPlane    = static_cast< T > ( values[5] );
}



template < class T > 
Frustum< T >::~Frustum()
{}



template < class T > 
void 
Frustum< T >::set( const T _left, const T _right, const T _bottom, 
    const T _top, const T _near, const T _far )
{
    left = _left;
    right = _right;
    bottom = _bottom;
    top = _top;
    nearPlane = _near;
    farPlane = _far;
}



template < class T > 
void 
Frustum<T>::adjustNear( const T _near )
{
    if( _near == nearPlane )
        return;

    const T ratio_2 = 0.5 * _near / nearPlane;

    const T hMiddle = (right + left) * 0.5;
    const T width_2 = (right - left) * ratio_2;
    right = hMiddle + width_2;
    left  = hMiddle - width_2; 

    const T vMiddle  = (top + bottom) * 0.5;
    const T height_2 = (top - bottom) * ratio_2;
    top    = vMiddle + height_2;
    bottom = vMiddle - height_2;

    nearPlane = _near;
}



template < class T > 
Matrix4<T> Frustum<T>::computeMatrix() const
{
    Matrix4<T> matrix;

    matrix.m00 = 2.0 * nearPlane / (right - left);
    matrix.m01 = 0.0;
    matrix.m02 = (right + left) / (right - left);
    matrix.m03 = 0.0;
    
    matrix.m10 = 0.0;
    matrix.m11 = 2.0 * nearPlane / (top - bottom);
    matrix.m12 = (top + bottom) / (top - bottom);
    matrix.m13 = 0.0;

    matrix.m20 = 0.0;
    matrix.m21 = 0.0;
    // NOTE: Some glFrustum man pages say wrongly '(far + near) / (far - near)'
    matrix.m22 = -(farPlane + nearPlane) / (farPlane - nearPlane);
    matrix.m23 = -2.0 * farPlane * nearPlane / (farPlane - nearPlane);

    matrix.m30 = 0.0;
    matrix.m31 = 0.0;
    matrix.m32 = -1.0;
    matrix.m33 =  0.0;

    return matrix;
}



template < typename T > 
Matrix4< T > 
Frustum< T >::computeOrthoMatrix() const
{
	return Matrix4< T > (  
        2.0 / ( right - left ), 
        0.0f , 
        0.0f, 
        -( right + left ) / ( right - left ), 

        0.0f, 
        2.0 / ( top-bottom ),
        0.0f,

        -( top + bottom ) / ( top - bottom ), 
        0.0f, 
        0.0f, 
        2.0 / ( nearPlane-farPlane ), 

        -( farPlane + nearPlane ) / ( farPlane - nearPlane ),
        0.0f, 
        0.0f, 
        0.0f, 
        1.0f 
       );
}	


}; //namespace vmml

#endif
