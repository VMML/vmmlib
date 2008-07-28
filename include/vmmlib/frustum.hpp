/* 
* VMMLib - Vector & Matrix Math Lib
*  
* @author Stefan Eilemann
*
* @license revised BSD license, check LICENSE
*/ 

#ifndef __VMML__FRUSTUM__HPP__
#define __VMML__FRUSTUM__HPP__

#include <vmmlib/matrix.hpp>

// - declaration -

namespace vmml
{

template< class T > 
class frustum
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
    frustum(); // warning: components NOT initialised ( for performance )
    frustum( const T left, const T right, const T bottom, const T top,
             const T nearPlane, const T farPlane );
        
    //the pointer 'values' must be a valid 6 component c array of the resp. type
    frustum( const float* values );
    frustum( const double* values );

    ~frustum();

    void set( const T _left, const T _right, const T _bottom, 
        const T _top, const T _near, const T _far );
        
    // set the frustum using the same parameters as gluPerspective. 
    void setPerspective( T fieldOfViewY, T aspectRatio, T nearPlane_, T farPlane );
             
    matrix< 4, 4, T > computeMatrix() const;
    matrix< 4, 4, T > computeOrthoMatrix() const;

    // 'move' the frustum. this function changes the nearPlane, and adjusts the
    // other parameters in a way that the 'perspective pyramid' stays the same.
    void       adjustNear( const T nearPlane );

    friend std::ostream& operator << ( std::ostream& os, const frustum& frustum_ )
    {
        const std::ios::fmtflags flags = os.flags();
        const int                prec  = os.precision();

        os.setf( std::ios::right, std::ios::adjustfield );
        os.precision( 5 );
        os << "[" << std::setw(10) << frustum_.left << " " 
           << std::setw(10) << frustum_.right  << " " 
           << std::setw(10) << frustum_.bottom << " " 
           << std::setw(10) << frustum_.top    << " " 
           << std::setw(10) << frustum_.nearPlane   << " " 
           << std::setw(10) << frustum_.farPlane    << "]";
        os.precision( prec );
        os.setf( flags );
        return os;
    };  

    static const frustum DEFAULT;
};


#if 0
typedef frustum< float >  Frustumf;
typedef frustum< double > Frustumd;

typedef frustum< float >  frustumf;
typedef frustum< double > frustumd;
#endif


} // namespace vmml

// - implementation - //

namespace vmml
{

template< typename T > 
const frustum< T > frustum< T >::DEFAULT( -1.0, 1.0, -1.0, 1.0, 0.1, 100.0 );



template < class T > 
frustum< T >::frustum() 
{} 



template < class T > 
frustum<T>::frustum( const T _left, const T _right, const T _bottom, 
                     const T _top, const T _near, const T _far )
    : left( _left ),
      right( _right ),
      bottom( _bottom ),
      top( _top ),
      nearPlane( _near ),
      farPlane( _far )
{} 



template < class T > 
frustum< T >::frustum( const float* values )
{
    assert( values && 
            "frustum: Nullpointer argument as source for initialisation!" );
    left   = static_cast< T > ( values[0] );
    right  = static_cast< T > ( values[1] );
    bottom = static_cast< T > ( values[2] );
    top    = static_cast< T > ( values[3] );
    nearPlane   = static_cast< T > ( values[4] );
    farPlane    = static_cast< T > ( values[5] );
}



template < class T > 
frustum< T >::frustum( const double* values )
{
    assert( values &&
            "frustum: Nullpointer argument as source for initialisation!" );
    left   = static_cast< T > ( values[0] );
    right  = static_cast< T > ( values[1] );
    bottom = static_cast< T > ( values[2] );
    top    = static_cast< T > ( values[3] );
    nearPlane   = static_cast< T > ( values[4] );
    farPlane    = static_cast< T > ( values[5] );
}



template < class T > 
frustum< T >::~frustum()
{}



template < class T > 
void 
frustum< T >::set( const T _left, const T _right, const T _bottom, 
    const T _top, const T _near, const T _far )
{
    left = _left;
    right = _right;
    bottom = _bottom;
    top = _top;
    nearPlane = _near;
    farPlane = _far;
}


// 'move' the frustum. this function changes the nearPlane, and adjusts the
// other parameters in a way that the 'perspective pyramid' stays the same.
template < class T > 
void
frustum<T>::adjustNear( const T newNear )
{
	if( newNear == nearPlane )
		return;

	const T ratio = newNear / nearPlane;
	right  *= ratio;
	left   *= ratio;
	top    *= ratio;
	bottom *= ratio;
	nearPlane = newNear;
}



// set the frustum using the same parameters as gluPerspective. 
template < class T > 
void
frustum<T>::setPerspective( T fieldOfViewY, T aspectRatio, T nearPlane_, T farPlane_ )
{
    nearPlane   = nearPlane_;
    farPlane    = farPlane_;
    
    top         = tan( 0.5 * fieldOfViewY * M_PI / 180.0 ) * 0.5;
    bottom      = - top;
    
    left        = bottom * aspectRatio;
    right       = top * aspectRatio;    
}



template < class T > 
matrix< 4, 4, T >
frustum<T>::computeMatrix() const
{
    matrix< 4, 4, T > M;

    M( 0,0 ) = 2.0 * nearPlane / (right - left);
    M( 0,1 ) = 0.0;
    M( 0,2 ) = (right + left) / (right - left);
    M( 0,3 ) = 0.0;
    
    M( 1,0 ) = 0.0;
    M( 1,1 ) = 2.0 * nearPlane / (top - bottom);
    M( 1,2 ) = (top + bottom) / (top - bottom);
    M( 1,3 ) = 0.0;

    M( 2,0 ) = 0.0;
    M( 2,1 ) = 0.0;
    // NOTE: Some glfrustum man pages say wrongly '(far + near) / (far - near)'
    M( 2,2 ) = -(farPlane + nearPlane) / (farPlane - nearPlane);
    M( 2,3 ) = -2.0 * farPlane * nearPlane / (farPlane - nearPlane);

    M( 3,0 ) = 0.0;
    M( 3,1 ) = 0.0;
    M( 3,2 ) = -1.0;
    M( 3,3 ) =  0.0;

    return M;
}



template < typename T > 
matrix< 4, 4, T >
frustum< T >::computeOrthoMatrix() const
{
    matrix< 4, 4, T > M;

    M( 0,0 ) = 2.0 / ( right - left );
    M( 0,1 ) = 0.0;
    M( 0,2 ) = 0.0;
    M( 0,3 ) = -( right + left ) / ( right - left );
    
    M( 1,0 ) = 0.0;
    M( 1,1 ) = 2.0 / ( top-bottom );
    M( 1,2 ) = 0.0f;
    M( 1,3 ) = -( top + bottom ) / ( top - bottom );

    M( 2,0 ) = 0.0;
    M( 2,1 ) = 0.0;
    M( 2,2 ) = -2.0 / ( farPlane-nearPlane );
    M( 2,3 ) = -( farPlane + nearPlane ) / ( farPlane - nearPlane );

    M( 3,0 ) = 0.0;
    M( 3,1 ) = 0.0;
    M( 3,2 ) = 0.0;
    M( 3,3 ) = 1.0f;

    return M;
}	


} //namespace vmml

#endif
