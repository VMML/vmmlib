/*
 * Copyright (c) 2006-2014, Visualization and Multimedia Lab,
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
#ifndef __VMML__AXIS_ALIGNED_BOUNDING_BOX__HPP__
#define __VMML__AXIS_ALIGNED_BOUNDING_BOX__HPP__

#include <vmmlib/vector.hpp>
#include <limits>

namespace vmml
{
/**
 * An axis-aligned bounding box.
 *
 * An empty bounding box has undefined, implementation-specific values. Read
 * operations (getMin(), getMax(), getDimension(), isIn(), etc.) have undefined
 * semantics on an empty bounding box. set() and merge() operations will define
 * the bounding box correctly.
 */
template< typename T > class AxisAlignedBoundingBox
{
public:
    /** Create an empty bounding box. */
    AxisAlignedBoundingBox();
    AxisAlignedBoundingBox( const vector< 3, T >& pMin,
                            const vector< 3, T >& pMax );
    AxisAlignedBoundingBox( const vector< 4, T >& sphere );
    AxisAlignedBoundingBox( T cx, T cy, T cz, T size );

    bool isIn( const vector< 3, T >& pos );
    bool isIn2d( const vector< 3, T >& pos ); // only x and y components are checked
    bool isIn( const vector< 4, T >& sphere );

    void set( const vector< 3, T >& pMin, const vector< 3, T >& pMax );
    void set( T cx, T cy, T cz, T size );
    void setMin( const vector< 3, T >& pMin );
    void setMax( const vector< 3, T >& pMax );

    const vector< 3, T >& getMin() const;
    const vector< 3, T >& getMax() const;
    vector< 3, T >& getMin();
    vector< 3, T >& getMax();

    void merge( const AxisAlignedBoundingBox< T >& aabb );
    void merge( const vector< 3, T >& point );

    void setEmpty();
    bool isEmpty() const;

    AxisAlignedBoundingBox operator*( const T value ) const;
    AxisAlignedBoundingBox operator/( const T value ) const;
    AxisAlignedBoundingBox operator+( const T value ) const;
    AxisAlignedBoundingBox operator-( const T value ) const;

    void operator*=( const T value );
    void operator/=( const T value );
    void operator+=( const T value );
    void operator-=( const T value );

    template< class U >
    bool operator==( const AxisAlignedBoundingBox< U >& other ) const;
    template< class U >
    bool operator!=( const AxisAlignedBoundingBox< U >& other ) const;

    vector< 3, T > getCenter() const;
    vector< 3, T > getDimension() const;

    static AxisAlignedBoundingBox< T > makeUnitBox();

protected:
    vector< 3, T > _min;
    vector< 3, T > _max;
};

template< typename T >
inline std::ostream& operator << ( std::ostream& os,
                                   const AxisAlignedBoundingBox< T >& aabb )
{
    return os << aabb.getMin() << " - " << aabb.getMax();
}

template< typename T >
AxisAlignedBoundingBox< T >::AxisAlignedBoundingBox()
    : _min( std::numeric_limits< T >::max( ))
    , _max( -std::numeric_limits< T >::max( ))
{}

template< typename T >
AxisAlignedBoundingBox< T >::AxisAlignedBoundingBox( const vector< 3, T >& pMin,
                                                     const vector< 3, T >& pMax)
    : _min( pMin )
    , _max( pMax )
{}

template< typename T >
AxisAlignedBoundingBox< T >::AxisAlignedBoundingBox( const vector< 4, T >& sphere )
{
    _max = _min = sphere.getCenter();
    _max += sphere.getRadius();
    _min -= sphere.getRadius();
}

template< typename T >
AxisAlignedBoundingBox< T >::AxisAlignedBoundingBox( T cx, T cy, T cz, T size )
{
    _max = _min = vector< 3, T >( cx, cy, cz );
    _max += size;
    _min -= size;
}

template< typename T >
inline bool AxisAlignedBoundingBox< T >::isIn( const vector< 4, T >& sphere )
{
    vector< 3, T > sv ( sphere.getCenter() );
    sv += sphere.getRadius();
    if ( sv.x() > _max.x() || sv.y() > _max.y() || sv.z() > _max.z() )
        return false;
    sv -= sphere.getRadius() * 2.0f;
    if ( sv.x() < _min.x() || sv.y() < _min.y() || sv.z() < _min.z() )
        return false;
    return true;
}

template< typename T >
inline bool AxisAlignedBoundingBox< T >::isIn( const vector< 3, T >& pos )
{
    if ( pos.x() > _max.x() || pos.y() > _max.y() || pos.z() > _max.z() ||
         pos.x() < _min.x() || pos.y() < _min.y() || pos.z() < _min.z( ))
    {
        return false;
    }
    return true;
}

template< typename T >
inline bool AxisAlignedBoundingBox< T >::isIn2d( const vector< 3, T >& pos )
{
    if ( pos.x() > _max.x() || pos.y() > _max.y() || pos.x() < _min.x() ||
         pos.y() < _min.y( ))
    {
        return false;
    }
    return true;
}

template< typename T >
inline void AxisAlignedBoundingBox< T >::set( const vector< 3, T >& pMin,
    const vector< 3, T >& pMax )
{
    _min = pMin;
    _max = pMax;
}

template< typename T >
inline void AxisAlignedBoundingBox< T >::set( T cx, T cy, T cz, T size )
{
    vector< 3, T > center( cx, cy, cz );
    _min = center - size;
    _max = center + size;
}

template< typename T >
inline void AxisAlignedBoundingBox< T >::setMin( const vector< 3, T >& pMin )
{
    _min = pMin;
}

template< typename T >
inline void AxisAlignedBoundingBox< T >::setMax( const vector< 3, T >& pMax )
{
    _max = pMax;
}

template< typename T >
inline const vector< 3, T >& AxisAlignedBoundingBox< T >::getMin() const
{
    return _min;
}

template< typename T >
inline const vector< 3, T >& AxisAlignedBoundingBox< T >::getMax() const
{
    return _max;
}

template< typename T > inline vector< 3, T >& AxisAlignedBoundingBox< T >::getMin()
{
    return _min;
}

template< typename T > inline vector< 3, T >& AxisAlignedBoundingBox< T >::getMax()
{
    return _max;
}

template< typename T > AxisAlignedBoundingBox< T >
AxisAlignedBoundingBox< T >::operator*( const T value ) const
{
    AxisAlignedBoundingBox result = *this;
    result *= value;
    return result;
}

template< typename T > AxisAlignedBoundingBox< T >
AxisAlignedBoundingBox< T >::operator/( const T value ) const
{
    AxisAlignedBoundingBox result = *this;
    result /= value;
    return result;
}

template< typename T > AxisAlignedBoundingBox< T >
AxisAlignedBoundingBox< T >::operator+( const T value ) const
{
    AxisAlignedBoundingBox result = *this;
    result += value;
    return result;
}

template< typename T > AxisAlignedBoundingBox< T >
AxisAlignedBoundingBox< T >::operator-( const T value ) const
{
    AxisAlignedBoundingBox result = *this;
    result -= value;
    return result;
}

template< typename T >
void AxisAlignedBoundingBox< T >::operator*=( const T value )
{
    _min *= value;
    _max *= value;
}

template< typename T >
void AxisAlignedBoundingBox< T >::operator/=( const T value )
{
    _min /= value;
    _max /= value;
}

template< typename T >
void AxisAlignedBoundingBox< T >::operator+=( const T value )
{
    _min += value;
    _max += value;
}

template< typename T >
void AxisAlignedBoundingBox< T >::operator-=( const T value )
{
    _min -= value;
    _max -= value;
}

template< typename T > template< class U > bool
AxisAlignedBoundingBox< T >::operator==( const AxisAlignedBoundingBox< U >& other )
    const
{
    return _min == other._min && _max == other._max;
}

template< typename T > template< class U > bool
AxisAlignedBoundingBox< T >::operator!=( const AxisAlignedBoundingBox< U >& other )
    const
{
    return _min != other._min || _max != other._max;
}

template< typename T >
vector< 3, T > AxisAlignedBoundingBox< T >::getCenter() const
{
    return _min + ( ( _max - _min ) * 0.5f );
}

template< typename T >
vector< 3, T > AxisAlignedBoundingBox< T >::getDimension() const
{
    return _max - _min;
}

template< typename T >
void AxisAlignedBoundingBox< T >::merge( const AxisAlignedBoundingBox<T>& aabb )
{
    const vector< 3, T >& min = aabb.getMin();
    const vector< 3, T >& max = aabb.getMax();

    if ( min.x() < _min.x() )
        _min.x() = min.x();
    if ( min.y() < _min.y() )
        _min.y() = min.y();
    if ( min.z() < _min.z() )
        _min.z() = min.z();

    if ( max.x() > _max.x() )
        _max.x() = max.x();
    if ( max.y() > _max.y() )
        _max.y() = max.y();
    if ( max.z() > _max.z() )
        _max.z() = max.z();
}

template< typename T >
void AxisAlignedBoundingBox< T >::merge( const vector< 3, T >& point )
{
    if ( point.x() < _min.x() )
        _min.x() = point.x();
    if ( point.y() < _min.y() )
        _min.y() = point.y();
    if ( point.z() < _min.z() )
        _min.z() = point.z();

    if ( point.x() > _max.x() )
        _max.x() = point.x();
    if ( point.y() > _max.y() )
        _max.y() = point.y();
    if ( point.z() > _max.z() )
        _max.z() = point.z();
}

template< typename T >inline
void AxisAlignedBoundingBox< T >::setEmpty()
{
    _min = std::numeric_limits< T >::max();
    _max = -std::numeric_limits< T >::max();
}


template< typename T > inline bool AxisAlignedBoundingBox< T >::isEmpty() const
{
    return ( _min.x() >=  _max.x() || _min.y() >=  _max.y() ||
             _min.z() >=  _max.x( ));
}

template< typename T >
AxisAlignedBoundingBox< T > AxisAlignedBoundingBox< T >::makeUnitBox()
{
    return AxisAlignedBoundingBox( vector< 3, T >::ZERO, vector< 3, T >::ONE );
}

}; //namespace vmml

#endif
