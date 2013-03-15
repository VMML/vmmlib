/*
 * VMMLib - Filter classes
 *
 * @author Jafet Villafranca
 *
 * Implements low pass filtering on a templated data type, which needs to
 * implement multiplication by a scalar float for smoothing.
 *
 */

#ifndef __VMML__LOWPASS_FILTER__HPP__
#define __VMML__LOWPASS_FILTER__HPP__

#include <deque>

namespace vmml
{
template< size_t M, typename T > class lowpass_filter
{
public:
    lowpass_filter( const float F ) : _smooth_factor(F) {}
    ~lowpass_filter() {}

    T get();
    T add( const T& value );
    void set_smooth_factor( const float& f );

private:
    std::deque< T > _data;
    float _smooth_factor;

}; // class lowpass_filter


template< size_t M, typename T > T lowpass_filter< M, T >::get()
{
    if( _data.empty( ))
        return T();

    typename std::deque< T >::const_iterator i = _data.begin();
    T filtered = *i;
    double weight = _smooth_factor;

    for( ++i ; i != _data.end(); ++i )
    {
        filtered = filtered * (1 - weight) + (*i) * weight;
        weight *= _smooth_factor;
    }

    return filtered;
}

template< size_t M, typename T > T lowpass_filter< M, T >::add( const T& value )
{
    _data.push_front( value );

    while( _data.size() > M )
        _data.pop_back();

    return get();
}

template< size_t M, typename T >
void lowpass_filter< M, T>::set_smooth_factor( const float& f )
{
    _smooth_factor = f;
}

} // namespace vmml

#endif
