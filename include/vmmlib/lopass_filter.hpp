/*
 * VMMLib - Filter classes
 *
 * @author Jafet Villafranca
 *
 * lopass_filter is a vector data structure in which it is possible to store
 * data values to get an smoothed output value
 *
 */

#ifndef __VMML__LOPASS_FILTER__HPP__
#define __VMML__LOPASS_FILTER__HPP__

#include <deque>

namespace vmml
{
template< size_t M, typename T > class lopass_filter
{
public:
    lopass_filter( const float F ) : _smooth_factor(F) {}
    ~lopass_filter() {}

    T filter();
    T add( const T& value );
    void set_smooth_factor( const float& f );

private:
    std::deque< T > _data;
    float _smooth_factor;

}; // class lopass_filter


template< size_t M, typename T > T lopass_filter< M, T >::filter()
{
    if( _data.empty( ))
        return T(0);

    typename std::deque< T >::const_iterator i = _data.begin();
    T filtered = *i;
    double weight = 1.0;

    for( ++i ; i != _data.end(); ++i )
    {
        filtered = filtered * (1 - weight) + (*i) * weight;
        weight *= _smooth_factor;
    }

    return filtered;
}

template< size_t M, typename T > T lopass_filter< M, T >::add( const T& value )
{
    _data.push_front( value );

    while( _data.size() > M )
        _data.pop_back();

    return filter();
}

template< size_t M, typename T >
void lopass_filter< M, T>::set_smooth_factor( const float& f )
{
    _smooth_factor = f;
}

} // namespace vmml

#endif
