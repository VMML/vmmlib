/*
 * VMMLib - Filter classes
 *
 * @author Jafet Villafranca
 *
 * lopass_filter is a vector data structure in which it is possible to store
 * data values in order to get an smoothed output value
 *
 */

#ifndef __VMML__LOPASS_FILTER__HPP__
#define __VMML__LOPASS_FILTER__HPP__

#include <deque>

namespace vmml
{
    template< size_t M, typename T >
    class lopass_filter
    {
    public:
        lopass_filter( float F )
            : _smooth_factor(F)
        {}

        ~lopass_filter() {}

        const T filter();

        const T add_value( const T& value );

    private:
        std::deque< T > _data;
        const float _smooth_factor;

    }; // class lopass_filter


    template< size_t M, typename T >
    const T
    lopass_filter< M, T >::filter()
    {
        T filtered;
        double weight = 1.0;

        for (size_t i = _data.size(); i --> 0;)
        {
            filtered = filtered * (1 - weight) + _data[i] * weight;
            weight *= _smooth_factor;
        }

        return filtered;

    }

    template< size_t M, typename T >
    const T
    lopass_filter< M, T >::add_value( const T& value )
    {
        _data.push_back(value);

        if (_data.size() > M)
            _data.pop_front();

        return filter();

    }


} // namespace vmml

#endif

