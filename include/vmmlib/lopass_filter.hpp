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

#include <vmmlib/vector.hpp>

namespace vmml
{
    template< size_t M, typename T >
    class lopass_filter
    {
    public:
        lopass_filter( vector< M, T > v, float F )
            : _smooth_factor(F)
        {
            _data.set(v);
        }

        ~lopass_filter() {}

        const T filter();

    private:
        vector< M, T > _data;
        const float _smooth_factor;

    }; // class lopass_filter


    template< size_t M, typename T >
    const T
    lopass_filter< M, T >::filter()
    {
        T filtered;
        double weight = 1.0;

        typedef typename vector< M, T >::const_iterator const_v_iterator;

        const_v_iterator it;
        for( it = _data.end()-1; it != _data.begin()-1; --it )
        {
            filtered = filtered * (1 - weight) + *it * weight;
            weight *= _smooth_factor;
        }

        return filtered;

    }


} // namespace vmml

#endif

