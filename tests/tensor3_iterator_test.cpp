#include "tensor3_iterator_test.hpp"

#include <stdint.h>
#include <vmmlib/tensor3_iterator.hpp>
#include <vmmlib/tensor3.hpp>

#include <sstream>

namespace vmml
{

	bool
	tensor3_iterator_test::run()
	{
        bool global_ok = true;
		bool ok = false;

		tensor3< 2, 3, 4, unsigned short >  t3;

		t3.fill_increasing_values();

        typedef tensor3< 2, 3, 4, unsigned short > myt3;

        std::vector< unsigned short > t3_iter_order;
        std::vector< unsigned short > hand_iter_order;


		for( myt3::iterator i = t3.begin(); i != t3.end(); ++i )
            t3_iter_order.push_back( *i );

        for( size_t index = 0; index < 4; ++index )
        {
            matrix< 2, 3, unsigned short >& m = t3.get_frontal_slice_fwd(index);

            for( tensor3< 2, 3, 4, unsigned short >::matrix_iterator i =
                     m.begin(); i != m.end(); ++i )
            {
                hand_iter_order.push_back( *i );
            }
        }

        TEST(t3_iter_order == hand_iter_order);

		log( "tensor3 iterator  ", ok  );


		return global_ok;
	}


} // namespace vmml
