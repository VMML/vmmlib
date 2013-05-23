///* 
// * VMMLib - Tensor Classes
// *  
// * @author Rafael Ballester
// * 
// */
//
//#ifndef __VMML__TENSOR4_ITERATOR__HPP__
//#define __VMML__TENSOR4_ITERATOR__HPP__
//
//#include <vmmlib/matrix.hpp>
//#include <vmmlib/tensor3.hpp>
//
//namespace vmml
//{
//	
//
//template< typename T >
//class tensor4_iterator
//{
//public:
//	
//	typedef typename T::value_type              value_type;
//	typedef typename T::pointer                 pointer;
//	typedef typename T::reference               reference;
//    typedef typename std::forward_iterator_tag  iterator_category;
//    typedef size_t                              difference_type;
//    
//	typedef typename matrix< T::ROWS, T::COLS, typename T::value_type >::iterator matrix_iterator;
//    
//    typedef typename T::front_slice_type  slice_type;
//	
//	tensor4_iterator() : _tensor4( 0 ), _matrix_index( 0 ) {};
//	
//	tensor4_iterator( T& t_, bool begin_ ) : _tensor4( &t_ ), _matrix_index( 0 )
//    {
//        if ( begin_ )
//        {
//            _matrix_index       = 0;
//            slice_type& slice_  = _tensor4->get_frontal_slice_fwd( _matrix_index );
//            
//            _matrix_it          = slice_.begin();
//            _matrix_it_end      = slice_.end();
//        }
//        else
//        {
//// TODO
//        }
//    }
//	
//	// TODO
//    
//protected:
//	matrix_iterator		_matrix_it;
//	matrix_iterator		_matrix_it_end;
//	T*					_tensor4;
//	size_t				_matrix_index;
//
//}; //end tensor4_iterator class
//
//// TODO const_iterator
//
//
//
//
//}// end vmml namespace
//
//#endif
//
