#ifndef __VMML__MATRIX_PSEUDOINVERSE__HPP__
#define __VMML__MATRIX_PSEUDOINVERSE__HPP__


#include <vmmlib/matrix.hpp>
#include <vmmlib/vector.hpp>
#include <cstddef>
#include <functional>
#include <vmmlib/lapack_svd.hpp>
//#include <vmmlib/enable_if.hpp>

/* 
  *** computes the pseudo inverse of a non-square matrix ***
 - the pseudo inverse is computed by the help of SVD
 - the tolerance for the significant singular values is optionally set
 - implementation works only for matrices with more rows than columns or quadratic 
   matrices. use a transposed input matrix for matrices with more columns than rows
 */

namespace vmml
{
	
	template< typename T >
	class compute_pseudoinverse
	{
		//TODO: Add restriction for matrices with M >= N only to template
				
	public:		
		/// do pseudo inverse for M >= N ///
		void operator()( const T& input, T& pseudoinverse_transposed, 
						typename T::value_type tolerance = std::numeric_limits< typename T::value_type >::epsilon() )
		{
			//const size_t M 			= T::ROWS;
			//const size_t N 			= T::COLS;
			//typedef T::value_type myfloat;		

			if ( T::ROWS < T::COLS ) {
				VMMLIB_ERROR( "matrix compute_pseudoinverse - number of matrix rows have to be greater or equal to number of matrix columns.", VMMLIB_HERE );
			}
			
						
			// perform an SVD on the matrix to get the singular values
			lapack_svd< T::ROWS, T::COLS, typename T::value_type > svd;
			matrix< T::ROWS, T::COLS, typename T::value_type > U;
			vector< T::COLS, typename T::value_type > sigmas;
			matrix< T::COLS, T::COLS, typename T::value_type > Vt;
			bool ok = svd.compute( input, U, sigmas, Vt ); 
			
			if ( !ok ) {
				VMMLIB_ERROR( "matrix compute_pseudoinverse - problem with lapack svd.", VMMLIB_HERE );
			}
		

			// get the number of significant singular, i.e., values which are above the tolerance value
			typename vector< T::COLS, typename T::value_type >::const_iterator it = sigmas.begin() , it_end = sigmas.end();
			size_t num_sigmas = 0;
			for( ; it != it_end; ++it )
			{
				if ( *it >= tolerance )
					++num_sigmas;
				else 
					return;
			}
			
			//compute inverse with all the significant inverse singular values
			matrix< T::COLS, T::ROWS, typename T::value_type > result;
			result.zero();
			matrix< T::COLS, T::ROWS, typename T::value_type > tmp;
			if ( num_sigmas >= 1 ) {
				
				it = sigmas.begin();
				for( size_t i = 0 ;  i < num_sigmas && it != it_end; ++it, ++i ) {
					matrix< 1, T::COLS, typename T::value_type > vt_i;
					Vt.get_sub_matrix(vt_i, i);
					matrix< T::COLS, 1, typename T::value_type > v_i = transpose(vt_i);
					matrix< T::ROWS, 1, typename T::value_type > u_i;
					U.get_sub_matrix(u_i, 0, i); 
					matrix< 1, T::ROWS, typename T::value_type > ut_i = transpose(u_i);
					
					//build outer product of v1_i and ut_i
					tmp.multiply(v_i, ut_i);
					
					//sigma value inverted: 1 / *it;
					tmp *= (1 / *it );
					result += tmp;
					
				}				
				pseudoinverse_transposed = transpose( result );
				
			} else {
				pseudoinverse_transposed.zero(); //return matrix with zeros
			}
		}
		
	}; //end compute_pseudoinverse class
}// end vmml namespace

#endif
