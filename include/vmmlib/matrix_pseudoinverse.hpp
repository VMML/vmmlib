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

			
			if ( T::ROWS < T::COLS ) {
				VMMLIB_ERROR( "matrix compute_pseudoinverse - number of matrix rows have to be greater or equal to number of matrix columns.", VMMLIB_HERE );
			}
			
						
			// perform an SVD on the matrix to get the singular values
			lapack_svd< T::ROWS, T::COLS, float_t > svd;
			matrix< T::ROWS, T::COLS, float_t > U;
			vector< T::COLS, float_t > sigmas;
			matrix< T::COLS, T::COLS, float_t > Vt;
			matrix< T::ROWS, T::COLS, float_t > input_data;
			input_data.cast_from( input );
			
			bool ok = svd.compute( input_data, U, sigmas, Vt ); 
						
			
			if ( !ok ) {
				VMMLIB_ERROR( "matrix compute_pseudoinverse - problem with lapack svd.", VMMLIB_HERE );
			}
			/*std::cout << "U: " << std::endl << U << std::endl
			<< " sigmas: " << std::endl << sigmas << std::endl
			<< " Vt: " << std::endl << Vt << std::endl;*/

			// get the number of significant singular, i.e., values which are above the tolerance value
			typename vector< T::COLS, float_t >::const_iterator it = sigmas.begin() , it_end = sigmas.end();
			size_t num_sigmas = 0;
			for( ; it != it_end; ++it )
			{
				if ( *it >= tolerance )
					++num_sigmas;
				else 
					return;
			}
			
			//compute inverse with all the significant inverse singular values
			matrix< T::COLS, T::ROWS, float_t > result;
			matrix< T::COLS, T::ROWS, typename T::value_type> pseudoinverse;
			result.zero();
			matrix< T::COLS, T::ROWS, float_t > tmp;
			double sigma_inv = 0;
			if ( num_sigmas >= 1 ) {
				
				it = sigmas.begin();
				for( size_t i = 0 ;  i < num_sigmas && it != it_end; ++it, ++i ) {
					matrix< 1, T::COLS, float_t > vt_i;
					Vt.get_sub_matrix(vt_i, i);
					matrix< T::COLS, 1, float_t > v_i = transpose(vt_i);
					matrix< T::ROWS, 1, float_t > u_i;
					U.get_sub_matrix(u_i, 0, i); 
					matrix< 1, T::ROWS, float_t > ut_i = transpose(u_i);
					
					//build outer product of v1_i and ut_i
					tmp.multiply(v_i, ut_i);
					
					//sigma value inverted: 1 / *it;
					//sigma_inv = (1 / *it );
					tmp *= *it ;
					result += tmp;
					
				}
				pseudoinverse.cast_from( result );
				pseudoinverse_transposed = transpose( pseudoinverse );
				
			} else {
				pseudoinverse_transposed.zero(); //return matrix with zeros
			}
		}
		
	}; //end compute_pseudoinverse class
}// end vmml namespace

#endif
