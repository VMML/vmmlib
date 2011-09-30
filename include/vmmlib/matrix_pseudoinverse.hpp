#ifndef __VMML__MATRIX_PSEUDOINVERSE__HPP__
#define __VMML__MATRIX_PSEUDOINVERSE__HPP__


#include <vmmlib/matrix.hpp>
#include <vmmlib/vector.hpp>
#include <cstddef>
#include <functional>
#include <vmmlib/lapack_svd.hpp>
#include <vmmlib/blas_dgemm.hpp>
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
		
		typedef double float_t;
		
		/// do pseudo inverse for M >= N ///
		void operator()( const T& input, T& pseudoinverse_transposed, 
						typename T::value_type tolerance = std::numeric_limits< typename T::value_type >::epsilon() )
		{

			
			if ( T::ROWS < T::COLS ) {
				VMMLIB_ERROR( "matrix compute_pseudoinverse - number of matrix rows have to be greater or equal to number of matrix columns.", VMMLIB_HERE );
			}
			
			typedef matrix< T::ROWS, T::COLS, float_t > matrix_mn_t;
			typedef matrix< T::COLS, T::ROWS, float_t > matrix_nm_t;
			typedef matrix< T::COLS, T::COLS, float_t > matrix_nn_t;
			typedef vector< T::COLS, float_t > vec_n_t;
			typedef vector< T::ROWS, float_t > vec_m_t;
			typedef matrix< T::COLS, T::ROWS, typename T::value_type> pinv_type;
			typedef lapack_svd< T::ROWS, T::COLS, float_t > svd_type;
			typedef blas_dgemm< T::COLS, 1, T::ROWS, float_t > blas_type;
						
			// perform an SVD on the matrix to get the singular values
			svd_type* svd = new svd_type;
			matrix_mn_t* U = new matrix_mn_t;
			vec_n_t sigmas;
			matrix_nn_t* Vt = new matrix_nn_t;
			matrix_mn_t* input_data =  new matrix_mn_t;
			input_data->cast_from( input );
			
			bool ok = svd->compute( *input_data, *U, sigmas, *Vt ); 
						
			
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
			matrix_nm_t* result = new matrix_nm_t;
			pinv_type* pseudoinverse = new pinv_type;
			result->zero();
			matrix_nm_t* tmp = new matrix_nm_t;
			sigmas.reciprocal();
			//double sigma_inv = 0;

			vec_n_t vt_i;
			vec_m_t u_i;
			blas_type blas_dgemm1;

			if ( num_sigmas >= 1 ) {
				
				it = sigmas.begin();
				for( size_t i = 0 ;  i < num_sigmas && it != it_end; ++it, ++i ) 
				{
					Vt->get_row( i, vt_i);
					U->get_column( i, u_i );

					blas_dgemm1.compute_vv_outer( vt_i, u_i, *tmp );
					
					//sigma value inverted: 1 / *it;
					//sigma_inv = (1 / *it );
					*tmp *= *it ;
					*result += *tmp;
					
				}
				pseudoinverse->cast_from( *result );
				pseudoinverse->transpose_to( pseudoinverse_transposed );
				
			} else {
				pseudoinverse_transposed.zero(); //return matrix with zeros
			}
			
			delete result;
			delete tmp;
			delete Vt;
			delete U;
			delete pseudoinverse;
			delete svd;
			delete input_data;
		}
		
	}; //end compute_pseudoinverse class
}// end vmml namespace

#endif
