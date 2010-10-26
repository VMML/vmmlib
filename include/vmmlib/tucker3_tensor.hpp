/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 * @author Jonas Boesch
 *
 * The tucker3 tensor class is consists of the same components (core tensor, basis matrices u1-u3) as the tucker3 model described in:
 * - Tucker,  “Some mathematical notes on three-mode factor analysis”, Psychometrika, vol. 31, no. 3, pp. 279–311., 1966 Sep.
 * - De Lathauwer L., De Moor B., Vandewalle J., ``A multilinear singular value decomposition'', 
 * SIAM J. Matrix Anal. Appl., vol. 21, no. 4, Apr. 2000, pp. 1253-1278.
 * - De Lathauwer L., De Moor B., Vandewalle J., ``On the Best rank-1 and Rank-$(R_1, R_2, ..., R_N)$ Approximation and Applications of Higher-Order Tensors'', 
 * SIAM J. Matrix Anal. Appl., vol. 21, no. 4, Apr. 2000, pp. 1324-1342.
 * - T. G. Kolda and B. W. Bader. Tensor Decompositions and Applications. 
 * SIAM Review, Volume 51, Number 3, Pages 455-500, September 2009.
 * 
 */

#ifndef __VMML__TUCKER3_TENSOR__HPP__
#define __VMML__TUCKER3_TENSOR__HPP__

#include <vmmlib/tensor3.hpp>
#include <vmmlib/tensor3_iterator.hpp>
#include <vmmlib/lapack_svd.hpp>
#include <vmmlib/matrix_pseudoinverse.hpp>

namespace vmml
{
	
	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T_value = float, typename T_coeff = float >
	class tucker3_tensor
	{
public:    
    tucker3_tensor( tensor3< R1, R2, R3, T_coeff >& core, matrix< I1, R1, T_coeff >& U1, matrix< I2, R2, T_coeff >& U2, matrix< I3, R3, T_coeff >& U3 );
	tucker3_tensor();

	typedef tensor3< I1, I2, I3, T_value > t3_type;
	typedef typename t3_type::iterator t3_iterator;
	typedef typename t3_type::const_iterator t3_const_iterator;
		
	typedef tensor3< I1, I2, I3, T_coeff > t3_coeff_type;
	typedef typename t3_coeff_type::iterator t3_coeff_iterator;
	typedef typename t3_coeff_type::const_iterator t3_coeff_const_iterator;
		
	typedef tensor3< R1, R2, R3, T_coeff > t3_core_type;
	typedef typename t3_core_type::iterator t3_core_iterator;
	typedef typename t3_core_type::const_iterator t3_core_const_iterator;
		
	typedef matrix< I1, R1, T_coeff > u1_type;
	typedef typename u1_type::iterator u1_iterator;
	typedef typename u1_type::const_iterator u1_const_iterator;

	typedef matrix< I2, R2, T_coeff > u2_type;
	typedef typename u2_type::iterator u2_iterator;
	typedef typename u2_type::const_iterator u2_const_iterator;
	
	typedef matrix< I3, R3, T_coeff > u3_type;
	typedef typename u3_type::iterator u3_iterator;
	typedef typename u3_type::const_iterator u3_const_iterator;
	
	//matrix types for inverted (pseudo-inverted) u1-u3
	typedef matrix< R1, I1, T_coeff > u1_inv_type;
	typedef matrix< R2, I2, T_coeff > u2_inv_type;
	typedef matrix< R3, I3, T_coeff > u3_inv_type;
		
	typedef matrix< I1, I2*I3, T_coeff > mode1_matricization_type;
	typedef matrix< I2, I1*I3, T_coeff > mode2_matricization_type;
	typedef matrix< I3, I1*I2, T_coeff > mode3_matricization_type;

	static const size_t SIZE = R1*R2*R3 + I1*R1 + I2*R2 + I3*R3;
		
	void set_core( const t3_core_type& core )  { _core = core; } ;
	void set_u1( const u1_type& U1 ) { _u1 = U1; } ;
	void set_u2( const u2_type& U2 ) { _u2 = U2; } ;
	void set_u3( const u3_type& U3 ) { _u3 = U3; } ;
	
    void get_core( t3_core_type& data_ ) const { data_ = _core; } ;
	void get_u1( u1_type& U1 ) const { U1 = _u1; } ;
	void get_u2( u2_type& U2 ) const { U2 = _u2; } ;
	void get_u3( u3_type& U3 ) const { U3 = _u3; } ;
		
	void export_to( std::vector< T_coeff >& data_ ) const;
	void import_from( const std::vector< T_coeff >& data_ );	
	
	void reconstruct( t3_type& data_ ) const;
	void decompose( const t3_type& data_ ); 
		
	/* derive core
	   implemented accodring to core = data x_1 U1_pinv x_2 U2_pinv x_3 U3_pinv, 
	   where x_1 ... x_3 are n-mode products and U1_pinv ... U3_pinv are inverted basis matrices
	   the inversion is done with a matrix pseudoinverse computation
	 */
    void derive_core( const t3_type& data_, t3_core_type& core_, const u1_type& U1_, const u2_type& U2_, const u3_type& U3_ );
	//faster: but only if basis matrices are orthogonal
	void derive_core_orthogonal_bases( const t3_type& data_, t3_core_type& core_, const u1_type& U1_, const u2_type& U2_, const u3_type& U3_ );

	/*	higher-order singular value decomposition (HOSVD) with full rank decomposition (also known as Tucker decomposition). 
		see: De Lathauer et al, 2000a: A multilinear singular value decomposition. 
		the hosvd can be computed (a) with n-mode PCA, i.e., an eigenvalue decomposition on the covariance matrix of every mode's matricization, and 
		(b) by performing a 2D SVD on the matricization of every mode. Matrix matricization means that a tensor I1xI2xI3 is unfolded/sliced into one matrix
		with the dimensions I1xI2I3, which corresponds to a matrizitation alonge mode I1.
		other known names for HOSVD: n-mode SVD, 3-mode factor analysis (3MFA, tucker3), 3M-PCA, n-mode PCA, higher-order SVD
	 */
	void hosvd( const t3_type& data_ );
	void hosvd_on_eigs( const t3_type& data_ );
	void hosvd_mode1( const t3_coeff_type& data_, u1_type& U1_ ) const;
	void hosvd_mode2( const t3_coeff_type& data_, u2_type& U2_ ) const;
	void hosvd_mode3( const t3_coeff_type& data_, u3_type& U3_ ) const;

		
	/*	higher-order orthogonal iteration (HOOI) is a truncated HOSVD decompositions, i.e., the HOSVD components are of lower-ranks. An optimal rank-reduction is 
		performed with an alternating least-squares (ALS) algorithm, which minimizes the error between the approximated and orignal tensor based on the Frobenius norm
		see: De Lathauwer et al, 2000b; On the best rank-1 and rank-(RRR) approximation of higher-order tensors.
		the HOOI can be computed based on (a) n-mode PCA, i.e., an eigenvalue decomposition on the covariance matrix of every mode's matriciziation, and 
		(b) by performing a 2D SVD on the matricization of every mode. Matrix matricization means that a tensor I1xI2xI3 is unfolded/sliced into one matrix
		with the dimensions I1xI2I3, which corresponds to a matrizitation alonge mode I1.
	 */
	void hooi( const t3_type& data_ );
		
	void optimize_mode1( const t3_coeff_type& data_, tensor3< I1, R2, R3, T_coeff >& projection_, const u2_type& U2_, const u3_type& U3_ ) const;
	void optimize_mode2( const t3_coeff_type& data_, tensor3< R1, I2, R3, T_coeff >& projection_, const u1_type& U1_, const u3_type& U3_ ) const;		
	void optimize_mode3( const t3_coeff_type& data_, tensor3< R1, R2, I3, T_coeff >& projection_, const u1_type& U1_, const u2_type& U2_ ) const;
	
		
	void tucker_als( const t3_type& data_ );	
		
	template< size_t K1, size_t K2, size_t K3>
	void reduce_ranks( const tucker3_tensor< K1, K2, K3, I1, I2, I3, T_value, T_coeff >& other ); //call TuckerJI.reduce_ranks(TuckerKI) K1 -> R1, K2 -> R2, K3 -> R3

	template< size_t K1, size_t K2, size_t K3>
	void subsampling( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor  );

	template< size_t K1, size_t K2, size_t K3>
	void subsampling_on_average( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor  );

	template< size_t K1, size_t K2, size_t K3>
	void region_of_interest( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, 
							const size_t& start_index1, const size_t& end_index1, 
							const size_t& start_index2, const size_t& end_index2, 
							const size_t& start_index3, const size_t& end_index3);
	
	
private:
    t3_core_type _core ;
	u1_type _u1 ;
	u2_type _u2 ;
	u3_type _u3 ;
	
}; // class tucker3_tensor


#define VMML_TEMPLATE_STRING    	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T_value, typename T_coeff >
#define VMML_TEMPLATE_CLASSNAME     tucker3_tensor< R1, R2, R3, I1, I2, I3, T_value, T_coeff >


VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tucker3_tensor( t3_core_type& core, u1_type& U1, u2_type& U2, u3_type& U3 )
: _core(core), _u1(U1), _u2(U2), _u3(U3)
{	
}
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tucker3_tensor( )
: _core(t3_core_type()), _u1(u1_type()), _u2(u2_type()), _u3(u3_type())
{
}
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruct( t3_type& data_ ) const
{
	t3_coeff_type data;
	data.convert_from_type( data_ );
	data.full_tensor3_matrix_multiplication( _core, _u1, _u2, _u3 );
	data_.convert_from_type( data );
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::decompose( const t3_type& data_ )
{
	tucker_als( data_ );
		
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::tucker_als( const t3_type& data_ )
{
	hooi( data_ );
	
	//derive core
	derive_core_orthogonal_bases( data_, _core, _u1, _u2, _u3 );
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode1( const t3_coeff_type& data_, u1_type& U1_ ) const
{
	mode1_matricization_type u; // -> u1
	data_.lateral_matricization( u);
		
	vector< I2*I3, T_coeff > lambdas;
	lapack_svd< I1, I2*I3, T_coeff > svd;
	if( svd.compute_and_overwrite_input( u, lambdas ))
		u.get_sub_matrix( U1_ );
	else 
		U1_.zero();
}

	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode2( const t3_coeff_type& data_, u2_type& U2_ ) const
{
	mode2_matricization_type u; // -> u2
	data_.frontal_matricization( u);
	
	vector< I1*I3, T_coeff > lambdas;
	lapack_svd< I2, I1*I3, T_coeff > svd;
	if( svd.compute_and_overwrite_input( u, lambdas ))
		u.get_sub_matrix( U2_ );
	else 
		U2_.zero();
}



VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode3( const t3_coeff_type& data_, u3_type& U3_ ) const
{
	mode3_matricization_type u; //-> u3
	data_.horizontal_matricization( u);
	
	vector< I1*I2, T_coeff > lambdas;
	lapack_svd< I3, I1*I2, T_coeff > svd;
	if( svd.compute_and_overwrite_input( u, lambdas ))
		u.get_sub_matrix( U3_ );
	else 
		U3_.zero();
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd( const t3_type& data_ )
{	
	t3_coeff_type data;
	data.convert_from_type( data_ );
	
	hosvd_mode1( data, _u1 );
	hosvd_mode2( data, _u2 );
	hosvd_mode3( data, _u3 );
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_on_eigs( const t3_type& data_ )
{
	//matricization along each mode (backward matricization after Lathauwer et al. 2000a)
	mode1_matricization_type m_lateral; // -> u1
	mode2_matricization_type m_frontal; // -> u2
	mode3_matricization_type m_horizontal; //-> u3
	data_.lateral_matricization( m_lateral);
	data_.frontal_matricization( m_frontal);
	data_.horizontal_matricization( m_horizontal);
	
	//std::cout << "tensor input for tucker, method1: " << std::endl << tensor_ << std::endl;
	
	//2-mode PCA for each matricization A_n: (1) covariance matrix, (2) SVD
	//covariance matrix S_n for each matrizitation A_n
	matrix< I1, I1, T_coeff > s1;
	matrix< I2, I2, T_coeff > s2;
	matrix< I3, I3, T_coeff > s3;
	
	s1.multiply(m_lateral, transpose(m_lateral));
	s2.multiply(m_frontal, transpose(m_frontal));
	s3.multiply(m_horizontal, transpose(m_horizontal));
	
	/*std::cout << "covariance matrix s1: " << std::endl << s1 << std::endl 
	 << "covariance matrix s2: " << s2 << std::endl
	 << "covariance matrix s3: " << s3 << std::endl;*/
	
	//eigenvalue decomposition for each covariance matrix
	
	
}

	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode1( const t3_coeff_type& data_, tensor3< I1, R2, R3, T_coeff >& projection_, const u2_type& U2_, const u3_type& U3_ ) const
{
#if 0
	//compute pseudo inverse for matrices u2,u3
	u2_type u2_pinv_t ;
	u3_type u3_pinv_t ;
	
	compute_pseudoinverse<  u2_type > compute_pinv_u2;
	compute_pinv_u2( U2_, u2_pinv_t );	
	compute_pseudoinverse<  u3_type > compute_pinv_u3;
	compute_pinv_u3( U3_, u3_pinv_t );
	
	u2_inv_type u2_pinv = transpose( u2_pinv_t );
	u3_inv_type u3_pinv = transpose( u3_pinv_t );
#endif
	
	u2_inv_type u2_pinv = transpose( U2_ );
	u3_inv_type u3_pinv = transpose( U3_ );
	
	//backward cyclic matricization (after Lathauwer et al., 2000a)
	tensor3< I1, R2, I3, T_coeff > tmp;
	tmp.multiply_frontal( data_, u2_pinv );
	projection_.multiply_horizontal( tmp, u3_pinv );
}
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode2( const t3_coeff_type& data_, tensor3< R1, I2, R3, T_coeff >& projection_, const u1_type& U1_, const u3_type& U3_ ) const
{
#if 0
	//compute pseudo inverse for matrices u2,u3
	u1_type u1_pinv_t ;
	u3_type u3_pinv_t ;
	
	compute_pseudoinverse<  u1_type > compute_pinv_u1;
	compute_pinv_u1( U1_, u1_pinv_t );	
	compute_pseudoinverse<  u3_type > compute_pinv_u3;
	compute_pinv_u3( U3_, u3_pinv_t );
	
	u1_inv_type u1_pinv = transpose( u1_pinv_t );
	u3_inv_type u3_pinv = transpose( u3_pinv_t );
#endif
	
	u1_inv_type u1_pinv = transpose( U1_ );
	u3_inv_type u3_pinv = transpose( U3_ );
	
	//backward cyclic matricization (after Lathauwer et al., 2000a)
	tensor3< R1, I2, I3, T_coeff > tmp;
	tmp.multiply_lateral( data_, u1_pinv );
	projection_.multiply_horizontal( tmp, u3_pinv );
}	

	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode3( const t3_coeff_type& data_, tensor3< R1, R2, I3, T_coeff >& projection_, const u1_type& U1_, const u2_type& U2_ ) const
{
#if 0
	//compute pseudo inverse for matrices u2,u3
	u1_type u1_pinv_t ;
	u2_type u2_pinv_t ;
	
	compute_pseudoinverse<  u1_type > compute_pinv_u1;
	compute_pinv_u1( U1_, u1_pinv_t );
	compute_pseudoinverse<  u2_type > compute_pinv_u2;
	compute_pinv_u2( U2_, u2_pinv_t );	
	
	u1_inv_type u1_pinv = transpose( u1_pinv_t );
	u2_inv_type u2_pinv = transpose( u2_pinv_t );
#endif
	
	u1_inv_type u1_pinv = transpose( U1_ );
	u2_inv_type u2_pinv = transpose( U2_ );

	//backward cyclic matricization (after Lathauwer et al., 2000a)
	tensor3< R1, I2, I3, T_coeff > tmp;
	tmp.multiply_lateral( data_, u1_pinv );
	projection_.multiply_frontal( tmp, u2_pinv );
}
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hooi( const t3_type& data_ )
{
	//intialize basis matrices
	hosvd( data_ );
	
	t3_coeff_type data;
	data.convert_from_type( data_ );
		
	//compute best rank-(R1, R2, R3) approximation (Lathauwer et al., 2000b)
	t3_type approximated_data;
	reconstruct( approximated_data );
	float_t f_norm = approximated_data.frobenius_norm();
	float_t max_f_norm = data_.frobenius_norm();
	//std::cout << "frobenius norm original: " << max_f_norm << std::endl;
	
	float_t last_f_norm = f_norm;
	float_t improvement = max_f_norm - f_norm;
	float_t min_improvement = 0.1;
	size_t i = 0;
	size_t max_iterations = 3;
	
	while( improvement > min_improvement && i < max_iterations )
	{
		
		//optimize for mode 1
		tensor3< I1, R2, R3, T_coeff > projection1; 
		optimize_mode1( data, projection1, _u2, _u3);
		hosvd_mode1( projection1, _u1 );
		
		//optimize for mode 2
		tensor3< R1, I2, R3, T_coeff > projection2; 
		optimize_mode2( data, projection2, _u1, _u3);
		hosvd_mode2( projection2, _u2 );
		
		//optimize for mode 3
		tensor3< R1, R2, I3, T_coeff > projection3; 
		optimize_mode3( data, projection3, _u1, _u2);
		hosvd_mode3( projection3, _u3);
		
		set_u1( _u1 );
		set_u2( _u2 );
		set_u3( _u3 );
		derive_core_orthogonal_bases(data_, _core, _u1, _u2, _u3);
		set_core( _core );
		
		reconstruct( approximated_data );
		f_norm = approximated_data.frobenius_norm();
		improvement = f_norm - last_f_norm;
		last_f_norm = f_norm;
		
		//std::cout << "iteration '" << i << "': frobenius norm: " << std::setprecision(8) << f_norm << ", improvement: " << improvement << std::endl;
		
		++i;
	}
	
	//std::cout << "number of iterations: " << i << std::endl;
	
	derive_core_orthogonal_bases(data_, _core, _u1, _u2, _u3);
	
}	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::derive_core_orthogonal_bases( const t3_type& data_, t3_core_type& core_, const u1_type& U1_, const u2_type& U2_, const u3_type& U3_ )
{
	matrix< R1, I1, T_coeff > u1_inv = transpose( U1_ );
	matrix< R2, I2, T_coeff > u2_inv = transpose( U2_ );
	matrix< R3, I3, T_coeff > u3_inv = transpose( U3_ );
	
	t3_coeff_type data = data_;
	core_.full_tensor3_matrix_multiplication( data, u1_inv, u2_inv, u3_inv );
}
	
	
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::derive_core( const t3_type& data_, t3_core_type& core_, const u1_type& U1_, const u2_type& U2_, const u3_type& U3_ )
{

	//compute pseudo inverse for matrices u1-u3
	u1_type u1_pinv_t ;
	u2_type u2_pinv_t ;
	u3_type u3_pinv_t ;
	
	
	compute_pseudoinverse<  u1_type > compute_pinv_u1;
	compute_pinv_u1( U1_, u1_pinv_t );
	compute_pseudoinverse<  u2_type > compute_pinv_u2;
	compute_pinv_u2( U2_, u2_pinv_t );	
	compute_pseudoinverse<  u3_type > compute_pinv_u3;
	compute_pinv_u3( U3_, u3_pinv_t );
	
	matrix< R1, I1, T_coeff > u1_pinv = transpose( u1_pinv_t );
	matrix< R2, I2, T_coeff > u2_pinv = transpose( u2_pinv_t );
	matrix< R3, I3, T_coeff > u3_pinv = transpose( u3_pinv_t );
		
	core_.full_tensor3_matrix_multiplication( data_, u1_pinv, u2_pinv, u3_pinv );


	//previous version of compute core	
	//	for( size_t R3 = 0; R3 < R3; ++R3 ) 
	//	{
	//		for( size_t R1 = 0; R1 < R1; ++R1 ) 
	//		{
	//			for( size_t R2 = 0; R2 < R2; ++R2 ) 
	//			{
	//				float_t sum_i1_i2_i3 = 0.0;
	//				for( size_t i3 = 0; i3 < I3; ++i3 ) 
	//				{
	//					for( size_t i1 = 0; i1 < I1; ++i1 ) 
	//					{
	//						for( size_t i2 = 0; i2 < I2; ++i2 ) 
	//						{
	//							sum_i1_i2_i3 += U1_.at( i1, R1 ) * U2_.at( i2, R2 ) * U3_.at( i3, R3 ) * data_.at( i1, i2, i3 );
	//						}
	//					}
	//				}
	//				core_.at( R1, R2, R3 ) = sum_i1_i2_i3;
	//			}
	//		}
	//	}
	//		

}
	
	



VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::reduce_ranks( const tucker3_tensor< K1, K2, K3, I1, I2, I3, T_value, T_coeff >& other )
//TuckerJI.rank_recuction(TuckerKI) K1 -> R1, K2 -> R2, K3 -> R3; I1, I2, I3 stay the same
{
	assert(R1 <= K1);
	assert(R2 <= K2);
	assert(R3 <= K3);	
		
	//reduce basis matrices
	matrix< I1, K1, T_coeff > u1;
	other.get_u1( u1);
	for( size_t r1 = 0; r1 < R1; ++r1 ) 
	{
		_u1.set_column( r1, u1.get_column( r1 ));
	}
	
	matrix< I2, K2, T_coeff > u2;
	other.get_u2( u2 );
	for( size_t r2 = 0; r2 < R2; ++r2) 
	{
		_u2.set_column( r2, u2.get_column( r2 ));
	}
	
	matrix< I3, K3, T_coeff > u3;
	other.get_u3( u3 );
	for( size_t r3 = 0; r3 < R3; ++r3) 
	{
		_u3.set_column( r3, u3.get_column( r3 ));
	}
	
	//reduce core
	tensor3<K1, K2, K3, T_coeff > other_core;
	other.get_core( other_core );

	for( size_t r3 = 0; r3 < R3; ++r3 ) 
	{
		for( size_t r1 = 0; r1 < R1; ++r1 ) 
		{
			for( size_t r2 = 0; r2 < R2; ++r2 ) 
			{
				_core.at( r1, r2, r3 ) = other_core.at( r1, r2, r3 );
			}
		}
	}

}



VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::subsampling( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor )
{
	assert(I1 <= K1);
	assert(I1 <= K2);
	assert(I1 <= K3);	
	
	//subsample basis matrices
	matrix< K1, R1, T_coeff > u1;
	other.get_u1( u1 );
	for( size_t i1 = 0, i = 0; i1 < K1; i1 += factor, ++i ) 
	{
		_u1.set_row( i, u1.get_row( i1 ));
	}
	
	matrix< K2, R2, T_coeff > u2;
	other.get_u2( u2 );
	for( size_t i2 = 0,  i = 0; i2 < K2; i2 += factor, ++i) 
	{
		_u2.set_row( i, u2.get_row( i2 ));
	}
	
	matrix< K3, R3, T_coeff > u3 ;
	other.get_u3( u3);
	for( size_t i3 = 0,  i = 0; i3 < K3; i3 += factor, ++i) 
	{
		_u3.set_row( i, u3.get_row( i3 ));
	}
	
	other.get_core(_core );
}


VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::subsampling_on_average( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor )
{
	assert(I1 <= K1);
	assert(I1 <= K2);
	assert(I1 <= K3);	
	
	
	//subsample basis matrices
	matrix< K1, R1, T_coeff > u1;
	other.get_u1( u1 );
	for( size_t i1 = 0, i = 0; i1 < K1; i1 += factor, ++i )
	{
		vector< R1, T_coeff > tmp_row = u1.get_row( i1 );
		T_coeff num_items_averaged = 1;
		for( size_t j = i1+1; (j < (factor+i1)) & (j < K1); ++j, ++num_items_averaged )
			tmp_row += u1.get_row( j );

		tmp_row /= num_items_averaged;
		_u1.set_row( i, tmp_row);
	}
	
	matrix< K2, R2, T_coeff > u2;
	other.get_u2( u2 );
	for( size_t i2 = 0,  i = 0; i2 < K2; i2 += factor, ++i) 
	{
		vector< R2, T_coeff > tmp_row = u2.get_row( i2 );
		T_coeff num_items_averaged = 1;
		for( size_t j = i2+1; (j < (factor+i2)) & (j < K2); ++j, ++num_items_averaged )
			tmp_row += u2.get_row( j );

		tmp_row /= num_items_averaged;
		_u2.set_row( i, u2.get_row( i2 ));
	}
	
	matrix< K3, R3, T_coeff > u3;
	other.get_u3( u3);
	for( size_t i3 = 0,  i = 0; i3 < K3; i3 += factor, ++i) 
	{
		vector< R3, T_coeff > tmp_row = u3.get_row( i3 );
		T_coeff num_items_averaged = 1;
		for( size_t j = i3+1; (j < (factor+i3)) & (j < K3); ++j, ++num_items_averaged )
			tmp_row += u3.get_row( j );
		
		tmp_row /= num_items_averaged;
		_u3.set_row( i, u3.get_row( i3 ));
	}
  	
	 other.get_core( _core );
}




VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::region_of_interest( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, 
											const size_t& start_index1, const size_t& end_index1, 
											const size_t& start_index2, const size_t& end_index2, 
											const size_t& start_index3, const size_t& end_index3)
{
	assert(I1 <= K1);
	assert(I1 <= K2);
	assert(I1 <= K3);
	assert(start_index1 < end_index1);
	assert(start_index2 < end_index2);
	assert(start_index3 < end_index3);
	assert(end_index1 < K1);
	assert(end_index2 < K2);
	assert(end_index3 < K3);
	
	//region_of_interes of basis matrices
	matrix< K1, R1, T_coeff > u1;
	other.get_u1( u1 );
	for( size_t i1 = start_index1,  i = 0; i1 < end_index1; ++i1, ++i ) 
	{
		_u1.set_row( i, u1.get_row( i1 ));
	}
	
	matrix< K2, R2, T_coeff> u2;
	other.get_u2( u2 );
	for( size_t i2 = start_index2,  i = 0; i2 < end_index2; ++i2, ++i) 
	{
		_u2.set_row( i, u2.get_row( i2 ));
	}
	
	matrix< K3, R3, T_coeff > u3;
	other.get_u3( u3 );
	for( size_t i3 = start_index3,  i = 0; i3 < end_index3; ++i3, ++i) 
	{
		_u3.set_row( i, u3.get_row( i3 ));
	}
	
	other.get_core( _core );
	
}
	
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::export_to( std::vector< T_coeff >& data_ ) const
{
	u1_const_iterator  it = _u1.begin(),
	it_end = _u1.end();
	for( ; it != it_end; ++it )
	{
		data_.push_back( *it );
	}
	
	u2_const_iterator  u2_it = _u2.begin(),
	u2_it_end = _u2.end();
	for( ; u2_it != u2_it_end; ++u2_it )
	{
		data_.push_back( *u2_it );
	}

	u3_const_iterator  u3_it = _u3.begin(),
	u3_it_end = _u3.end();
	for( ; u3_it != u3_it_end; ++u3_it )
	{
		data_.push_back( *u3_it );
	}
	
	t3_core_const_iterator  it_core = _core.begin(),
	it_core_end = _core.end();
	for( ; it_core != it_core_end; ++it_core )
	{
		data_.push_back( *it_core );
	}
}
	
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_from( const std::vector< T_coeff >& data_ )
{
	size_t i = 0; //iterator over data_
	size_t data_size = (size_t) data_.size();
	
    if ( data_size != SIZE  )
        VMMLIB_ERROR( "import_from: the input data must have the size R1xR2xR3 + R1xI1 + R2xI2 + R3xI3 ", VMMLIB_HERE );
	
	u1_iterator  it = _u1.begin(),
	it_end = _u1.end();
	for( ; it != it_end; ++it, ++i )
	{
		*it = data_.at(i);
	}
	
	u2_iterator  u2_it = _u2.begin(),
	u2_it_end = _u2.end();
	for( ; u2_it != u2_it_end; ++u2_it, ++i )
	{
		*u2_it = data_.at(i);
	}
	
	u3_iterator  u3_it = _u3.begin(),
	u3_it_end = _u3.end();
	for( ; u3_it != u3_it_end; ++u3_it, ++i )
	{
		*u3_it = data_.at(i);
	}
	
	t3_core_iterator  it_core = _core.begin(),
	it_core_end = _core.end();
	for( ; it_core != it_core_end; ++it_core, ++i )
	{
		*it_core = data_.at(i);
	}
	
}	
	


#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
} // namespace vmml

#endif
