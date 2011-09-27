/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * the higher-order power method (HOPM) is also known as CP-ALS (ALS: alternating least squares)
 * CP stands for Candecomp/Parafac (1970)
 * references:
 * - Carroll & Chang, 1970: Analysis of Individual Differences in Multidimensional Scaling via an N-way generalization of ``Eckart--Young'' decompositions, Psychometrika.
 * - Harshman, 1970: Foundations of the PARAFAC procedure: Models and conditions for an 'explanatory' multi-modal factor analysis, UCLA Working Papers in Phonetics.
 * - De Lathauwer, De Moor, Vandewalle, 2000: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * 
 */

#ifndef __VMML__T3_HOPM__HPP__
#define __VMML__T3_HOPM__HPP__

#include <vmmlib/t3_hosvd.hpp>
#include <vmmlib/matrix_pseudoinverse.hpp>
#include <vmmlib/blas_dgemm.hpp>

namespace vmml
{
	
	template< size_t R, size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_hopm
	{
	public:    
		typedef tensor3< I1, I2, I3, T > t3_type;
		
		typedef vector< R, T > lambda_type;

		typedef matrix< I1, R, T > u1_type;
		typedef matrix< I2, R, T > u2_type;
		typedef matrix< I3, R, T > u3_type;
		
		typedef matrix< R, I1, T > u1_inv_type;
		typedef matrix< R, I2, T > u2_inv_type;
		typedef matrix< R, I3, T > u3_inv_type;

		typedef matrix< I1, I2*I3, T > u1_unfolded_type;
		typedef matrix< I2, I1*I3, T > u2_unfolded_type;
		typedef matrix< I3, I1*I2, T > u3_unfolded_type;

		typedef matrix< R, R , T > m_r2_type;

		typedef typename lambda_type::iterator lvalue_iterator;
		typedef typename lambda_type::const_iterator lvalue_const_iterator;
		typedef std::pair< T, size_t >  lambda_pair_type;
		
		//higher-order power method (lathauwer et al., 2000b)
		template< typename T_init >
		static void als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_, T_init init, const size_t max_iterations_ = 100 );
		static void reconstruct( t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, const lambda_type& lambdas_ );
		
		// init functors 
		struct init_hosvd
		{
			inline void operator()( const t3_type& data_, u2_type& u2_, u3_type& u3_ )
			{
				t3_hosvd< R, R, R, I1, I2, I3, T >::apply_mode2( data_, u2_ );
				t3_hosvd< R, R, R, I1, I2, I3, T >::apply_mode3( data_, u3_ );
			}
		};
		
		struct init_random
		{
			inline void operator()( const t3_type& data_, u2_type& u2_, u3_type& u3_ )
			{
				int seed = time( NULL );
				u2_.set_random( seed );
				u3_.set_random( rand() );
			}
		};

		struct init_dct
		{
			inline void operator()( const t3_type& data_, u2_type& u2_, u3_type& u3_ )
			{
				u2_.set_dct();
				u3_.set_dct();
			}
		};
		
		
	protected:
		
		static void optimize_mode1( const t3_type& data_, u1_type& u1, const u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ );
		static void optimize_mode2( const t3_type& data_, const u1_type& u1_, u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ );		
		static void optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ );
		
		template< size_t J, size_t K, size_t L >
		static void optimize( const matrix< J, K*L, T >& unfolding_, 
							 matrix< J, R, T >& uj_, 
							 const matrix< K, R, T >& uk_, const matrix< L, R, T >& ul_,
							 vector< R, T>& lambdas_ 
							 );	
		
		static void sort_dec( u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ );
		
		// comparison functor 
		struct lambda_compare
		{
			inline bool operator()( const lambda_pair_type& a, const lambda_pair_type& b )
			{
				return fabs( a.first ) > fabs( b.first );
			}
		};
	
	};
	


#define VMML_TEMPLATE_STRING        template< size_t R, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME     t3_hopm< R, I1, I2, I3, T >


	
VMML_TEMPLATE_STRING
template< typename T_init>
void 
VMML_TEMPLATE_CLASSNAME::als( const t3_type& data_, 
							 u1_type& u1_, u2_type& u2_, u3_type& u3_, 
							 lambda_type& lambdas_, 
							 T_init init,
							 const size_t max_iterations_ )
{
	t3_type* approximated_data = new t3_type;
	t3_type* residual_data = new t3_type;
	residual_data->zero();
	
	double approx_norm = 0;
	double max_f_norm = data_.frobenius_norm();
	double normresidual  = 0;
	double fit = 0;
	if (max_f_norm == 0 )
		fit = 1;
	double fitchange = 1;
	double fitold = fit;
	double fitchange_tolerance = 1.0e-4;
	
	//intialize u1-u3
	//inital guess not needed for u1 since it will be computed in the first optimization step
	init( data_, u2_, u3_ );

#if CP_LOG
	std::cout << "CP ALS: HOPM (for tensor3) " << std::endl;
#endif	
	
	size_t i = 0;
	//size_t max_iterations = 100;
	while( (fitchange >= fitchange_tolerance) && ( i < max_iterations_ ) ) //do until converges
	{
		fitold = fit;
		optimize_mode1( data_, u1_, u2_, u3_, lambdas_ );
		optimize_mode2( data_, u1_, u2_, u3_, lambdas_ );
		optimize_mode3( data_, u1_, u2_, u3_, lambdas_ );
		
		//Reconstruct cptensor and measure norm of approximation
		reconstruct( *approximated_data, u1_, u2_, u3_, lambdas_ );
		approx_norm = approximated_data->frobenius_norm();
		*residual_data = data_ - *approximated_data;
		normresidual = residual_data->frobenius_norm();
		fit = 1 - ( normresidual / max_f_norm ); 
		fitchange = fabs(fitold - fit);
		
#if CP_LOG
		std::cout << "iteration '" << i 
		<< "', fit: " << fit 
		<< ", fitdelta: " << fitchange 
		<< ", frobenius norm: " << approx_norm << std::endl;		
#endif
		++i;
	} // end ALS

	//sort lambdas by decreasing magnitude
	sort_dec( u1_, u2_, u3_, lambdas_ );
	
	delete residual_data;
	delete approximated_data;
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode1( const t3_type& data_, u1_type& u1_, const u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ )
{	
	u1_unfolded_type* unfolding = new u1_unfolded_type; // -> u1
	//data_.horizontal_unfolding_bwd( *unfolding ); //lathauwer
	data_.frontal_unfolding_fwd( *unfolding ); 
	
	optimize( *unfolding, u1_, u2_, u3_, lambdas_ );
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode2( const t3_type& data_, const u1_type& u1_, u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_  )
{
	u2_unfolded_type* unfolding = new u2_unfolded_type; // -> u2
	data_.frontal_unfolding_bwd( *unfolding ); //lathauwer
	//data_.horizontal_unfolding_fwd( *unfolding );
	
	optimize( *unfolding, u2_, u1_, u3_, lambdas_ );
}	


VMML_TEMPLATE_STRING
void  
VMML_TEMPLATE_CLASSNAME::optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ )
{
	u3_unfolded_type* unfolding = new u3_unfolded_type; //-> u3
	//data_.horizontal_unfolding_bwd( *unfolding );//lathauwer
	data_.lateral_unfolding_fwd( *unfolding );
	
	optimize( *unfolding, u3_, u1_, u2_, lambdas_ );
}
	
VMML_TEMPLATE_STRING
template< size_t J, size_t K, size_t L >
void
VMML_TEMPLATE_CLASSNAME::optimize( 
								  const matrix< J, K*L, T >& unfolding_, 
								  matrix< J, R, T >& uj_, 
								  const matrix< K, R, T >& uk_, const matrix< L, R, T >& ul_,
								  vector< R, T>& lambdas_ 
								  )	
{
	typedef matrix< K*L, R, T > krp_matrix_type;
	krp_matrix_type* krp_prod  = new krp_matrix_type;
	ul_.khatri_rao_product( uk_, *krp_prod );	
	matrix< J, R, T >* u_new = new matrix< J, R, T >;
	
	blas_dgemm< J, K*L, R, T>* blas_dgemm1 = new blas_dgemm< J, K*L, R, T>;
	blas_dgemm1->compute( unfolding_, *krp_prod, *u_new );
	delete blas_dgemm1;	
	
	//square matrix of U_l and U_k
	m_r2_type* uk_r = new m_r2_type;
	m_r2_type* ul_r = new m_r2_type;
	blas_dgemm< R, K, R, T> blas_dgemm2;
	blas_dgemm2.compute_t( uk_, *uk_r );
	blas_dgemm< R, L, R, T> blas_dgemm3;
	blas_dgemm3.compute_t( ul_, *ul_r );
	
	uk_r->multiply_piecewise( *ul_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *uk_r, *pinv_t );
	
	blas_dgemm< J, R, R, T> blas_dgemm4;
	blas_dgemm4.compute_bt( *u_new, *pinv_t, uj_ );
	
	*u_new = uj_;
	u_new->multiply_piecewise( *u_new ); //2 norm
	u_new->columnwise_sum( lambdas_ );
	lambdas_.sqrt_elementwise();
	lambda_type* tmp = new lambda_type;
	*tmp = lambdas_;
	tmp->reciprocal();
	m_r2_type* diag_lambdas = new m_r2_type;
	diag_lambdas->diag( *tmp );
	
	blas_dgemm4.compute( uj_, *diag_lambdas, uj_ );
	
	delete krp_prod;
	delete uk_r;
	delete ul_r;
	delete pinv_t;
	delete u_new;
	delete diag_lambdas;
	delete tmp;
}
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruct( t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, const lambda_type& lambdas_ ) 
{
	u1_inv_type* u1_t = new u1_inv_type;
	u2_inv_type* u2_t = new u2_inv_type;
	u3_inv_type* u3_t = new u3_inv_type;
	typedef matrix<  R, I2 * I3, T > m_temp_type;
	m_temp_type* temp =  new m_temp_type; 
	
	*u1_t = transpose( u1_ );
	*u2_t = transpose( u2_ );
	*u3_t = transpose( u3_ );
	
	data_.reconstruct_CP( lambdas_, *u1_t, *u2_t, *u3_t, *temp );
	
	delete temp;
	delete u1_t;
	delete u2_t;
	delete u3_t;
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::sort_dec( u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ ) 
{
	//keep copy of original matrices
	u1_type *orig_u1 = new u1_type( u1_ );
	u2_type *orig_u2 = new u2_type( u2_ );
	u3_type *orig_u3 = new u3_type( u3_ );
	lambda_type sorted_lvalues;
	
	//(1) store permutations of the lambdas (According to larges magnitude
	std::vector< lambda_pair_type > lambda_permut;
	
	lvalue_const_iterator it = lambdas_.begin(), it_end = lambdas_.end();
	size_t counter = 0;
	for( ; it != it_end; ++it, ++counter )
	{
		 lambda_permut.push_back(  lambda_pair_type ( *it, counter ) );
	}
    
    std::sort(
			   lambda_permut.begin(),
			   lambda_permut.end(), 
			   lambda_compare()
			  );
	
	//(2) sort the matrix vectors according to lambda permutations and set sorted lambdas
	typename std::vector< lambda_pair_type >::const_iterator it2 = lambda_permut.begin(), it2_end = lambda_permut.end();
	lvalue_iterator lvalues_it = lambdas_.begin();
	for( counter = 0; it2 != it2_end; ++it2, ++counter, ++lvalues_it )
	{
		*lvalues_it = it2->first;
		u1_.set_column( counter, orig_u1->get_column( it2->second ));
		u2_.set_column( counter, orig_u2->get_column( it2->second ));
		u3_.set_column( counter, orig_u3->get_column( it2->second ));
	}	
		
	delete orig_u1;
	delete orig_u2;
	delete orig_u3;
}	
	
	
	
	
	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
}//end vmml namespace

#endif
