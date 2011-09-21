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

enum init_mode {
	init_hosvd_e,
	init_rand_e,
	init_dct_e
}; 

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

		typedef matrix< I1, I2*I3, T > mode1_matricization_type;
		typedef matrix< I2, I1*I3, T > mode2_matricization_type;
		typedef matrix< I3, I1*I2, T > mode3_matricization_type;	

		typedef matrix< R, R , T > m_r2_type;

		typedef typename lambda_type::iterator lvalue_iterator;
		typedef typename lambda_type::const_iterator lvalue_const_iterator;
		typedef std::pair< T, size_t >  lambda_pair_type;
		
		//higher-order power method (lathauwer et al., 2000b)
		static void als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_, init_mode init_mode_, const size_t max_iterations_ = 100 );
		static void reconstruct( t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, const lambda_type& lambdas_ );

		// only the factor matrices u2, u3 are initialized since we get u1 from the first optimization iteration
		static void init( const t3_type& data_, u2_type& u2_, u3_type& u3_, init_mode init_mode_ );
		static void init_hosvd( const t3_type& data_, u2_type& u2_, u3_type& u3_ );
		static void init_random( const t3_type& data_, u2_type& u2_, u3_type& u3_ );
		static void init_dct( const t3_type& data_, u2_type& u2_, u3_type& u3_ );

		
	protected:
				
		template< size_t M, size_t N >
        static void fill_random_2d( int seed, matrix< M, N, T >& u );
		
		template< size_t M, size_t N >
        static void fill_dct_matrix( matrix< M, N, T >& u );
		
		static void optimize_mode1( const t3_type& data_, u1_type& u1, const u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ );
		static void optimize_mode2( const t3_type& data_, const u1_type& u1_, u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ );		
		static void optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ );
		
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
	


#define VMML_TEMPLATE_STRING    	template< size_t R, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME     t3_hopm< R, I1, I2, I3, T >


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::init( const t3_type& data_, u2_type& u2_, u3_type& u3_, init_mode init_mode_  )
{	
	switch ( init_mode_ )
	{
		case 0:
			init_hosvd( data_, u2_, u3_ );
            break;
		case 1:
			init_random( data_, u2_, u3_ );
            break;
		case 2:
			init_dct( data_, u2_, u3_ );
           break;
		default:
			init_hosvd( data_, u2_, u3_ );
	}
}	
	


VMML_TEMPLATE_STRING
template< size_t M, size_t N >
void 
VMML_TEMPLATE_CLASSNAME::fill_dct_matrix( matrix< M, N, T >& u )
{
	double weight = 0.0f;
	double	fill_value = 0.0f;
	for( size_t row = 0; row < M; ++row )
	{
		//todo: weight = (row > 1) ? (1/M) : (2/M);
		weight = 1 / static_cast< double >(M);
		for( size_t col = 0; col < N; ++col )
		{
			fill_value = (2 * col + 1) * row * (M_PI / (2*M));
			fill_value = cos( fill_value );
			fill_value *= weight * fill_value;
			u.at( row, col ) = static_cast< T >( fill_value )  ;
		}
	}
}		
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::init_dct( const t3_type& data_, u2_type& u2_, u3_type& u3_ )
{	
	fill_dct_matrix( u2_ );
	fill_dct_matrix( u3_ );
}		

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::init_hosvd( const t3_type& data_, u2_type& u2_, u3_type& u3_ )
{	
	t3_hosvd< R, R, R, I1, I2, I3, T >::apply_mode2( data_, u2_ );
	t3_hosvd< R, R, R, I1, I2, I3, T >::apply_mode3( data_, u3_ );
}	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::init_random( const t3_type& data_, u2_type& u2_, u3_type& u3_ )
{	
	int seed = time(NULL);
	fill_random_2d(seed, u2_ );
	fill_random_2d(rand(), u3_ );
}	
	
VMML_TEMPLATE_STRING
template< size_t M, size_t N >
void 
VMML_TEMPLATE_CLASSNAME::fill_random_2d( int seed, matrix< M, N, T >& u)
{
	double fillValue = 0.0f;
	srand(seed);
	for( size_t row = 0; row < M; ++row )
	{
		for( size_t col = 0; col < N; ++col )
		{
			fillValue = rand();
			fillValue /= RAND_MAX;
			u.at( row, col ) = -1.0 + 2.0 * static_cast< double >( fillValue )  ;
		}
	}
}	
	
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::als( const t3_type& data_, 
							 u1_type& u1_, u2_type& u2_, u3_type& u3_, 
							 lambda_type& lambdas_, 
							 init_mode init_mode_,
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
	init( data_, u2_, u3_, init_mode_ );

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
	mode1_matricization_type* unfolding = new mode1_matricization_type; // -> u1
	//data_.horizontal_unfolding_bwd( *unfolding ); //lathauwer
	data_.frontal_unfolding_fwd( *unfolding ); 
	
	typedef matrix< I2*I3, R, T > krp_matrix_type;
	krp_matrix_type* u1_krp  = new krp_matrix_type;
	*u1_krp =u3_.khatri_rao_product( u2_ );	
	u1_type* u_new = new u1_type;
	
	blas_dgemm< I1, I2*I3, R, T>* blas_dgemm1 = new blas_dgemm< I1, I2*I3, R, T>;
	blas_dgemm1->compute( *unfolding, *u1_krp, *u_new );
	delete blas_dgemm1;	
	
	//square matrix of u2 and u3
	m_r2_type* u2_r = new m_r2_type;
	m_r2_type* u3_r = new m_r2_type;
	blas_dgemm< R, I2, R, T> blas_dgemm2;
	blas_dgemm2.compute_t( u2_, *u2_r );
	blas_dgemm< R, I3, R, T> blas_dgemm3;
	blas_dgemm3.compute_t( u3_, *u3_r );
	
	u2_r->multiply_piecewise( *u3_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *u2_r, *pinv_t );
	
	blas_dgemm< I1, R, R, T> blas_dgemm4;
	blas_dgemm4.compute_bt( *u_new, *pinv_t, u1_ );
	
	*u_new = u1_;
	u_new->multiply_piecewise( *u_new ); //2 norm
	u_new->columnwise_sum( lambdas_ );
	lambdas_.sqrt_elementwise();
	lambda_type* tmp = new lambda_type;
	*tmp = lambdas_;
	tmp->reciprocal();
	m_r2_type* diag_lambdas = new m_r2_type;
	diag_lambdas->diag( *tmp );
	
	blas_dgemm4.compute( u1_, *diag_lambdas, u1_ );
	
	delete unfolding;
	delete u1_krp;
	delete u2_r;
	delete u3_r;
	delete pinv_t;
	delete u_new;
	delete diag_lambdas;
	delete tmp;
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode2( const t3_type& data_, const u1_type& u1_, u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_  )
{
	mode2_matricization_type* unfolding = new mode2_matricization_type; // -> u2
	//data_.frontal_unfolding_bwd( *unfolding ); //lathauwer
	data_.frontal_unfolding_bwd( *unfolding );
	
	typedef matrix< I1*I3, R, T > krp_matrix_type;
	krp_matrix_type* u2_krp  = new krp_matrix_type;
	*u2_krp = u3_.khatri_rao_product( u1_ );	
	u2_type* u_new = new u2_type;
	
	blas_dgemm< I2, I1*I3, R, T>* blas_dgemm1 = new blas_dgemm< I2, I1*I3, R, T>;
	blas_dgemm1->compute( *unfolding, *u2_krp, *u_new );
	delete blas_dgemm1;	
	
	//square matrix of u1 and u3
	m_r2_type* u1_r = new m_r2_type;
	m_r2_type* u3_r = new m_r2_type;
	blas_dgemm< R, I1, R, T> blas_dgemm2;
	blas_dgemm2.compute_t( u1_, *u1_r );
	blas_dgemm< R, I3, R, T> blas_dgemm3;
	blas_dgemm3.compute_t( u3_, *u3_r );

	u1_r->multiply_piecewise( *u3_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *u1_r, *pinv_t );
	
	blas_dgemm< I2, R, R, T> blas_dgemm4;
	blas_dgemm4.compute_bt( *u_new, *pinv_t, u2_ );
	
	//normalize with lambdas
	*u_new = u2_;
	u_new->multiply_piecewise( *u_new ); //2 norm
	u_new->columnwise_sum( lambdas_ );
	lambdas_.sqrt_elementwise();
	lambda_type* tmp = new lambda_type;
	*tmp = lambdas_;
	tmp->reciprocal();
	m_r2_type* diag_lambdas = new m_r2_type;
	diag_lambdas->diag( *tmp );
	
	blas_dgemm4.compute( u2_, *diag_lambdas, u2_ );
	
	delete unfolding;
	delete u2_krp;
	delete u1_r;
	delete u3_r;
	delete pinv_t;
	delete u_new;
	delete diag_lambdas;
	delete tmp;
}	


VMML_TEMPLATE_STRING
void  
VMML_TEMPLATE_CLASSNAME::optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ )
{
	mode3_matricization_type* unfolding = new mode3_matricization_type; //-> u3
	//data_.horizontal_unfolding_bwd( *unfolding );//lathauwer
	data_.lateral_unfolding_fwd( *unfolding );
	
	typedef matrix< I1*I2, R, T > krp_matrix_type;
	krp_matrix_type* u3_krp  = new krp_matrix_type;
	*u3_krp = u2_.khatri_rao_product( u1_ );	
	u3_type* u_new = new u3_type;
	
	blas_dgemm< I3, I1*I2, R, T>* blas_dgemm1 = new blas_dgemm< I3, I1*I2, R, T>;
	blas_dgemm1->compute( *unfolding, *u3_krp, *u_new );
	delete blas_dgemm1;	
	
	//square matrix of u1 and u3
	m_r2_type* u1_r = new m_r2_type;
	m_r2_type* u2_r = new m_r2_type;
	blas_dgemm< R, I1, R, T> blas_dgemm2;
	blas_dgemm2.compute_t( u1_, *u1_r );
	blas_dgemm< R, I2, R, T> blas_dgemm3;
	blas_dgemm3.compute_t( u2_, *u2_r );
	
	u1_r->multiply_piecewise( *u2_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *u1_r, *pinv_t );
	
	blas_dgemm< I3, R, R, T> blas_dgemm4;
	blas_dgemm4.compute_bt( *u_new, *pinv_t, u3_ );
	
	//normalize with lambdas
	*u_new = u3_;
	u_new->multiply_piecewise( *u_new ); //2 norm
	u_new->columnwise_sum( lambdas_ );
	lambdas_.sqrt_elementwise();
	lambda_type* tmp = new lambda_type;
	*tmp = lambdas_;
	tmp->reciprocal();
	m_r2_type* diag_lambdas = new m_r2_type;
	diag_lambdas->diag( *tmp );
	
	blas_dgemm4.compute( u3_, *diag_lambdas, u3_ );
	
	delete unfolding;
	delete u3_krp;
	delete u1_r;
	delete u2_r;
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
