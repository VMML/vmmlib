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

		
		//higher-order power method (lathauwer et al., 2000b)
		static void als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_, const size_t max_iterations_ = 100 );
		static void reconstruct( t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, const lambda_type& lambdas_ );
		
	protected:
		
		static void optimize_mode1( const t3_type& data_, u1_type& u1, const u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ );
		static void optimize_mode2( const t3_type& data_, const u1_type& u1_, u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ );		
		static void optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ );
		
		
	};
	


#define VMML_TEMPLATE_STRING    	template< size_t R, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME     t3_hopm< R, I1, I2, I3, T >


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_, const size_t max_iterations_ )
{
	t3_type approximated_data;
	t3_type residual_data;
	residual_data.zero();
	
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
	//hosvd_mode1( data_, _u1 ); inital guess not needed for u1 since it will be computed in the first optimization step
	t3_hosvd< R, R, R, I1, I2, I3, T >::apply_mode2( data_, u2_ );
	t3_hosvd< R, R, R, I1, I2, I3, T >::apply_mode3( data_, u3_ );
	
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
		reconstruct( approximated_data, u1_, u2_, u3_, lambdas_ );
		approx_norm = approximated_data.frobenius_norm();
		residual_data = data_ - approximated_data;
		normresidual = residual_data.frobenius_norm();
		fit = 1 - ( normresidual / max_f_norm ); 
		fitchange = fabs(fitold - fit);
		
#if CP_LOG
		std::cout << "iteration '" << i << "', fit: " << fit 
		<< ", fitdelta: " << fitchange 
		<< ", frobenius norm: " << approx_norm << std::endl;		
#endif
		++i;
	} // end ALS
	
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode1( const t3_type& data_, u1_type& u1_, const u2_type& u2_, const u3_type& u3_, lambda_type& lambdas_ )
{	
	mode1_matricization_type* unfolding = new mode1_matricization_type; // -> u1
	//data_.horizontal_unfolding_bwd( *unfolding ); //lathauwer
	data_.frontal_unfolding_fwd( *unfolding ); //lathauwer
	
	//std::cout << "u2:\n" << u2_ << std::endl;

	typedef matrix< I2*I3, R, T > krp_matrix_type;
	krp_matrix_type* u1_krp  = new krp_matrix_type;
	*u1_krp =u3_.khatri_rao_product( u2_ );	
	u1_type* u_new = new u1_type;
	u_new->multiply( *unfolding, *u1_krp );
	
	m_r2_type* u2_r = new m_r2_type;
	m_r2_type* u3_r = new m_r2_type;
	
	u2_r->multiply( transpose( u2_ ), u2_ );
	u3_r->multiply( transpose( u3_ ), u3_ );
	u2_r->multiply_piecewise( *u3_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *u2_r, *pinv_t );
	
	u1_.multiply( *u_new, transpose(*pinv_t) );
	
	*u_new = u1_;
	u_new->multiply_piecewise( *u_new ); //2 norm
	u_new->columnwise_sum( lambdas_ );
	lambdas_.sqrt_elementwise();
	lambda_type* tmp = new lambda_type;
	*tmp = lambdas_;
	tmp->reciprocal();
	m_r2_type* diag_lambdas = new m_r2_type;
	diag_lambdas->diag( *tmp );
	
	*u_new = u1_;
	u1_.multiply( *u_new, *diag_lambdas ); 
	
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
	u_new->multiply( *unfolding, *u2_krp );
	
	m_r2_type* u1_r = new m_r2_type;
	m_r2_type* u3_r = new m_r2_type;
	
	u1_r->multiply( transpose( u1_ ), u1_ );
	u3_r->multiply( transpose( u3_ ), u3_ );
	u1_r->multiply_piecewise( *u3_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *u1_r, *pinv_t );
	
	u2_.multiply( *u_new, transpose(*pinv_t) );
	
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
	
	*u_new = u2_;
	u2_.multiply( *u_new, *diag_lambdas );
	
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
	u_new->multiply( *unfolding, *u3_krp );
	
	typedef matrix< R, R , T > m_r2_type;
	m_r2_type* u1_r = new m_r2_type;
	m_r2_type* u2_r = new m_r2_type;
	
	u1_r->multiply( transpose( u1_ ), u1_ );
	u2_r->multiply( transpose( u2_ ), u2_ );
	u1_r->multiply_piecewise( *u2_r );
	
	m_r2_type* pinv_t = new m_r2_type;
	compute_pseudoinverse< m_r2_type > compute_pinv;
	compute_pinv( *u1_r, *pinv_t );
	
	u3_.multiply( *u_new, transpose(*pinv_t) );
	
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
	
	*u_new = u3_;
	u3_.multiply( *u_new, *diag_lambdas );
	
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

	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
}//end vmml namespace

#endif
