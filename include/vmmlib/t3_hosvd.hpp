/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * The tucker3 tensor class is consists of the same components (core tensor, basis matrices u1-u3) as the tucker3 model described in:
 * - Tucker, 1966: Some mathematical notes on three-mode factor analysis, Psychometrika.
 * - De Lathauwer, De Moor, Vandewalle, 2000a: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * 
 */

#ifndef __VMML__T3_HOSVD__HPP__
#define __VMML__T3_HOSVD__HPP__

#include <vmmlib/tensor3.hpp>
#include <vmmlib/lapack_svd.hpp>
#include <vmmlib/lapack_sym_eigs.hpp>
#include <vmmlib/blas_dgemm.hpp>

namespace vmml
{
	
	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_hosvd
	{
	public:    
		
		typedef double T_svd;
		
		typedef tensor3< I1, I2, I3, T > t3_type;
		
		typedef matrix< I1, R1, T > u1_type;
		typedef matrix< I2, R2, T > u2_type;
		typedef matrix< I3, R3, T > u3_type;
		
		typedef matrix< I1, I2*I3, T > u1_unfolded_type;
		typedef matrix< I2, I1*I3, T > u2_unfolded_type;
		typedef matrix< I3, I1*I2, T > u3_unfolded_type;
		
		typedef matrix< I1, I1, T > u1_cov_type;
		typedef matrix< I2, I2, T > u2_cov_type;
		typedef matrix< I3, I3, T > u3_cov_type;
		
		/*	higher-order singular value decomposition (HOSVD) with full rank decomposition (also known as Tucker decomposition). 
		 see: De Lathauer et al, 2000a: A multilinear singular value decomposition. 
		 the hosvd can be computed (a) with n-mode PCA, i.e., an eigenvalue decomposition on the covariance matrix of every mode's matricization, and 
		 (b) by performing a 2D SVD on the matricization of every mode. Matrix matricization means that a tensor I1xI2xI3 is unfolded/sliced into one matrix
		 with the dimensions I1xI2I3, which corresponds to a matrizitation alonge mode I1.
		 other known names for HOSVD: n-mode SVD, 3-mode factor analysis (3MFA, tucker3), 3M-PCA, n-mode PCA, higher-order SVD
		 */

		//FIXME: hosvd_mode1 -> fixme names, e.g., apply_mode1 -> choose either svd_mode1, eigs_mode1
		
		static void apply_mode1( const t3_type& data_, u1_type& u1_ );
		static void apply_mode2( const t3_type& data_, u2_type& u2_ );
		static void apply_mode3( const t3_type& data_, u3_type& u3_ );
		
		static void apply_all( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ );
		static void hosvd( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ );
		static void hoeigs( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ );
		
	protected:
		
		//hosvd
        template< size_t M, size_t N, size_t R >
        static void get_svd_u_red( const matrix< M, N, T >& data_, matrix< M, R, T >& u_ );
		
        static void svd_mode1( const t3_type& data_, u1_type& u1_ );
        static void svd_mode2( const t3_type& data_, u2_type& u2_ );
        static void svd_mode3( const t3_type& data_, u3_type& u3_ );
 		
		//hosvd on eigenvalue decomposition = hoeigs
		template< size_t N, size_t R  >
        static void get_eigs_u_red( const matrix< N, N, T >& data_, matrix< N, R, T >& u_ );
		
		static void eigs_mode1( const t3_type& data_, u1_type& u1_ );
        static void eigs_mode2( const t3_type& data_, u2_type& u2_ );
        static void eigs_mode3( const t3_type& data_, u3_type& u3_ );
		
	}; //end hosvd class
	
#define VMML_HOSVD_TEMPLATE_STRING        template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_HOSVD_TEMPLATE_CLASSNAME     t3_hosvd< R1, R2, R3, I1, I2, I3, T >


	
VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::hosvd( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
{
	svd_mode1( data_, u1_ );
	svd_mode2( data_, u2_ );
	svd_mode3( data_, u3_ );
}
	
VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::hoeigs( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
{
	eigs_mode1( data_, u1_ );
	eigs_mode2( data_, u2_ );
	eigs_mode3( data_, u3_ );
}
	
VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::apply_all( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
{
	apply_mode1( data_, u1_ );
	apply_mode2( data_, u2_ );
	apply_mode3( data_, u3_ );
}	
	
VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::apply_mode1( const t3_type& data_, u1_type& u1_ )
{
#ifdef HOEIGS
	eigs_mode1( data_, u1_ );
#else
	svd_mode1( data_, u1_ );
#endif
}	


VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::apply_mode2( const t3_type& data_, u2_type& u2_ )
{
#ifdef HOEIGS
	eigs_mode2( data_, u2_ );
#else
	svd_mode2( data_, u2_ );
#endif
}



VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::apply_mode3( const t3_type& data_, u3_type& u3_ )
{
#ifdef HOEIGS
	eigs_mode3( data_, u3_ );
#else
	svd_mode3( data_, u3_ );
#endif
}
	
/* SVD along mode 1, 2, and 3*/ 
	
VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::svd_mode1( const t3_type& data_, u1_type& u1_ )
{
	u1_unfolded_type* u = new u1_unfolded_type; // -> u1
	data_.lateral_unfolding_bwd( *u );
	
	get_svd_u_red( *u, u1_ );
	
	delete u;
}	

VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::svd_mode2( const t3_type& data_, u2_type& u2_ )
{
	u2_unfolded_type* u = new u2_unfolded_type; // -> u1
	data_.frontal_unfolding_bwd( *u );
	
	get_svd_u_red( *u, u2_ );
	
	delete u;
}

VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::svd_mode3( const t3_type& data_, u3_type& u3_ )
{
	u3_unfolded_type* u = new u3_unfolded_type; // -> u1
	data_.horizontal_unfolding_bwd( *u );
	
	get_svd_u_red( *u, u3_ );
	
	delete u;
}
	
/* EIGS for mode 1, 2 and 3*/ 	
	
VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::eigs_mode1( const t3_type& data_, u1_type& u1_ )
{
	//unfolding / matricization
	u1_unfolded_type* m_lateral = new u1_unfolded_type; // -> u1
	data_.lateral_unfolding_bwd( *m_lateral );
	
	//covariance matrix of unfolded data
	u1_cov_type* cov  = new u1_cov_type;
	m_lateral->symmetric_covariance( *cov );
	delete m_lateral;
	
	//compute x largest magnitude eigenvalues; x = R
	get_eigs_u_red( *cov, u1_ );
	
	delete cov;
}	

VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::eigs_mode2( const t3_type& data_, u2_type& u2_ )
{
	//unfolding / matricization
	u2_unfolded_type* m_frontal = new u2_unfolded_type; // -> u2
	data_.frontal_unfolding_bwd( *m_frontal );
	
	//covariance matrix of unfolded data
	u2_cov_type* cov  = new u2_cov_type;
	m_frontal->symmetric_covariance( *cov );
	delete m_frontal;
	
	//compute x largest magnitude eigenvalues; x = R
	get_eigs_u_red( *cov, u2_ );
	
	delete cov;
}	

VMML_HOSVD_TEMPLATE_STRING
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::eigs_mode3( const t3_type& data_, u3_type& u3_)
{
	//unfolding / matricization
	u3_unfolded_type* m_horizontal = new u3_unfolded_type; // -> u3
	data_.horizontal_unfolding_bwd( *m_horizontal );
	
	//covariance matrix of unfolded data
	u3_cov_type* cov  = new u3_cov_type;
#if 1
	m_horizontal->symmetric_covariance( *cov );
#else	
	blas_dgemm< I3, I1*I2, I3, T>* blas_cov = new blas_dgemm< I3, I1*I2, I3, T>;
	blas_cov->compute( *m_horizontal, *cov );
#endif
	
	delete m_horizontal;
	
	//compute x largest magnitude eigenvalues; x = R
	get_eigs_u_red( *cov, u3_ );
	
	delete cov;
}	
	
	
	
/* helper methods for SVD and EIGS*/	
	
VMML_HOSVD_TEMPLATE_STRING
template< size_t N, size_t R >
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::get_eigs_u_red( const matrix< N, N, T >& data_, matrix< N, R, T >& u_ )
{
	typedef matrix< N, N, T_svd > cov_matrix_type;
	typedef vector< R, T_svd > eigval_type;
	typedef	matrix< N, R, T_svd > eigvec_type;
	//typedef	matrix< N, R, T_coeff > coeff_type;
	
	//compute x largest magnitude eigenvalues; x = R
	eigval_type* eigxvalues =  new eigval_type;
	eigvec_type* eigxvectors = new eigvec_type; 
	
	lapack_sym_eigs< N, T_svd >  eigs;
	cov_matrix_type* data = new cov_matrix_type;
	data->cast_from( data_ );
	if( eigs.compute_x( *data, *eigxvectors, *eigxvalues) ) {
		
		/*if( _is_quantify_coeff ){
			coeff_type* evec_quant = new coeff_type; 
			T min_value = 0; T max_value = 0;
			u_.cast_from( *eigxvectors );
			u_.quantize( *evec_quant, min_value, max_value );
			evec_quant->dequantize( u_, min_value, max_value );
			delete evec_quant;
		} else */ if ( sizeof( T ) != 4 ){
			u_.cast_from( *eigxvectors );
		} else {
			u_ = *eigxvectors;
		}
		
	} else {
		u_.zero();
	}
	
	delete eigxvalues;
	delete eigxvectors;
	delete data;
	
}

VMML_HOSVD_TEMPLATE_STRING
template< size_t M, size_t N, size_t R >
void 
VMML_HOSVD_TEMPLATE_CLASSNAME::get_svd_u_red( const matrix< M, N, T >& data_, matrix< M, R, T >& u_ )
{
	typedef	matrix< M, N, T_svd > svd_m_type;
	//FIXME: typedef	matrix< M, N, T_coeff > coeff_type;
	typedef	matrix< M, N, T > m_type;
	typedef vector< N, T_svd > lambdas_type;
	
	svd_m_type* u_double = new svd_m_type; 
	u_double->cast_from( data_ );
	
	//FIXME:: coeff_type* u_quant = new coeff_type; 
	m_type* u_out = new m_type; 
	
	lambdas_type* lambdas  = new lambdas_type;
	lapack_svd< M, N, T_svd >* svd = new lapack_svd<  M, N, T_svd >();
	if( svd->compute_and_overwrite_input( *u_double, *lambdas )) {
		
		/*		if( _is_quantify_coeff ){
		 T min_value = 0; T max_value = 0;
		 u_comp->cast_from( *u_double );
		 //FIXME: u_comp->quantize( *u_quant, min_value, max_value );
		 //FIXME: u_quant->dequantize( *u_internal, min_value, max_value );
		 } else */ if ( sizeof( T ) != 4 ){
			 u_out->cast_from( *u_double );
		 } else {
			 *u_out = *u_double;
		 }
		
		u_out->get_sub_matrix( u_ );
		
	} else {
		u_.zero();
	}
	
	delete lambdas;
	delete svd;
	delete u_double;
	//delete u_quant;
	delete u_out;
}	

#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME

}//end vmml namespace

#endif