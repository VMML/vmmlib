/* 
 * VMMLib - CUDA Tensor Classes
 *  
 * @author Susanne Suter
 *
 * HOSVD based on covariance matrix (using CUBLAS DGEMM) of unfolded input and successive EIGS on covariance matrix
 *
 * The Tucker3 tensor class is consists of the same components (core tensor, basis matrices u1-u3) as the tucker3 model described in:
 * - Tucker, 1966: Some mathematical notes on three-mode factor analysis, Psychometrika.
 * - Kroonenberg & De Leeuw, 1980: Principal component analysis of three-mode data by means of alternating least squares algorithms. Psychometrika. (TUCKALS)
 * - De Lathauwer, De Moor, Vandewalle, 2000a: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * - Bader & Kolda, 2006: Algorithm 862: Matlab tensor classes for fast algorithm prototyping. ACM Transactions on Mathematical Software.
 * 
 * TODO:: use CULA Tools for EIGS
 * 
 */

#ifndef __VMML__T3_CUBLAS_HOSVD__HPP__
#define __VMML__T3_CUBLAS_HOSVD__HPP__

#include <vmmlib/tensor3.hpp>
#include <vmmlib/lapack_svd.hpp>
#include <vmmlib/lapack_sym_eigs.hpp>
#include <vmmlib/cublas_dgemm.cu>


namespace vmml
{
	
	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3 >
	class t3_cublas_hosvd
	{
	public:    
		
		typedef double T_svd;
		
		typedef tensor3< I1, I2, I3, float > t3_type;
		
		typedef matrix< I1, R1, float > u1_type;
		typedef matrix< I2, R2, float > u2_type;
		typedef matrix< I3, R3, float > u3_type;
		
		typedef matrix< I1, I2*I3, float > u1_unfolded_type;
		typedef matrix< I2, I1*I3, float > u2_unfolded_type;
		typedef matrix< I3, I1*I2, float > u3_unfolded_type;
		
		typedef matrix< I1, I1, float > u1_cov_type;
		typedef matrix< I2, I2, float > u2_cov_type;
		typedef matrix< I3, I3, float > u3_cov_type;
		
		/*	higher-order singular value decomposition (HOSVD) with full rank decomposition (also known as Tucker decomposition). 
		 see: De Lathauer et al, 2000a: A multilinear singular value decomposition. 
		 the hosvd can be computed (a) with n-mode PCA, i.e., an eigenvalue decomposition on the covariance matrix of every mode's matricization, and 
		 other known names for HOSVD: n-mode SVD, 3-mode factor analysis (3MFA, tucker3), 3M-PCA, n-mode PCA, higher-order SVD
		 */

		static void apply_mode1( const t3_type& data_, u1_type& u1_ );
		static void apply_mode2( const t3_type& data_, u2_type& u2_ );
		static void apply_mode3( const t3_type& data_, u3_type& u3_ );
		
		static void apply_all( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ );
		static void hoeigs( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ );
		
	protected:
		
		//hosvd on eigenvalue decomposition = hoeigs
		template< size_t N, size_t R  >
        static void get_eigs_u_red( const matrix< N, N, float >& data_, matrix< N, R, float >& u_ );
		
		static void eigs_mode1( const t3_type& data_, u1_type& u1_ );
        static void eigs_mode2( const t3_type& data_, u2_type& u2_ );
        static void eigs_mode3( const t3_type& data_, u3_type& u3_ );
		
	}; //end hosvd class
	
#define VMML_TEMPLATE_STRING        template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3 >
#define VMML_TEMPLATE_CLASSNAME     t3_cublas_hosvd< R1, R2, R3, I1, I2, I3 >


	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hoeigs( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
{
	eigs_mode1( data_, u1_ );
	eigs_mode2( data_, u2_ );
	eigs_mode3( data_, u3_ );
}
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::apply_all( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
{
	apply_mode1( data_, u1_ );
	apply_mode2( data_, u2_ );
	apply_mode3( data_, u3_ );
}	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::apply_mode1( const t3_type& data_, u1_type& u1_ )
{
	eigs_mode1( data_, u1_ );
}	


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::apply_mode2( const t3_type& data_, u2_type& u2_ )
{
	eigs_mode2( data_, u2_ );
}



VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::apply_mode3( const t3_type& data_, u3_type& u3_  )
{
	eigs_mode3( data_, u3_ );
}
	
	
/* EIGS for mode 1, 2 and 3*/ 	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::eigs_mode1( const t3_type& data_, u1_type& u1_ )
{
	//unfolding / matricization
	u1_unfolded_type* m_lateral = new u1_unfolded_type; // -> u1
	data_.lateral_unfolding_bwd( *m_lateral );
	
	//covariance matrix of unfolded data
	u1_cov_type* cov  = new u1_cov_type;
	cublas_dgemm< I1, I2*I3, I1, float >* blas_cov = new cublas_dgemm< I1, I2*I3, I1, float >;
	blas_cov->compute( *m_lateral, *cov );
	delete blas_cov;
	delete m_lateral;

	//compute x largest magnitude eigenvalues; x = R
	get_eigs_u_red( *cov, u1_ );
	
	delete cov;
}	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::eigs_mode2( const t3_type& data_, u2_type& u2_ )
{
	//unfolding / matricization
	u2_unfolded_type* m_frontal = new u2_unfolded_type; // -> u2
	data_.frontal_unfolding_bwd( *m_frontal );
	
	//covariance matrix of unfolded data
	u2_cov_type* cov  = new u2_cov_type;
	cublas_dgemm< I2, I1*I3, I2, float >* blas_cov = new cublas_dgemm< I2, I1*I3, I2, float >;
	blas_cov->compute( *m_frontal, *cov );
	delete blas_cov;
	delete m_frontal;
	
	//compute x largest magnitude eigenvalues; x = R
	get_eigs_u_red( *cov, u2_ );
	
	delete cov;
}	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::eigs_mode3( const t3_type& data_, u3_type& u3_)
{
	//unfolding / matricization
	u3_unfolded_type* m_horizontal = new u3_unfolded_type; // -> u3
	data_.horizontal_unfolding_bwd( *m_horizontal );
	
	//covariance matrix of unfolded data
	u3_cov_type* cov  = new u3_cov_type;
	cublas_dgemm< I3, I1*I2, I3, float >* blas_cov = new cublas_dgemm< I3, I1*I2, I3, float >;
	blas_cov->compute( *m_horizontal, *cov );
	delete blas_cov;
	delete m_horizontal;
	
	//compute x largest magnitude eigenvalues; x = R
	get_eigs_u_red( *cov, u3_ );
	
	delete cov;
}	
	
	
	
/* helper methods for EIGS*/	
	
VMML_TEMPLATE_STRING
template< size_t N, size_t R >
void 
VMML_TEMPLATE_CLASSNAME::get_eigs_u_red( const matrix< N, N, float >& data_, matrix< N, R, float >& u_ )
{
	typedef matrix< N, N, T_svd > cov_matrix_type;
	typedef vector< R, T_svd > eigval_type;
	typedef	matrix< N, R, T_svd > eigvec_type;
	//typedef	matrix< N, R, T_coeff > coeff_type;
	
	//compute x largest magnitude eigenvalues; x = R
	eigval_type* eigxvalues =  new eigval_type;
	eigvec_type* eigxvectors = new eigvec_type; 
	
	
	//TODO: could try to use CULA TOOLS
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
		} else */ 
		u_ = *eigxvectors;
		
	} else {
		u_.zero();
	}
	
	delete eigxvalues;
	delete eigxvectors;
	delete data;
	
}

#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME

}//end vmml namespace

#endif

