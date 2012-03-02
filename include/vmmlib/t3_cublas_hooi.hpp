/* 
 * VMMLib - CUDA Tensor Classes
 *  
 * @author Susanne Suter
 *
 * The higher-order orthogonal iteration (HOOI) is also known as Tucker-ALS (Tuck-ALS)
 * The t3_cublas_hooi implements a HOOI for a third-order tensor using CUBLAS
 * references:
 * - Tucker, 1966: Some mathematical notes on three-mode factor analysis, Psychometrika.
 * - De Lathauwer, De Moor, Vandewalle, 2000a: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - De Lathauwer, De Moor, Vandewalle, 2000b: On the Best rank-1 and Rank-(R_1, R_2, ..., R_N) Approximation and Applications of Higher-Order Tensors, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * - Bader & Kolda, 2006: Algorithm 862: Matlab tensor classes for fast algorithm prototyping. ACM Transactions on Mathematical Software.
 *
 */

#ifndef __VMML__T3_CUBLAS_HOOI__HPP__
#define __VMML__T3_CUBLAS_HOOI__HPP__


#include <vmmlib/t3_cublas_hosvd.hpp>
#include <vmmlib/t3_cublas_ttm.hpp>
#include <vmmlib/matrix_pseudoinverse.hpp>

namespace vmml
{
	
	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3 >
	class t3_cublas_hooi
	{
	public:    
		
		typedef tensor3< I1, I2, I3, float > t3_type;
		
		typedef tensor3< R1, R2, R3, float > t3_core_type;
		
		typedef matrix< I1, R1, float > u1_type;
		typedef matrix< I2, R2, float > u2_type;
		typedef matrix< I3, R3, float > u3_type;
		
		typedef matrix< R1, I1, float > u1_t_type;
		typedef matrix< R2, I2, float > u2_t_type;
		typedef matrix< R3, I3, float > u3_t_type;

		/*	higher-order orthogonal iteration (HOOI) is a truncated HOSVD decompositions, i.e., the HOSVD components are of lower-ranks. An optimal rank-reduction is 
		 performed with an alternating least-squares (ALS) algorithm, which minimizes the error between the approximated and orignal tensor based on the Frobenius norm
		 see: De Lathauwer et al, 2000b; On the best rank-1 and rank-(RRR) approximation of higher-order tensors.
		 the HOOI can be computed based on (a) n-mode PCA, i.e., an eigenvalue decomposition on the covariance matrix of every mode's matriciziation, and 
		 (b) by performing a 2D SVD on the matricization of every mode. Matrix matricization means that a tensor I1xI2xI3 is unfolded/sliced into one matrix
		 with the dimensions I1xI2I3, which corresponds to a matrizitation alonge mode I1.
		 */
		template< typename T_init>
		static void als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, t3_core_type& core_, T_init init );

		//core not needed
		template< typename T_init>
		static void als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, T_init init );
		

		/* derive core
		 implemented accodring to core = data x_1 U1_pinv x_2 U2_pinv x_3 U3_pinv, 
		 where x_1 ... x_3 are n-mode products and U1_pinv ... U3_pinv are inverted basis matrices
		 the inversion is done with a matrix pseudoinverse computation
		 */
		static void derive_core(  const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_ );
		//faster: but only if basis matrices are orthogonal
		static void derive_core_orthogonal_bases(  const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_  );
		
		// init functors 
		struct init_hosvd
		{
			inline void operator()( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
			{
				t3_cublas_hosvd< R1, R2, R3, I1, I2, I3 >::apply_mode1( data_, u1_ );
				t3_cublas_hosvd< R1, R2, R3, I1, I2, I3 >::apply_mode2( data_, u2_ );
				t3_cublas_hosvd< R1, R2, R3, I1, I2, I3 >::apply_mode3( data_, u3_ );
			}
		};
		
		struct init_random
		{
			inline void operator()( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_ )
			{
				srand( time(NULL) );
				u1_.set_random();
				u2_.set_random();
				u3_.set_random();
				
				u1_ /= u1_.frobenius_norm();
				u2_ /= u2_.frobenius_norm();
				u3_ /= u3_.frobenius_norm();
			}
		};
		
		
		
	protected:
	
			
        static void optimize_mode1( const t3_type& data_, const u2_type& u2_, const u3_type& u3_, tensor3< I1, R2, R3, float >& projection_ );
        static void optimize_mode2( const t3_type& data_, const u1_type& u1_, const u3_type& u3_, tensor3< R1, I2, R3, float >& projection_ );		
        static void optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, tensor3< R1, R2, I3, float >& projection_ );
		
	};//end class t3_cublas_hooi
	
#define VMML_TEMPLATE_STRING        template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3 >
#define VMML_TEMPLATE_CLASSNAME     t3_cublas_hooi< R1, R2, R3, I1, I2, I3 >

	
VMML_TEMPLATE_STRING
template< typename T_init>
void 
VMML_TEMPLATE_CLASSNAME::als( const t3_type& data_, 
							 u1_type& u1_, u2_type& u2_, u3_type& u3_, 
							 T_init init )
{
	t3_core_type core;
	core.zero();
	als( data_, u1_, u2_, u3_, core, init );
}	
	
	
	
VMML_TEMPLATE_STRING
template< typename T_init>
void 
VMML_TEMPLATE_CLASSNAME::als( const t3_type& data_, 
								  u1_type& u1_, u2_type& u2_, u3_type& u3_, 
								  t3_core_type& core_,
								  T_init init )
{
	//intialize basis matrices
	init( data_, u1_, u2_, u3_ );
	
	//derve core from initialized matrices
	derive_core_orthogonal_bases( data_, u1_, u2_, u3_, core_ );
	
	//removed initial fit to save computing time
	//compute best rank-(R1, R2, R3) approximation (Lathauwer et al., 2000b)
	//t3_type approximated_data;
	//t3_cublas_ttm::full_tensor3_matrix_multiplication( core_, u1_, u2_, u3_, approximated_data );
	
	double f_norm = 0;//approximated_data.frobenius_norm();
	double max_f_norm = data_.frobenius_norm();
	double normresidual  =  0; //sqrt( (max_f_norm * max_f_norm) - (f_norm * f_norm));
	double fit = 0;
	/*if ( (max_f_norm != 0) && (max_f_norm > f_norm) ) 
	{
		fit = 1 - (normresidual / max_f_norm);
	} else { 
		fit = 1;
	}*/
	
	double fitchange = 1;
	double fitold = fit;
	double fitchange_tolerance = 1.0e-4;
	
	tensor3< I1, R2, R3, float > projection1; 
	tensor3< R1, I2, R3, float > projection2; 
	tensor3< R1, R2, I3, float > projection3; 
	
#if TUCKER_LOG
	std::cout << "HOOI ALS (for tensor3) " << std::endl 
	<< "initial fit: " << fit  << ", "
	<< "frobenius norm original: " << max_f_norm << std::endl;
#endif	
	size_t i = 0;
	size_t max_iterations = 10;
	while( (fitchange >= fitchange_tolerance) && (i < max_iterations) )
	{
		fitold = fit;
		
		//optimize modes
		optimize_mode1( data_, u2_, u3_, projection1 );
		t3_cublas_hosvd< R1, R2, R3, I1, R2, R3 >::apply_mode1( projection1, u1_ );
		
		optimize_mode2( data_, u1_, u3_, projection2 );
		t3_cublas_hosvd< R1, R2, R3, R1, I2, R3 >::apply_mode2( projection2, u2_ );
		
		optimize_mode3( data_, u1_, u2_, projection3 );
		t3_cublas_hosvd< R1, R2, R3, R1, R2, I3 >::apply_mode3( projection3, u3_ );
		
		t3_cublas_ttm::multiply_horizontal_bwd( projection3, transpose( u3_ ), core_ );
		f_norm = core_.frobenius_norm();
		normresidual  = sqrt( max_f_norm * max_f_norm - f_norm * f_norm);
		fit = 1 - (normresidual / max_f_norm);
		fitchange = fabs(fitold - fit);
		
#if TUCKER_LOG
		std::cout << "iteration '" << i << "', fit: " << fit 
		<< ", fitdelta: " << fitchange 
		<< ", frobenius norm: " << f_norm << std::endl;		
#endif
		++i;
	}
}	



VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode1( const t3_type& data_, const u2_type& u2_, const u3_type& u3_, tensor3< I1, R2, R3, float >& projection_ )
{
	u2_t_type* u2_inv = new u2_t_type;
	u3_t_type* u3_inv = new u3_t_type;
	u2_.transpose_to( *u2_inv );
	u3_.transpose_to( *u3_inv );
	
	//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
	tensor3< I1, R2, I3, float >* tmp  = new tensor3< I1, R2, I3, float >;
	t3_cublas_ttm::multiply_frontal_bwd( data_, *u2_inv, *tmp );
	t3_cublas_ttm::multiply_horizontal_bwd( *tmp, *u3_inv, projection_ );
	
	delete u2_inv;
	delete u3_inv;
	delete tmp;
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode2( const t3_type& data_, const u1_type& u1_, const u3_type& u3_, tensor3< R1, I2, R3, float >& projection_ )
{
	u1_t_type* u1_inv = new u1_t_type();
	u3_t_type* u3_inv = new u3_t_type();
	u1_.transpose_to( *u1_inv );
	u3_.transpose_to( *u3_inv );
	
	
	//backward cyclic matricization (after Lathauwer et al., 2000a)
	tensor3< R1, I2, I3, float >* tmp = new tensor3< R1, I2, I3, float >();
	t3_cublas_ttm::multiply_lateral_bwd( data_, *u1_inv, *tmp );
	t3_cublas_ttm::multiply_horizontal_bwd( *tmp, *u3_inv, projection_ );
	
	delete u1_inv;
	delete u3_inv;
	delete tmp;
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode3( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, tensor3< R1, R2, I3, float >& projection_ )
{
	u1_t_type* u1_inv = new u1_t_type;
	u2_t_type* u2_inv = new u2_t_type;
	u1_.transpose_to( *u1_inv );
	u2_.transpose_to( *u2_inv );
	
	//backward cyclic matricization (after Lathauwer et al., 2000a)
	tensor3< R1, I2, I3, float >* tmp = new tensor3< R1, I2, I3, float >();
	t3_cublas_ttm::multiply_lateral_bwd( data_, *u1_inv, *tmp );
	t3_cublas_ttm::multiply_frontal_bwd( *tmp, *u2_inv, projection_ );
	
	delete u1_inv;
	delete u2_inv;
	delete tmp;
}
	
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::derive_core_orthogonal_bases( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_ )
{
	u1_t_type* u1_inv = new u1_t_type;
	u2_t_type* u2_inv = new u2_t_type;
	u3_t_type* u3_inv = new u3_t_type;
	
	u1_.transpose_to( *u1_inv );
	u2_.transpose_to( *u2_inv );
	u3_.transpose_to( *u3_inv );
	
	t3_cublas_ttm::full_tensor3_matrix_multiplication( data_, *u1_inv, *u2_inv, *u3_inv, core_ );
	
	delete u1_inv;
	delete u2_inv;
	delete u3_inv;
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::derive_core( const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_ )
{
	//compute pseudo inverse for matrices u1-u3
	u1_type* u1_pinv_t = new u1_type;
	u2_type* u2_pinv_t = new u2_type;
	u3_type* u3_pinv_t = new u3_type;
	
	compute_pseudoinverse<  u1_type > compute_pinv_u1;
	compute_pinv_u1( u1_, *u1_pinv_t );
	compute_pseudoinverse<  u2_type > compute_pinv_u2;
	compute_pinv_u2( u2_, *u2_pinv_t );	
	compute_pseudoinverse<  u3_type > compute_pinv_u3;
	compute_pinv_u3( u3_, *u3_pinv_t );
	
	u1_t_type* u1_pinv = new u1_t_type;
	u2_t_type* u2_pinv = new u2_t_type;
	u3_t_type* u3_pinv = new u3_t_type;
	
	u1_pinv_t->transpose_to( *u1_pinv );
	u2_pinv_t->transpose_to( *u2_pinv );
	u3_pinv_t->transpose_to( *u3_pinv );
		
	t3_cublas_ttm::full_tensor3_matrix_multiplication( data_, *u1_pinv, *u2_pinv, *u3_pinv, core_ );
	
	delete u1_pinv;
	delete u2_pinv;
	delete u3_pinv;
	delete u1_pinv_t;
	delete u2_pinv_t;
	delete u3_pinv_t;
}

	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
	
}//end vmml namespace

#endif

