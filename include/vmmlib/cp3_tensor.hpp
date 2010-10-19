/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 * @author Jonas Boesch
 *
 * The cp3 tensor class is consists of three basis matrices u1-u3 and R lambda values for a given rank-R approximation
 * CP stands for Candecomp/Parafac (1970)
 * - Carroll, J. and Chang, Jih-Jie: Analysis of Individual Differences in Multidimensional Scaling via an N-way generalization of ``Eckart--Young'' decompositions,
 * Psychometrika 35: 283–319, 1970
 * - R. A. Harshman: Foundations of the PARAFAC procedure: Models and conditions for an 'explanatory' multi-modal factor analysis,
 * UCLA Working Papers in Phonetics, Vol. 16, No. 1., 1970
 * - De Lathauwer L., De Moor B., Vandewalle J., ``On the Best rank-1 and Rank-$(R_1, R_2, ..., R_N)$ Approximation and Applications of Higher-Order Tensors'', 
 * SIAM J. Matrix Anal. Appl., vol. 21, no. 4, Apr. 2000, pp. 1324-1342.
 * - T. G. Kolda and B. W. Bader. Tensor Decompositions and Applications. 
 * SIAM Review, Volume 51, Number 3, Pages 455-500, September 2009.
 * 
 */

#ifndef __VMML__CP3_TENSOR__HPP__
#define __VMML__CP3_TENSOR__HPP__

#include <vmmlib/tensor3.hpp>
#include <vmmlib/lapack_svd.hpp>
#include <vmmlib/matrix_pseudoinverse.hpp>

namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, size_t R, typename T = float >
	class cp3_tensor
	{
	public:    
		cp3_tensor( matrix< I1, R, T >& U1, matrix< I2, R, T >& U2, matrix< I3, R, T >& U3, vector< R, T >& lambdas_ );
		
		void set_lambdas( const vector< R, T >& lambdas_ )  { _lambdas = lambdas_; } ;
		void set_u1( const matrix< I1, R, T >& U1 ) { _u1 = U1; } ;
		void set_u2( const matrix< I2, R, T >& U2 ) { _u2 = U2; } ;
		void set_u3( const matrix< I3, R, T >& U3 ) { _u3 = U3; } ;
		
		vector< R, T > get_lambdas() const { return _lambdas; } ;
		matrix< I1, R, T > get_u1() const { return _u1; } ;
		matrix< I2, R, T > get_u2() const { return _u2; } ;
		matrix< I3, R, T > get_u3() const { return _u3; } ;
		
		void decomposition( const tensor3< I1, I2, I3, T >& data_ ); 
		void reconstruction( tensor3< I1, I2, I3, T >& data_ ) const;
		
		void cp_als( const tensor3< I1, I2, I3, T >& data_ );
		void hoii( const tensor3< I1, I2, I3, T >& data_ );
		
		void hosvd_mode1( const tensor3< I1, I2, I3, T >& data_, matrix< I1, R, T >& U1_ ) const;
		void hosvd_mode2( const tensor3< I1, I2, I3, T >& data_, matrix< I2, R, T >& U2_ ) const;
		void hosvd_mode3( const tensor3< I1, I2, I3, T >& data_, matrix< I3, R, T >& U3_ ) const;
		
	private:
		vector< R, T > _lambdas ;
		matrix< I1, R, T > _u1 ;
		matrix< I2, R, T > _u2 ;
		matrix< I3, R, T > _u3 ;
		
	}; // class cp3_tensor
	
	
#define VMML_TEMPLATE_STRING    	template< size_t I1, size_t I2, size_t I3, size_t R, typename T >
#define VMML_TEMPLATE_CLASSNAME     cp3_tensor< I1, I2, I3, R, T >
	
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::cp3_tensor( matrix< I1, R, T >& U1, matrix< I2, R, T >& U2, matrix< I3, R, T >& U3, vector< R, T >& lambdas_ )
: _lambdas(lambdas_), _u1(U1), _u2(U2), _u3(U3)
{
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruction( tensor3< I1, I2, I3, T >& data_ ) const
{
	//data_.full_tensor3_matrix_multiplication( _core, _u1, _u2, _u3 );
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::decomposition( const tensor3< I1, I2, I3, T >& data_ )
{
	cp_als( data_ );
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::cp_als( const tensor3< I1, I2, I3, T >& data_ )
{
	hoii( data_ );
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hoii( const tensor3< I1, I2, I3, T >& data_ )
{
	
	//compute best rank-(R) approximation (Lathauwer et al., 2000b)
	tensor3< I1, I2, I3, T > approximated_data;
	reconstruction( approximated_data );
	double f_norm = approximated_data.compute_frobenius_norm();
	double max_f_norm = data_.compute_frobenius_norm();
	//std::cout << "frobenius norm original: " << max_f_norm << std::endl;
	
	double last_f_norm = f_norm;
	double improvement = max_f_norm - f_norm;
	double min_improvement = 0.1;
	size_t i = 0;
	size_t max_iterations = 3;
	
		
	//optimize for mode 1
	//intialize u1
	hosvd_mode1( data_, _u1 );
	std::cout << "initial u1: " << std::endl << _u1 << std::endl;
	while( improvement > min_improvement && i < max_iterations )
	{
		
		
		set_u1( _u1 );
		set_u2( _u2 );
		set_u3( _u3 );
		set_lambdas( _lambdas );
		
		reconstruction( approximated_data );
		f_norm = approximated_data.compute_frobenius_norm();
		improvement = f_norm - last_f_norm;
		last_f_norm = f_norm;
		
		//std::cout << "iteration '" << i << "': frobenius norm: " << std::setprecision(8) << f_norm << ", improvement: " << improvement << std::endl;
		
		++i;
	}
	//std::cout << "number of iterations: " << i << std::endl;
	
		/*tensor3< I1, J2, J3, T > projection1; 
		optimize_mode1( data_, projection1, _u2, _u3);
		hosvd_mode1( projection1, _u1 );
		
		//optimize for mode 2
		tensor3< J1, I2, J3, T > projection2; 
		optimize_mode2( data_, projection2, _u1, _u3);
		hosvd_mode2( projection2, _u2 );
		
		//optimize for mode 3
		tensor3< J1, J2, I3, T > projection3; 
		optimize_mode3( data_, projection3, _u1, _u2);
		hosvd_mode3( projection3, _u3);*/
		


}

	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode1( const tensor3< I1, I2, I3, T >& data_, matrix< I1, R, T >& U1_ ) const
{
	matrix< I1, I2*I3, T> m_lateral; // -> u1
	data_.lateral_matricization( m_lateral);
	
	//std::cout << "hosvd mode1, m_lateral: " << std::endl << m_lateral << std::endl;
	
	vector< I2*I3, double > lambdas_u1;
	
	lapack_svd< I1, I2*I3, double > svd1;
	if( svd1.compute_and_overwrite_input( m_lateral, lambdas_u1 ))
		m_lateral.get_sub_matrix( U1_ );
	else 
		U1_.zero();
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode2( const tensor3< I1, I2, I3, T >& data_, matrix< I2, R, T >& U2_ ) const
{
	matrix< I2, I1*I3, T> m_frontal; // -> u2
	data_.frontal_matricization( m_frontal);
	//std::cout << "hosvd mode2, m_frontal: " << std::endl << m_frontal << std::endl;
	
	vector< I1*I3, double > lambdas_u2;
	
	lapack_svd< I2, I1*I3, double > svd2;
	if( svd2.compute_and_overwrite_input( m_frontal, lambdas_u2 ))
		m_frontal.get_sub_matrix( U2_ );
	else 
		U2_.zero();
}



VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode3( const tensor3< I1, I2, I3, T >& data_, matrix< I3, R, T >& U3_ ) const
{
	matrix< I3, I1*I2, T> m_horizontal; //-> u3
	data_.horizontal_matricization( m_horizontal);
	//std::cout << "hosvd mode3, m_horizontal: " << std::endl << m_horizontal << std::endl;
	
	vector< I1*I2, double > lambdas_u3;
	lapack_svd< I3, I1*I2, double > svd3;
	if( svd3.compute_and_overwrite_input( m_horizontal, lambdas_u3 ))
		m_horizontal.get_sub_matrix( U3_ );
	else 
		U3_.zero();
}
	
		
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
		
	} // namespace vmml
	
#endif
		
