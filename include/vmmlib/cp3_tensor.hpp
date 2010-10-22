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
#include <vmmlib/tensor3_iterator.hpp>
#include <vmmlib/lapack_svd.hpp>
#include <vmmlib/matrix_pseudoinverse.hpp>

namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, size_t R, typename T = float >
	class cp3_tensor
	{
	public:    
		cp3_tensor( matrix< I1, R, T >& U1, matrix< I2, R, T >& U2, matrix< I3, R, T >& U3, vector< R, T >& lambdas_ );
		
		typedef matrix< I1, R, T > u1_type;
		typedef typename u1_type::iterator u1_iterator;
		typedef typename u1_type::const_iterator u1_const_iterator;
		
		typedef matrix< I2, R, T > u2_type;
		typedef typename u2_type::iterator u2_iterator;
		typedef typename u2_type::const_iterator u2_const_iterator;
		
		typedef matrix< I3, R, T > u3_type;
		typedef typename u3_type::iterator u3_iterator;
		typedef typename u3_type::const_iterator u3_const_iterator;
		
		//TODO typedef for m_lateral_type, m_frontal_type, and m_horizontal_type;
		

		void set_lambdas( const vector< R, T >& lambdas_ )  { _lambdas = lambdas_; } ;
		void set_u1( const u1_type& U1 ) { _u1 = U1; } ;
		void set_u2( const u2_type& U2 ) { _u2 = U2; } ;
		void set_u3( const u3_type& U3 ) { _u3 = U3; } ;
		
		void get_lambdas( vector< R, T >& data_ ) const { data_  = _lambdas; } ;
		void get_u1( u1_type& U1 ) const { U1 = _u1; } ;
		void get_u2( u2_type& U2 ) const { U2 = _u2; } ;
		void get_u3( u3_type& U3 ) const { U3 = _u3; } ;
		
		void export_to( std::vector< T >& data_ ) const;
		void import_from( std::vector< T >& data_ );	
		
		void decomposition( const tensor3< I1, I2, I3, T >& data_ ); 
		void reconstruction( tensor3< I1, I2, I3, T >& data_ ) const;
		
		void cp_als( const tensor3< I1, I2, I3, T >& data_ );
		
		//higher-order power method (lathauwer et al., 2000b)
		void hopm( const tensor3< I1, I2, I3, T >& data_ );
		
		void hosvd_mode1( const tensor3< I1, I2, I3, T >& data_, u1_type& U1_ ) const;
		void hosvd_mode2( const tensor3< I1, I2, I3, T >& data_, u2_type& U2_ ) const;
		void hosvd_mode3( const tensor3< I1, I2, I3, T >& data_, u3_type& U3_ ) const;
		
		void optimize_mode1( const tensor3< I1, I2, I3, T >& data_, u1_type& U1_optimized_, const u2_type& U2_, const u3_type& U3_ ) const;
		void optimize_mode2( const tensor3< I1, I2, I3, T >& data_, const u1_type& U1_, u2_type& U2_optimized_, const u3_type& U3_ ) const;		
		double optimize_mode3( const tensor3< I1, I2, I3, T >& data_, const u1_type& U1_, const u2_type& U2_, u3_type& U3_optimized_ ) const;
		
	private:
		vector< R, T > _lambdas ;
		u1_type _u1 ;
		u2_type _u2 ;
		u3_type _u3 ;
		
	}; // class cp3_tensor
	
	
#define VMML_TEMPLATE_STRING    	template< size_t I1, size_t I2, size_t I3, size_t R, typename T >
#define VMML_TEMPLATE_CLASSNAME     cp3_tensor< I1, I2, I3, R, T >
	
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::cp3_tensor( u1_type& U1, u2_type& U2, u3_type& U3, vector< R, T >& lambdas_ )
: _lambdas(lambdas_), _u1(U1), _u2(U2), _u3(U3)
{
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruction( tensor3< I1, I2, I3, T >& data_ ) const
{
	tensor3< R, R, R, T > core_diag;
	core_diag.diag( _lambdas );
	
	data_.full_tensor3_matrix_multiplication( core_diag, _u1, _u2, _u3 );
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
	hopm( data_ );
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hopm( const tensor3< I1, I2, I3, T >& data_ )
{
	
	//compute best rank-(R) approximation (Lathauwer et al., 2000b)
	tensor3< I1, I2, I3, T > approximated_data;
	reconstruction( approximated_data );
	double max_f_norm = data_.frobenius_norm();
	//std::cout << "frobenius norm original: " << max_f_norm << std::endl;
	
	double f_norm = approximated_data.frobenius_norm();
	double last_f_norm = f_norm;
	double improvement = max_f_norm - f_norm;
	double min_improvement = 0.0001;
	size_t i = 0;
	size_t max_iterations = 20;
	double lambda;
	
	//intialize u1-u3
	//hosvd_mode1( data_, _u1 ); inital guess not needed for u1 since it will be computed in the first optimization step
	hosvd_mode2( data_, _u2 );
	hosvd_mode3( data_, _u3 );
	
	//std::cout << "initial u2: " << std::endl << _u2 << std::endl;
	//std::cout << "initial u3: " << std::endl << _u3 << std::endl;
	
	//std::cout << " data: " << std::endl << data_ << std::endl;
	
	while( improvement > min_improvement && i < max_iterations )
	{
		
		//optimize u1
		optimize_mode1( data_, _u1, _u2, _u3);
		//std::cout << std::endl << " *** iteration: " << i << std::endl << " new u1: " << std::endl << _u1 << std::endl;

		//optimize u1
		optimize_mode2( data_, _u1, _u2, _u3);
		//std::cout << " new u2: " << std::endl << _u2 << std::endl;

		//optimize u1
		lambda = optimize_mode3( data_, _u1, _u2, _u3);
		//std::cout << " new u3: " << std::endl << _u3 <<  std::endl;
		
		
		set_u1(_u1); set_u2(_u2); set_u3(_u3);
		_lambdas.at(0) = lambda; //TODO: for all lambdas/ranks
		set_lambdas( _lambdas );
		
		reconstruction( approximated_data );
		f_norm = approximated_data.frobenius_norm();
		improvement = f_norm - last_f_norm;
		last_f_norm = f_norm;
						
		++i;
	}
	
	std::cout << "number of cp iterations: " << i << std::endl;
}

	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode1( const tensor3< I1, I2, I3, T >& data_, u1_type& U1_ ) const
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
VMML_TEMPLATE_CLASSNAME::hosvd_mode2( const tensor3< I1, I2, I3, T >& data_, u2_type& U2_ ) const
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
VMML_TEMPLATE_CLASSNAME::hosvd_mode3( const tensor3< I1, I2, I3, T >& data_, u3_type& U3_ ) const
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

	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode1( const tensor3< I1, I2, I3, T >& data_, u1_type& U1_optimized_, const u2_type& U2_, const u3_type& U3_ ) const
{	
	matrix< I1, I2*I3, T> unfolding; // -> u1
	data_.lateral_matricization( unfolding);
	
	matrix< I2*I3, R, T> u1_krp;
	u1_krp = U2_.khatri_rao_product( U3_ );	
	U1_optimized_.multiply( unfolding, u1_krp );
	
	//std::cout << "m_lateral " << std::endl << m_lateral << std::endl;
	//std::cout << "khatri-rao  " << std::endl << u1_krp << std::endl;
	//std::cout << "u1_ optimized " << std::endl << U1_optimized_ << std::endl;
	
	//normalize u with lambda (= norm)
	double lambda = U1_optimized_.frobenius_norm();
	U1_optimized_ *= (1 / lambda);
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode2( const tensor3< I1, I2, I3, T >& data_, const u1_type& U1_, u2_type& U2_optimized_, const u3_type& U3_ ) const
{
	matrix< I2, I1*I3, T> unfolding; // -> u2
	data_.frontal_matricization( unfolding);
	
	matrix< I1*I3, R, T> u2_krp;
	u2_krp = U1_.khatri_rao_product( U3_ );
	U2_optimized_.multiply( unfolding, u2_krp );
	
	
	//normalize u with lambda (= norm)
	double lambda = U2_optimized_.frobenius_norm();
	U2_optimized_ *= (1 / lambda);
}	


VMML_TEMPLATE_STRING
double  
VMML_TEMPLATE_CLASSNAME::optimize_mode3( const tensor3< I1, I2, I3, T >& data_, const u1_type& U1_, const u2_type& U2_,  u3_type& U3_optimized_ ) const
{
	matrix< I3, I1*I2, T> unfolding; //-> u3
	data_.horizontal_matricization( unfolding);
	
	matrix< I1*I2, R, T> u3_krp;
	u3_krp = U1_.khatri_rao_product( U2_ );
	U3_optimized_.multiply( unfolding, u3_krp );
	
	//normalize u with lambda (= norm)
	double lambda = U3_optimized_.frobenius_norm();
	U3_optimized_ *= (1 / lambda);
	
	return lambda;
}

	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::export_to( std::vector< T >& data_ ) const
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
	
	//TODO: iterate over lambdas
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_from( std::vector< T >& data_ )
{
	size_t i = 0; //iterator over data_
	
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
	
	//TODO: import lambdas
	
}	

	
	
	
		
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
		
	} // namespace vmml
	
#endif
		
