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

//TODO allocate data with new

namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, size_t R, typename T_value = float, typename T_coeff = float >
	class cp3_tensor
	{
	public:    		
		typedef tensor3< I1, I2, I3, T_value > t3_type;
		typedef typename t3_type::iterator t3_iterator;
		typedef typename t3_type::const_iterator t3_const_iterator;
		
		typedef tensor3< I1, I2, I3, T_coeff > t3_coeff_type;
		typedef typename t3_coeff_type::iterator t3_coeff_iterator;
		typedef typename t3_coeff_type::const_iterator t3_coeff_const_iterator;
		
		typedef matrix< I1, R, T_coeff > u1_type;
		typedef typename u1_type::iterator u1_iterator;
		typedef typename u1_type::const_iterator u1_const_iterator;
		
		typedef matrix< I2, R, T_coeff > u2_type;
		typedef typename u2_type::iterator u2_iterator;
		typedef typename u2_type::const_iterator u2_const_iterator;
		
		typedef matrix< I3, R, T_coeff > u3_type;
		typedef typename u3_type::iterator u3_iterator;
		typedef typename u3_type::const_iterator u3_const_iterator;
		
		typedef matrix< I1, I2*I3, T_coeff > mode1_matricization_type;
		typedef matrix< I2, I1*I3, T_coeff > mode2_matricization_type;
		typedef matrix< I3, I1*I2, T_coeff > mode3_matricization_type;		

		cp3_tensor(  u1_type& U1, u2_type& U2, u3_type& U3, vector< R, T_coeff >& lambdas_ );
		cp3_tensor();
		~cp3_tensor();
		
		void set_lambdas( const vector< R, T_coeff >& lambdas_ )  { *_lambdas = vector< R, T_coeff >( lambdas_); } ;
		void set_u1( u1_type& U1 ) { *_u1 = U1; } ;
		void set_u2( u2_type& U2 ) { *_u2 = U2; } ;
		void set_u3( u3_type& U3 ) { *_u3 = U3; } ;
		
		void get_lambdas( vector< R, T_coeff >& data_ ) const { data_  = *_lambdas; } ;
		void get_u1( u1_type& U1 ) const { U1 = *_u1; } ;
		void get_u2( u2_type& U2 ) const { U2 = *_u2; } ;
		void get_u3( u3_type& U3 ) const { U3 = *_u3; } ;
		
		void export_to( std::vector< T_coeff >& data_ ) const;
		void import_from( std::vector< T_coeff >& data_ );	
		
		void decompose( const t3_type& data_ ); 
		void reconstruct( t3_type& data_ ) const;
		
		void cp_als( const t3_type& data_ );
		
		//higher-order power method (lathauwer et al., 2000b)
		void hopm( const t3_type& data_ );
		
		void hosvd_mode1( const t3_type& data_, u1_type& U1_ ) const;
		void hosvd_mode2( const t3_type& data_, u2_type& U2_ ) const;
		void hosvd_mode3( const t3_type& data_, u3_type& U3_ ) const;
		
		void optimize_mode1( const t3_type& data_, u1_type& U1_optimized_, const u2_type& U2_, const u3_type& U3_ ) const;
		void optimize_mode2( const t3_type& data_, const u1_type& U1_, u2_type& U2_optimized_, const u3_type& U3_ ) const;		
		float_t optimize_mode3( const t3_type& data_, const u1_type& U1_, const u2_type& U2_, u3_type& U3_optimized_ ) const;
		
	protected:
		cp3_tensor( const cp3_tensor< R, I1, I1, I1, T_value, T_coeff >& other ) {};
		cp3_tensor< R, I1, I1, I1, T_value, T_coeff > operator=( const cp3_tensor< R, I1, I1, I1, T_value, T_coeff >& other ) { return *this; };
		
	private:
		vector< R, T_coeff >* _lambdas ;
		u1_type* _u1 ;
		u2_type* _u2 ;
		u3_type* _u3 ;
		
	}; // class cp3_tensor
	
	
#define VMML_TEMPLATE_STRING    	template< size_t I1, size_t I2, size_t I3, size_t R, typename T_value, typename T_coeff >
#define VMML_TEMPLATE_CLASSNAME     cp3_tensor< I1, I2, I3, R, T_value, T_coeff >
	
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::cp3_tensor( u1_type& U1, u2_type& U2, u3_type& U3, vector< R, T_coeff >& lambdas_ )
{
	set_lambdas(lambdas_);
	set_u1( U1);
	set_u2( U2);
	set_u3( U3);
}

VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::cp3_tensor()
{
	_lambdas = new vector< R, T_coeff>(); _lambdas->core();
	_u1 = new u1_type(); _u1->zero();
	_u2 = new u2_type(); _u2->zero();
	_u3 = new u3_type(); _u3->zero();
}
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::~cp3_tensor()
{
	delete _u1;
	delete _u2;
	delete _u3;
	delete _lambdas;
}
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruct( t3_type& data_ ) const
{
	tensor3< R, R, R, T_coeff > core_diag;
	core_diag.diag( *_lambdas );
	
	t3_coeff_type data;
	data.cast_from( data_ );
	data.full_tensor3_matrix_multiplication( core_diag, *_u1, *_u2, *_u3 );
	if( (sizeof(T_value) == 1) || (sizeof(T_value) == 2) ){
		data_.float_t_to_uint_t( data );
	} else {
		data_.cast_from( data );
	}
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::decompose( const t3_type& data_ )
{
	cp_als( data_ );
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::cp_als( const t3_type& data_ )
{
	hopm( data_ );
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hopm( const t3_type& data_ )
{
	
	//compute best rank-(R) approximation (Lathauwer et al., 2000b)
	t3_type approximated_data;
	reconstruct( approximated_data );
	float_t max_f_norm = data_.frobenius_norm();
	//std::cout << "frobenius norm original: " << max_f_norm << std::endl;
	
	float_t f_norm = approximated_data.frobenius_norm();
	float_t last_f_norm = f_norm;
	float_t improvement = max_f_norm - f_norm;
	float_t min_improvement = 0.0001;
	size_t i = 0;
	size_t max_iterations = 20;
	float_t lambda;
	
	//intialize u1-u3
	//hosvd_mode1( data_, _u1 ); inital guess not needed for u1 since it will be computed in the first optimization step
	hosvd_mode2( data_, *_u2 );
	hosvd_mode3( data_, *_u3 );
	
	//std::cout << "initial u2: " << std::endl << _u2 << std::endl;
	//std::cout << "initial u3: " << std::endl << _u3 << std::endl;
	
	//std::cout << " data: " << std::endl << data_ << std::endl;
	
	while( improvement > min_improvement && i < max_iterations )
	{
		
		//optimize u1
		optimize_mode1( data_, *_u1, *_u2, *_u3);
		//std::cout << std::endl << " *** iteration: " << i << std::endl << " new u1: " << std::endl << _u1 << std::endl;

		//optimize u1
		optimize_mode2( data_, *_u1, *_u2, *_u3);
		//std::cout << " new u2: " << std::endl << _u2 << std::endl;

		//optimize u1
		lambda = optimize_mode3( data_, *_u1, *_u2, *_u3);
		//std::cout << " new u3: " << std::endl << _u3 <<  std::endl;
		
		
		set_u1( *_u1 ); set_u2( *_u2 ); set_u3( *_u3 );
		_lambdas->at(0) = lambda; //TODO: for all lambdas/ranks
		set_lambdas( *_lambdas );
		
		reconstruct( approximated_data );
		f_norm = approximated_data.frobenius_norm();
		improvement = f_norm - last_f_norm;
		last_f_norm = f_norm;
						
		++i;
	}
	
	std::cout << "number of cp iterations: " << i << std::endl;
}

	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode1( const t3_type& data_, u1_type& U1_ ) const
{
	t3_coeff_type data;
	data.cast_from( data_ );
	mode1_matricization_type u; // -> u1
	data.lateral_matricization( u);
		
	vector< I2*I3, T_coeff > lambdas;
	lapack_svd< I1, I2*I3, T_coeff > svd;
	if( svd.compute_and_overwrite_input( u, lambdas ))
		u.get_sub_matrix( U1_ );
	else 
		U1_.zero();
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode2( const t3_type& data_, u2_type& U2_ ) const
{
	t3_coeff_type data;
	data.cast_from( data_ );
	mode2_matricization_type u; // -> u2
	data.frontal_matricization_bwd( u );
	
	vector< I1*I3, T_coeff > lambdas;
	lapack_svd< I2, I1*I3, T_coeff > svd;
	if( svd.compute_and_overwrite_input( u, lambdas ))
		u.get_sub_matrix( U2_ );
	else 
		U2_.zero();
}



VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::hosvd_mode3( const t3_type& data_, u3_type& U3_ ) const
{
	t3_coeff_type data;
	data.cast_from( data_ );
	mode3_matricization_type u; //-> u3
	data.horizontal_matricization_bwd( u );
	
	vector< I1*I2, T_coeff > lambdas;
	lapack_svd< I3, I1*I2, T_coeff > svd;
	if( svd.compute_and_overwrite_input( u, lambdas ))
		u.get_sub_matrix( U3_ );
	else 
		U3_.zero();
}

	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode1( const t3_type& data_, u1_type& U1_optimized_, const u2_type& U2_, const u3_type& U3_ ) const
{	
	t3_coeff_type data;
	data.cast_from( data_ );
	mode1_matricization_type unfolding; // -> u1
	data.lateral_matricization_bwd( unfolding );
	
	matrix< I2*I3, R, T_coeff> u1_krp;
	u1_krp = U2_.khatri_rao_product( U3_ );	
	U1_optimized_.multiply( unfolding, u1_krp );
	
	//std::cout << "m_lateral " << std::endl << m_lateral << std::endl;
	//std::cout << "khatri-rao  " << std::endl << u1_krp << std::endl;
	//std::cout << "u1_ optimized " << std::endl << U1_optimized_ << std::endl;
	
	//normalize u with lambda (= norm)
	float_t lambda = U1_optimized_.frobenius_norm();
	U1_optimized_ *= (1 / lambda);
}


VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::optimize_mode2( const t3_type& data_, const u1_type& U1_, u2_type& U2_optimized_, const u3_type& U3_ ) const
{
	t3_coeff_type data;
	data.cast_from( data_ );
	mode2_matricization_type unfolding; // -> u2
	data.frontal_matricization_bwd( unfolding );
	
	matrix< I1*I3, R, T_coeff> u2_krp;
	u2_krp = U1_.khatri_rao_product( U3_ );
	U2_optimized_.multiply( unfolding, u2_krp );
	
	
	//normalize u with lambda (= norm)
	float_t lambda = U2_optimized_.frobenius_norm();
	U2_optimized_ *= (1 / lambda);
}	


VMML_TEMPLATE_STRING
float_t  
VMML_TEMPLATE_CLASSNAME::optimize_mode3( const t3_type& data_, const u1_type& U1_, const u2_type& U2_,  u3_type& U3_optimized_ ) const
{
	t3_coeff_type data;
	data.cast_from( data_ );
	mode3_matricization_type unfolding; //-> u3
	data.horizontal_matricization_bwd( unfolding);
	
	matrix< I1*I2, R, T_coeff> u3_krp;
	u3_krp = U1_.khatri_rao_product( U2_ );
	U3_optimized_.multiply( unfolding, u3_krp );
	
	//normalize u with lambda (= norm)
	float_t lambda = U3_optimized_.frobenius_norm();
	U3_optimized_ *= (1 / lambda);
	
	return lambda;
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
	
	//TODO: iterate over lambdas
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_from( std::vector< T_coeff >& data_ )
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
		
