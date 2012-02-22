/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * Tensor times matrix multiplication for tensor3 (t3)
 * see e.g.:
 * - Bader & Kolda, 2006: Algorithm 862: Matlab tensor classes for fast algorithm prototyping. ACM Transactions on Mathematical Software.
 * 
 */

#ifndef __VMML__T3_TTM__HPP__
#define __VMML__T3_TTM__HPP__

#include <vmmlib/tensor3.hpp>
#include <vmmlib/blas_dgemm.hpp>

//TODO
enum dgemm_method {
	blas_e,
	cublas_e
}; 


namespace vmml
{
	
	class t3_ttm
	{
	public:    
		
		typedef float T_blas;
	
		//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
		template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
		static void full_tensor3_matrix_multiplication( const tensor3< J1, J2, J3, T >& core_, const matrix< I1, J1, T >& U1, const matrix< I2, J2, T >& U2, const matrix< I3, J3, T >& U3, tensor3< I1, I2, I3, T >& reco_ );
	
		template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
		static void full_tensor3_matrix_kronecker_mult( const tensor3< J1, J2, J3, T >& core_, const matrix< I1, J1, T >& U1, const matrix< I2, J2, T >& U2, const matrix< I3, J3, T >& U3, tensor3< I1, I2, I3, T >& reco_ );
		
		//tensor times matrix multiplication along different modes
		template< size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
		static void multiply_horizontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I3, J3, T >& other_slice_, tensor3< J1, J2, I3, T >& t3_res_ ); //output: tensor3< J1, J2, I3, T >  
		
		template< size_t I1, size_t J1, size_t J2, size_t J3, typename T > 
		static void multiply_lateral_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I1, J1, T >& other_slice_, tensor3< I1, J2, J3, T >& t3_res_ ); //output: tensor3< I1, J2, J3, T > 
		
		template< size_t I2, size_t J1, size_t J2, size_t J3, typename T > 
		static void multiply_frontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I2, J2, T >& other_slice_, tensor3< J1, I2, J3, T >& t3_res_ ); //output: tensor3< J1, I2, J3, T >
		
	protected:
		
			
	}; //end hosvd class
	
#define VMML_TEMPLATE_CLASSNAME     t3_ttm
	
	


	
template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
void
VMML_TEMPLATE_CLASSNAME::full_tensor3_matrix_multiplication(  const tensor3< J1, J2, J3, T >& core_, 
															const matrix< I1, J1, T >& U1, 
															const matrix< I2, J2, T >& U2, 
															const matrix< I3, J3, T >& U3,
															tensor3< I1, I2, I3, T >& reco_
															)
{
	tensor3< I1, J2, J3, T> t3_result_1;
	tensor3< I1, I2, J3, T> t3_result_2;
	
	//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
	multiply_lateral_bwd( core_, U1, t3_result_1 );
	multiply_frontal_bwd( t3_result_1, U2, t3_result_2 );
	multiply_horizontal_bwd( t3_result_2, U3, reco_ );
}

template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
void
VMML_TEMPLATE_CLASSNAME::full_tensor3_matrix_kronecker_mult(  const tensor3< J1, J2, J3, T >& core_, 
															const matrix< I1, J1, T >& U1, 
															const matrix< I2, J2, T >& U2, 
															const matrix< I3, J3, T >& U3,
															tensor3< I1, I2, I3, T >& reco_
															)
{
	//TODO should use blas
	
	matrix< J1, J2*J3, T>* core_unfolded = new matrix< J1, J2*J3, T>;
	core_.lateral_unfolding_bwd( *core_unfolded );
	matrix< I1, J2*J3, T>* tmp1 = new matrix< I1, J2*J3, T>;
	tmp1->multiply( U1, *core_unfolded );
	
	matrix< I2*I3, J2*J3, T>* kron_prod = new matrix< I2*I3, J2*J3, T>;
	U2.kronecker_product( U3, *kron_prod );
	
	matrix< I1, I2*I3, T>* res_unfolded = new matrix< I1, I2*I3, T>;
	res_unfolded->multiply( *tmp1, transpose(*kron_prod) );
	
	//std::cout << "reco2 result (matricized): " << std::endl << *res_unfolded << std::endl;
	
	//set this from res_unfolded
	size_t i2 = 0;
	for( size_t i = 0; i < (I2*I3); ++i, ++i2 )
	{
		if (i2 >= I2) {
			i2 = 0;
		}
		
		size_t i3 = i % I3;
		reco_.set_column( i2, i3, res_unfolded->get_column(i));
	}
	
	delete core_unfolded;
	delete kron_prod;
	delete tmp1;
	delete res_unfolded;
}

	
	
	
//tensor matrix multiplications

template< size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
void
VMML_TEMPLATE_CLASSNAME::multiply_horizontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I3, J3, T >& other_slice_, tensor3< J1, J2, I3, T >& t3_res_ )
{
	matrix< J3, J2, T >* slice = new matrix< J3, J2, T >;
	matrix< I3, J2, T >* slice_new = new matrix< I3, J2, T >;
	
	blas_dgemm< I3, J3, J2, T_blas >* blas_dgemm1 = new blas_dgemm< I3, J3, J2, T_blas >;
	matrix< J3, J2, T_blas >* slice_blas = new matrix< J3, J2, T_blas >;
	matrix< I3, J2, T_blas >* slice_new_blas = new matrix< I3, J2, T_blas >;
	matrix< I3, J3, T_blas >* other_slice_blas = new matrix< I3, J3, T_blas >;
	other_slice_blas->cast_from( other_slice_ );
	
	for ( size_t i1 = 0; i1 < J1; ++i1 )
	{
		other.get_horizontal_slice_bwd( i1, *slice );
		
		slice_blas->cast_from( *slice );
		blas_dgemm1->compute( *other_slice_blas, *slice_blas, *slice_new_blas );
		slice_new->cast_from( *slice_new_blas );
		//slice_new->multiply( other_slice_, *slice );
		
		t3_res_.set_horizontal_slice_bwd( i1, *slice_new );		
	}
	
	delete blas_dgemm1;	
	delete slice;
	delete slice_new;
	delete slice_blas;
	delete slice_new_blas;
	delete other_slice_blas;
}


template< size_t I1, size_t J1, size_t J2, size_t J3, typename T > 
void
VMML_TEMPLATE_CLASSNAME::multiply_lateral_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I1, J1, T >& other_slice_, tensor3< I1, J2, J3, T >& t3_res_ )
{
	matrix< J1, J3, T>* slice = new matrix< J1, J3, T>();
	matrix< I1, J3, T>* slice_new = new matrix< I1, J3, T>();
	
	blas_dgemm< I1, J1, J3, T_blas >* blas_dgemm1 = new blas_dgemm< I1, J1, J3, T_blas >;
	matrix< J1, J3, T_blas >* slice_blas = new matrix< J1, J3, T_blas >;
	matrix< I1, J3, T_blas >* slice_new_blas = new matrix< I1, J3, T_blas >;
	matrix< I1, J1, T_blas >* other_slice_blas = new matrix< I1, J1, T_blas >;
	other_slice_blas->cast_from( other_slice_ );
	
	for ( size_t i2 = 0; i2 < J2; ++i2 )
	{
		other.get_lateral_slice_bwd( i2, *slice );
		
		slice_blas->cast_from( *slice );
		blas_dgemm1->compute( *other_slice_blas, *slice_blas, *slice_new_blas );
		slice_new->cast_from( *slice_new_blas );
		
		//slice_new->multiply( other_slice_, *slice );
		t3_res_.set_lateral_slice_bwd( i2, *slice_new );		
	}
	delete blas_dgemm1;	
	delete slice;
	delete slice_new;
	delete slice_blas;
	delete slice_new_blas;
	delete other_slice_blas;
}



template< size_t I2, size_t J1, size_t J2, size_t J3, typename T > 
void
VMML_TEMPLATE_CLASSNAME::multiply_frontal_bwd( const tensor3< J1, J2, J3, T >& other, const matrix< I2, J2, T >& other_slice_, tensor3< J1, I2, J3, T >& t3_res_  )
{
	matrix< J2, J1, T>* slice = new matrix< J2, J1, T>();
	matrix< I2, J1, T>* slice_new = new matrix< I2, J1, T>();
	
	blas_dgemm< I2, J2, J1, T_blas >* blas_dgemm1 = new blas_dgemm< I2, J2, J1, T_blas >;
	matrix< J2, J1, T_blas >* slice_blas = new matrix< J2, J1, T_blas >;
	matrix< I2, J1, T_blas >* slice_new_blas = new matrix< I2, J1, T_blas >;
	matrix< I2, J2, T_blas >* other_slice_blas = new matrix< I2, J2, T_blas >;
	other_slice_blas->cast_from( other_slice_ );
	
	for ( size_t i3 = 0; i3 < J3; ++i3 )
	{
		other.get_frontal_slice_bwd( i3, *slice );
		
		slice_blas->cast_from( *slice );
		blas_dgemm1->compute( *other_slice_blas, *slice_blas, *slice_new_blas );
		slice_new->cast_from( *slice_new_blas );
		
		//slice_new->multiply( other_slice_, *slice );
		t3_res_.set_frontal_slice_bwd( i3, *slice_new );		
	}
	delete blas_dgemm1;	
	delete slice;
	delete slice_new;
	delete slice_blas;
	delete slice_new_blas;
	delete other_slice_blas;
}

	
	

#undef VMML_TEMPLATE_CLASSNAME
	
}//end vmml namespace

#endif

