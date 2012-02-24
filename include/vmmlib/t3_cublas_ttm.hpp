/* 
 * VMMLib - CUDA Tensor Classes
 *  
 * @author Susanne Suter
 *
 * Tensor times matrix multiplication for tensor3 (t3)
 * using CUBLAS 
 * see e.g.:
 * - Bader & Kolda, 2006: Algorithm 862: Matlab tensor classes for fast algorithm prototyping. ACM Transactions on Mathematical Software.
 * 
 */

#ifndef __VMML__T3_CUBLAS_TTM__HPP__
#define __VMML__T3_CUBLAS_TTM__HPP__

#include <vmmlib/tensor3.hpp>
#include <vmmlib/cublas_dgemm.cu>


namespace vmml
{
	
	class t3_cublas_ttm
	{
	public:    
		
		typedef float T_cublas;
		
		//floating point versions (without type casting
		//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
		template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3 > 
		static void full_tensor3_matrix_multiplication( const tensor3< J1, J2, J3, T_cublas >& t3_in_, const matrix< I1, J1, T_cublas >& U1, const matrix< I2, J2, T_cublas >& U2, const matrix< I3, J3, T_cublas >& U3, tensor3< I1, I2, I3, T_cublas >& t3_res_ );
		
		//tensor times matrix multiplication along different modes
		template< size_t I3, size_t J1, size_t J2, size_t J3 > 
		static void multiply_horizontal_bwd( const tensor3< J1, J2, J3, T_cublas >& t3_in_, const matrix< I3, J3, T_cublas >& in_slice_, tensor3< J1, J2, I3, T_cublas >& t3_res_ ); //output: tensor3< J1, J2, I3, T >  
		
		template< size_t I1, size_t J1, size_t J2, size_t J3  > 
		static void multiply_lateral_bwd( const tensor3< J1, J2, J3, T_cublas >& t3_in_, const matrix< I1, J1, T_cublas >& in_slice_, tensor3< I1, J2, J3, T_cublas >& t3_res_ ); //output: tensor3< I1, J2, J3, T > 
		
		template< size_t I2, size_t J1, size_t J2, size_t J3 > 
		static void multiply_frontal_bwd( const tensor3< J1, J2, J3, T_cublas >& t3_in_, const matrix< I2, J2, T_cublas >& in_slice_, tensor3< J1, I2, J3, T_cublas >& t3_res_ ); //output: tensor3< J1, I2, J3, T >

	protected:
		
		
	}; //end cublas ttm class
	
#define VMML_TEMPLATE_CLASSNAME     t3_cublas_ttm
	
	
	template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3 > 
	void
	VMML_TEMPLATE_CLASSNAME::full_tensor3_matrix_multiplication(  const tensor3< J1, J2, J3, T_cublas >& t3_in_, 
																const matrix< I1, J1, T_cublas >& U1, 
																const matrix< I2, J2, T_cublas >& U2, 
																const matrix< I3, J3, T_cublas >& U3,
																tensor3< I1, I2, I3, T_cublas >& t3_res_
																)
	{
		tensor3< I1, J2, J3, T_cublas > t3_result_1;
		tensor3< I1, I2, J3, T_cublas > t3_result_2;
		
		//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
		multiply_lateral_bwd( t3_in_, U1, t3_result_1 );
		multiply_frontal_bwd( t3_result_1, U2, t3_result_2 );
		multiply_horizontal_bwd( t3_result_2, U3, t3_res_ );
	}
	
	
	//tensor matrix multiplications
	
	template< size_t I3, size_t J1, size_t J2, size_t J3 > 
	void
	VMML_TEMPLATE_CLASSNAME::multiply_horizontal_bwd( const tensor3< J1, J2, J3, T_cublas >& t3_in_, 
													 const matrix< I3, J3, T_cublas >& in_slice_, 
													 tensor3< J1, J2, I3, T_cublas >& t3_res_ )
	{
#if 1
		typedef matrix< J3, J1*J2, T_cublas > unfolding_t; 
		typedef matrix< I3, J1*J2, T_cublas > unfolding_res_t; 
		typedef cublas_dgemm< I3, J3, J1 * J2, T_cublas > cublas_t;
		
		unfolding_t* unfolding = new unfolding_t;
		t3_in_.horizontal_unfolding_bwd( *unfolding );
		unfolding_res_t* unfolding_2 = new unfolding_res_t;
		cublas_t* multiplier = new cublas_t;
		multiplier->compute( in_slice_, *unfolding, *unfolding_2 );
		
		t3_res_.horizontal_folding_bwd( *unfolding_2 );
		
		delete unfolding;
		delete unfolding_2;
#else
		typedef matrix< J3, J2, T_cublas > slice_t;
		typedef matrix< I3, J2, T_cublas > slice_new_t;
		typedef cublas_dgemm< I3, J3, J2, T_cublas > cublas_t;
		
		slice_t* slice = new slice_t;
		slice_new_t* slice_new = new slice_new_t;
		
		cublas_t* multiplier = new cublas_t;
		
		for ( size_t i1 = 0; i1 < J1; ++i1 )
		{
			t3_in_.get_horizontal_slice_bwd( i1, *slice );
			
			multiplier->compute( in_slice_, *slice, *slice_new );
			
			t3_res_.set_horizontal_slice_bwd( i1, *slice_new );		
		}
		
		delete multiplier;	
		delete slice;
		delete slice_new;
#endif
	}
	
	
	template< size_t I1, size_t J1, size_t J2, size_t J3 > 
	void
	VMML_TEMPLATE_CLASSNAME::multiply_lateral_bwd( const tensor3< J1, J2, J3, T_cublas >& t3_in_, 
												  const matrix< I1, J1, T_cublas >& in_slice_, 
												  tensor3< I1, J2, J3, T_cublas >& t3_res_ )
	{
#if 1
		typedef matrix< J1, J3*J2, T_cublas > unfolding_t; 
		typedef matrix< I1, J3*J2, T_cublas > unfolding_res_t; 
		typedef cublas_dgemm< I1, J1, J3 * J2, T_cublas > cublas_t;
		
		unfolding_t* unfolding = new unfolding_t;
		t3_in_.lateral_unfolding_bwd( *unfolding );
		unfolding_res_t* unfolding_2 = new unfolding_res_t;
		cublas_t* multiplier = new cublas_t;
		multiplier->compute( in_slice_, *unfolding, *unfolding_2 );
		
		t3_res_.lateral_folding_bwd( *unfolding_2 );
		
		delete unfolding;
		delete unfolding_2;
#else
		typedef matrix< J1, J3, T_cublas > slice_t;
		typedef matrix< I1, J3, T_cublas > slice_new_t;
		typedef cublas_dgemm< I1, J1, J3, T_cublas > cublas_t;
		
		slice_t* slice = new slice_t;
		slice_new_t* slice_new = new slice_new_t;
		
		cublas_t* multiplier = new cublas_t;
		
		for ( size_t i2 = 0; i2 < J2; ++i2 )
		{
			t3_in_.get_lateral_slice_bwd( i2, *slice );
			
			multiplier->compute( in_slice_, *slice, *slice_new );
			
			t3_res_.set_lateral_slice_bwd( i2, *slice_new );		
		}
		
		delete multiplier;	
		delete slice;
		delete slice_new;
#endif
	}
	
	
	
	template< size_t I2, size_t J1, size_t J2, size_t J3 > 
	void
	VMML_TEMPLATE_CLASSNAME::multiply_frontal_bwd( const tensor3< J1, J2, J3, T_cublas >& t3_in_, 
												  const matrix< I2, J2, T_cublas >& in_slice_, 
												  tensor3< J1, I2, J3, T_cublas >& t3_res_  )
	{
#if 1
		typedef matrix< J2, J1*J3, T_cublas > unfolding_t; 
		typedef matrix< I2, J1*J3, T_cublas > unfolding_res_t; 
		typedef cublas_dgemm< I2, J2, J1 * J3, T_cublas > cublas_t;
		
		unfolding_t* unfolding = new unfolding_t;
		t3_in_.frontal_unfolding_bwd( *unfolding );
		unfolding_res_t* unfolding_2 = new unfolding_res_t;
		cublas_t* multiplier = new cublas_t;
		multiplier->compute( in_slice_, *unfolding, *unfolding_2 );
		
		t3_res_.frontal_folding_bwd( *unfolding_2 );
		
		delete unfolding;
		delete unfolding_2;
#else
		typedef matrix< J2, J1, T_cublas > slice_t;
		typedef matrix< I2, J1, T_cublas > slice_new_t;
		typedef cublas_dgemm< I2, J2, J1, T_cublas > cublas_t;
		
		slice_t* slice = new slice_t;
		slice_new_t* slice_new = new slice_new_t;
		
		cublas_t* multiplier = new cublas_t;
		
		for ( size_t i3 = 0; i3 < J3; ++i3 )
		{
			t3_in_.get_frontal_slice_bwd( i3, *slice );
			
			multiplier->compute( in_slice_, *slice, *slice_new );
			
			t3_res_.set_frontal_slice_bwd( i3, *slice_new );	
		}
		
		delete multiplier;	
		delete slice;
		delete slice_new;
#endif
	}
	
	
#undef VMML_TEMPLATE_CLASSNAME
	
}//end vmml namespace

#endif

