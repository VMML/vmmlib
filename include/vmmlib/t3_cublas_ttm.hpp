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
		
		//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
		template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
		static void full_tensor3_matrix_multiplication( const tensor3< J1, J2, J3, T >& t3_in_, const matrix< I1, J1, T >& U1, const matrix< I2, J2, T >& U2, const matrix< I3, J3, T >& U3, tensor3< I1, I2, I3, T >& t3_res_ );
		
		template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
		static void full_tensor3_matrix_kronecker_mult( const tensor3< J1, J2, J3, T >& t3_in_, const matrix< I1, J1, T >& U1, const matrix< I2, J2, T >& U2, const matrix< I3, J3, T >& U3, tensor3< I1, I2, I3, T >& t3_res_ );
		
		//tensor times matrix multiplication along different modes
		template< size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
		static void multiply_horizontal_bwd( const tensor3< J1, J2, J3, T >& t3_in_, const matrix< I3, J3, T >& in_slice_, tensor3< J1, J2, I3, T >& t3_res_ ); //output: tensor3< J1, J2, I3, T >  
		
		template< size_t I1, size_t J1, size_t J2, size_t J3, typename T > 
		static void multiply_lateral_bwd( const tensor3< J1, J2, J3, T >& t3_in_, const matrix< I1, J1, T >& in_slice_, tensor3< I1, J2, J3, T >& t3_res_ ); //output: tensor3< I1, J2, J3, T > 
		
		template< size_t I2, size_t J1, size_t J2, size_t J3, typename T > 
		static void multiply_frontal_bwd( const tensor3< J1, J2, J3, T >& t3_in_, const matrix< I2, J2, T >& in_slice_, tensor3< J1, I2, J3, T >& t3_res_ ); //output: tensor3< J1, I2, J3, T >
		
	protected:
		
		
	}; //end hosvd class
	
#define VMML_TEMPLATE_CLASSNAME     t3_cublas_ttm
	
	
	
	
	
	template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
	void
	VMML_TEMPLATE_CLASSNAME::full_tensor3_matrix_multiplication(  const tensor3< J1, J2, J3, T >& t3_in_, 
																const matrix< I1, J1, T >& U1, 
																const matrix< I2, J2, T >& U2, 
																const matrix< I3, J3, T >& U3,
																tensor3< I1, I2, I3, T >& t3_res_
																)
	{
		tensor3< I1, J2, J3, T> t3_result_1;
		tensor3< I1, I2, J3, T> t3_result_2;
		
		//backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
		multiply_lateral_bwd( t3_in_, U1, t3_result_1 );
		multiply_frontal_bwd( t3_result_1, U2, t3_result_2 );
		multiply_horizontal_bwd( t3_result_2, U3, t3_res_ );
	}
	
	template< size_t I1, size_t I2, size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
	void
	VMML_TEMPLATE_CLASSNAME::full_tensor3_matrix_kronecker_mult(  const tensor3< J1, J2, J3, T >& t3_in_, 
																const matrix< I1, J1, T >& U1, 
																const matrix< I2, J2, T >& U2, 
																const matrix< I3, J3, T >& U3,
																tensor3< I1, I2, I3, T >& t3_res_
																)
	{
		//TODO should use blas
		
		matrix< J1, J2*J3, T>* core_unfolded = new matrix< J1, J2*J3, T>;
		t3_in_.lateral_unfolding_bwd( *core_unfolded );
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
			t3_res_.set_column( i2, i3, res_unfolded->get_column(i));
		}
		
		delete core_unfolded;
		delete kron_prod;
		delete tmp1;
		delete res_unfolded;
	}
	
	
	
	
	//tensor matrix multiplications
	
	template< size_t I3, size_t J1, size_t J2, size_t J3, typename T > 
	void
	VMML_TEMPLATE_CLASSNAME::multiply_horizontal_bwd( const tensor3< J1, J2, J3, T >& t3_in_, 
													 const matrix< I3, J3, T >& in_slice_, 
													 tensor3< J1, J2, I3, T >& t3_res_ )
	{
		typedef matrix< J3, J2, T > slice_t;
		typedef matrix< I3, J2, T > slice_new_t;
		typedef matrix< J3, J2, T_cublas > cublas_slice_t;
		typedef matrix< I3, J2, T_cublas > cublas_slice_new_t;
		typedef matrix< I3, J3, T_cublas > cublas_in_slice_t;
		typedef cublas_dgemm< I3, J3, J2, T_cublas > cublas_t;
		
		slice_t* slice = new slice_t;
		slice_new_t* slice_new = new slice_new_t;
		
		cublas_slice_t* slice_blas = new cublas_slice_t;
		cublas_slice_new_t* slice_new_blas = new cublas_slice_new_t;
		cublas_in_slice_t* in_slice_blas = new cublas_in_slice_t;
		in_slice_blas->cast_from( in_slice_ );
		
		cublas_t* multiplier = new cublas_t;
		
		for ( size_t i1 = 0; i1 < J1; ++i1 )
		{
			t3_in_.get_horizontal_slice_bwd( i1, *slice );
			
			slice_blas->cast_from( *slice );
			
			multiplier->compute( *in_slice_blas, *slice_blas, *slice_new_blas );
			//slice_new->multiply( in_slice_, *slice );
			
			slice_new->cast_from( *slice_new_blas );
			
			t3_res_.set_horizontal_slice_bwd( i1, *slice_new );		
			
		}
				
		delete multiplier;	
		delete slice;
		delete slice_new;
		delete slice_blas;
		delete slice_new_blas;
		delete in_slice_blas;
	}
	
	
	template< size_t I1, size_t J1, size_t J2, size_t J3, typename T > 
	void
	VMML_TEMPLATE_CLASSNAME::multiply_lateral_bwd( const tensor3< J1, J2, J3, T >& t3_in_, 
												  const matrix< I1, J1, T >& in_slice_, 
												  tensor3< I1, J2, J3, T >& t3_res_ )
	{
		typedef matrix< J1, J3, T > slice_t;
		typedef matrix< I1, J3, T > slice_new_t;
		typedef matrix< J1, J3, T_cublas > cublas_slice_t;
		typedef matrix< I1, J3, T_cublas > cublas_slice_new_t;
		typedef matrix< I1, J1, T_cublas > cublas_in_slice_t;
		typedef cublas_dgemm< I1, J1, J3, T_cublas > cublas_t;
		
		slice_t* slice = new slice_t;
		slice_new_t* slice_new = new slice_new_t;
		
		cublas_slice_t* slice_blas = new cublas_slice_t;
		cublas_slice_new_t* slice_new_blas = new cublas_slice_new_t;
		cublas_in_slice_t* in_slice_blas = new cublas_in_slice_t;
		in_slice_blas->cast_from( in_slice_ );
		
		cublas_t* multiplier = new cublas_t;

		for ( size_t i2 = 0; i2 < J2; ++i2 )
		{
			t3_in_.get_lateral_slice_bwd( i2, *slice );
			
			slice_blas->cast_from( *slice );
			multiplier->compute( *in_slice_blas, *slice_blas, *slice_new_blas );
			slice_new->cast_from( *slice_new_blas );
			
			//slice_new->multiply( other_slice_, *slice );
			t3_res_.set_lateral_slice_bwd( i2, *slice_new );		
		}
		
		delete multiplier;	
		delete slice;
		delete slice_new;
		delete slice_blas;
		delete slice_new_blas;
		delete in_slice_blas;
	}
	
	
	
	template< size_t I2, size_t J1, size_t J2, size_t J3, typename T > 
	void
	VMML_TEMPLATE_CLASSNAME::multiply_frontal_bwd( const tensor3< J1, J2, J3, T >& t3_in_, 
												  const matrix< I2, J2, T >& in_slice_, 
												  tensor3< J1, I2, J3, T >& t3_res_  )
	{
		typedef matrix< J2, J1, T > slice_t;
		typedef matrix< I2, J1, T > slice_new_t;
		typedef matrix< J2, J1, T_cublas > cublas_slice_t;
		typedef matrix< I2, J1, T_cublas > cublas_slice_new_t;
		typedef matrix< I2, J2, T_cublas > cublas_in_slice_t;
		typedef cublas_dgemm< I2, J2, J1, T_cublas > cublas_t;
		
		slice_t* slice = new slice_t;
		slice_new_t* slice_new = new slice_new_t;
		
		cublas_slice_t* slice_blas = new cublas_slice_t;
		cublas_slice_new_t* slice_new_blas = new cublas_slice_new_t;
		cublas_in_slice_t* in_slice_blas = new cublas_in_slice_t;
		in_slice_blas->cast_from( in_slice_ );
		
		cublas_t* multiplier = new cublas_t;
		
		for ( size_t i3 = 0; i3 < J3; ++i3 )
		{
			t3_in_.get_frontal_slice_bwd( i3, *slice );
			
			slice_blas->cast_from( *slice );
			multiplier->compute( *in_slice_blas, *slice_blas, *slice_new_blas );
			slice_new->cast_from( *slice_new_blas );
			
			//slice_new->multiply( other_slice_, *slice );
			t3_res_.set_frontal_slice_bwd( i3, *slice_new );	
		}
		
		delete multiplier;	
		delete slice;
		delete slice_new;
		delete slice_blas;
		delete slice_new_blas;
		delete in_slice_blas;
	}
	
	
	
	
#undef VMML_TEMPLATE_CLASSNAME
	
}//end vmml namespace

#endif

