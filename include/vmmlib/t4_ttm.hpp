/* 
 * VMMLib - Tensor Classes
 *  
 * @author Rafael Ballester
 *
 * Tensor times matrix multiplication for tensor4 (t4)
 * using BLAS
 * see e.g.:
 * - Bader & Kolda, 2006: Algorithm 862: Matlab tensor classes for fast algorithm prototyping. ACM Transactions on Mathematical Software.
 * 
 */

#ifndef __VMML__T4_TTM__HPP__
#define __VMML__T4_TTM__HPP__

#include <vmmlib/tensor4.hpp>
#include <vmmlib/blas_dgemm.hpp>
#ifdef VMMLIB_USE_OPENMP
#  include <omp.h>
#endif

namespace vmml
{
	
	class t4_ttm
	{
	public:    

		template< size_t I1, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
		static void mode1_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I1, J1, T >& in_slice_, tensor4< I1, J2, J3, J4, T >& t4_res_ );
        
        template< size_t I2, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
		static void mode2_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I2, J2, T >& in_slice_, tensor4< J1, I2, J3, J4, T >& t4_res_ );
                
        template< size_t I3, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
		static void mode3_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I3, J3, T >& in_slice_, tensor4< J1, J2, I3, J4, T >& t4_res_ );
                        
        template< size_t I4, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
		static void mode4_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I4, J4, T >& in_slice_, tensor4< J1, J2, J3, I4, T >& t4_res_ );
		
	protected:
			
	};
	
#define VMML_TEMPLATE_CLASSNAME     t4_ttm

	template< size_t I1, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
    void
    VMML_TEMPLATE_CLASSNAME::mode1_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I1, J1, T >& in_slice_, tensor4< I1, J2, J3, J4, T >& t4_res_ ) {
        
    }
    
    template< size_t I2, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
	void
    VMML_TEMPLATE_CLASSNAME::mode2_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I2, J2, T >& in_slice_, tensor4< J1, I2, J3, J4, T >& t4_res_ ) {
        
    }
    
    template< size_t I3, size_t J1, size_t J2, size_t J3, size_t J4, typename T  > 
    void
    VMML_TEMPLATE_CLASSNAME::mode3_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I3, J3, T >& in_slice_, tensor4< J1, J2, I3, J4, T >& t4_res_ ) {
        
    }
    
    template< size_t I4, size_t J1, size_t J2, size_t J3, size_t J4, typename T  >
    void
    VMML_TEMPLATE_CLASSNAME::mode4_multiply_fwd( const tensor4< J1, J2, J3, J4, T >& t4_in_, const matrix< I4, J4, T >& in_slice_, tensor4< J1, J2, J3, I4, T >& t4_res_ ) {
        
    }

#undef VMML_TEMPLATE_CLASSNAME
	
}//end vmml namespace

#endif
