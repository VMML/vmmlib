/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * The tucker3 tensor class is consists of the same components (core tensor, basis matrices u1-u3) as the tucker3 model described in:
 * - Tucker,  “Some mathematical notes on three-mode factor analysis”, Psychometrika, vol. 31, no. 3, pp. 279–311., 1966 Sep.
 * - De Lathauwer L., De Moor B., Vandewalle J., ``A multilinear singular value decomposition'', 
 * SIAM J. Matrix Anal. Appl., vol. 21, no. 4, Apr. 2000, pp. 1253-1278.
 * - De Lathauwer L., De Moor B., Vandewalle J., ``On the Best rank-1 and Rank-$(R_1, R_2, ..., R_N)$ Approximation and Applications of Higher-Order Tensors'', 
 * SIAM J. Matrix Anal. Appl., vol. 21, no. 4, Apr. 2000, pp. 1324-1342.
 * - T. G. Kolda and B. W. Bader. Tensor Decompositions and Applications. 
 * SIAM Review, Volume 51, Number 3, Pages 455-500, September 2009.
 * 
 */

#ifndef __VMML__T3_HOPM__HPP__
#define __VMML__T3_HOPM__HPP__

#include <vmmlib/tensor3.hpp>

namespace vmml
{
	
	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T_value = float, typename T_coeff = double >
	class hopm
	{
	public:    
		
		typedef float T_internal;	
		
		
	};
	
	
}//end vmml namespace

#endif