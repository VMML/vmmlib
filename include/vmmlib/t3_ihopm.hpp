/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * iHOPM stands for incremental higher-order power method.
 * in other words, it is an incremental rank-r CP-ALS
 * 
 * CP stands for Candecomp/Parafac (1970); ALS for alternating least squares algorithm
 * references:
 * - Carroll & Chang, 1970: Analysis of Individual Differences in Multidimensional Scaling via an N-way generalization of ``Eckart--Young'' decompositions, Psychometrika.
 * - Harshman, 1970: Foundations of the PARAFAC procedure: Models and conditions for an 'explanatory' multi-modal factor analysis, UCLA Working Papers in Phonetics.
 * - De Lathauwer, De Moor, Vandewalle, 2000: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * 
 * incremental rank-r approximation:
 * - Zang & Golub, 2001: Rank-one approximation to higher order tensors, SIAM J. Matrix Anal. Appl.
 */



#ifndef __VMML__T3_IHOPM__HPP__
#define __VMML__T3_IHOPM__HPP__

#include <vmmlib/t3_hopm.hpp>

namespace vmml
{
	
	template< size_t R, size_t R_incr, size_t I1, size_t I2, size_t I3, typename T_val = float, typename T_res = double >
	class t3_ihopm
	{
	public:    
		
		typedef tensor3< I1, I2, I3, T_val > t3_type;
		typedef tensor3< I1, I2, I3, T_res > t3_res_type;
		
		typedef vector< R*R_incr, T_val > lambda_type;
		typedef vector< R*R_incr, T_res > lambda_res_type;
		
		typedef matrix< I1*R_incr, R, T_val > u1_type;
		typedef matrix< I2*R_incr, R, T_val > u2_type;
		typedef matrix< I3*R_incr, R, T_val > u3_type;
		
		typedef matrix< I1*R_incr, R, T_res > u1_res_type;
		typedef matrix< I2*R_incr, R, T_res > u2_res_type;
		typedef matrix< I3*R_incr, R, T_res > u3_res_type;
		
		//incremental cp als (zang&golub, 2001)
		static void incremental_als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ );
		
	};
	
	
	
#define VMML_TEMPLATE_STRING		template< size_t R, size_t R_incr, size_t I1, size_t I2, size_t I3, typename T_val, typename T_res >
#define VMML_TEMPLATE_CLASSNAME		t3_ihopm< R, R_incr, I1, I2, I3, T_val, T_res >
	
	
	VMML_TEMPLATE_STRING
	void 
	VMML_TEMPLATE_CLASSNAME::incremental_als( const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, lambda_type& lambdas_ )
	{
		t3_res_type* data = new t3_res_type;
		data->cast_from( data_ );
		t3_res_type* approx_data = new t3_res_type;
		approx_data->zero();
		t3_res_type* residual_data = new t3_res_type;
		residual_data->cast_from( data_ );
		
		lambda_res_type* lambdas_res = new lambda_res_type;
		u1_res_type* u1_res = new u1_res;
		u2_res_type* u2_res = new u2_res;
		u3_res_type* u3_res = new u3_res;
		
		for ( size_t i= 0; i < R_incr; ++i )
		{
/*			% Results of the given iteration with the current residual as input
			P = cp_als(res, R);
			
			Pfull.lambda((iter-1)*R+1:iter*R)=P.lambda;
			
			Pfull.u{1}(:,(iter-1)*R+1:iter*R)=P.u{1};
			Pfull.u{2}(:,(iter-1)*R+1:iter*R)=P.u{2};
			Pfull.u{3}(:,(iter-1)*R+1:iter*R)=P.u{3};
			% Storage of the data
			
			res = res - tensor(P);*/
			
			t3_hopm< R, I1, I2, I3, T_res >::als( residual_data, u1_res, u2_res, u3_res, lambdas_res );
			
			//set lambdas und us to appropriate position
			
				
			
			t3_hopm< R, I1, I2, I3, T_res >::reconstruct( approx_data, u1_res, u2_res, u3_res, lambdas_res );
			
			residual_data = residual_data - approx_data;
		
		}
		

		delete u1_res;
		delete u2_res;
		delete u3_res;
		delete lambdas_res;
		delete residual_data;
	}
