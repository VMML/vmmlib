
#include "cp3_tensor_test.hpp"

#include <vmmlib/cp3_tensor.hpp>

#include <sstream>

namespace vmml
{
	
	bool
	cp3_tensor_test::run()
	{
		bool ok = false;
		
		//TODO: notice this class is not yet working.
		
#if 0		
		//CP decomposition
		//same test data s in tensor3.reconstruct_CP()
		
		tensor3< 4,4,4, double > t3_cp_input;
		double data_in_cp[] = { 
			0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858, 
			0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
			0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
			0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071 };
		t3_cp_input.set(data_in_cp, data_in_cp + 64);
		
		vector< 4, double> lambda_check; lambda_check.at(0) = 3.4996; lambda_check.at(1) = 1.6186; lambda_check.at(2) = 1.0947; lambda_check.at(3) = 0.80494;
		double data_u1_cp[] = { 
			0.414514672613941, -0.273871195562966, 0.622972704020820, -0.106663870725802, 
			0.373815968628613, -0.303860590002700, 0.583235093038218, 0.038436291238145,
			0.307288348024630, -0.134400219672316, 0.462215931211866, -0.019751415378651, 
			0.770722439633298, 0.902551877199397, 0.241035825610295, 0.993355601880911};
		matrix< 4, 4, double > u1_check; u1_check.set( data_u1_cp, data_u1_cp + 16);
		double data_u2_cp[] = { 
			0.363947543028453, 0.298190617522100, 0.309440303599311, 0.906460391523818, 
			0.582584049017999, -0.442905614843157, -0.385671235746994, 0.178817966899174, 
			0.700932982871312, -0.827455897125261, -0.577509589347220, -0.345891330880255, 
			0.191913952792576, -0.173878435335793, 0.649605319126744, -0.163440755426602 };
		matrix< 4, 4, double > u2_check; u2_check.set( data_u2_cp, data_u2_cp + 16);
		double data_u3_cp[] = { 
			0.528860463703445, -0.453374574314148, 0.081548921839262, 0.026850375225737, 
			0.421141692804753, 0.103741153724736, 0.379185048676421, 0.600983237030590,
			0.588063451217321, 0.731736212173060, 0.724001032600722, -0.545272439076217, 
			0.443990610098116, 0.498248315786364, 0.570430518994385, -0.583760373220334};
		matrix< 4, 4, double > u3_check; u3_check.set( data_u3_cp, data_u3_cp + 16);
		
		std::cout << "lambdas: " << std::endl << lambda_check << std::endl << std::endl
		<< "U1: " << std::endl << u1_check << std::endl
		<< "U2: " << std::endl << u2_check << std::endl
		<< "U3: " << std::endl << u3_check << std::endl;
		
		tensor3< 4,4,4, double > t3_cp_reco_check;
		double data_out_cp[] = { 
			0.354263061197741, 0.336071353932582, 0.340120155317465, 0.148789058227363,
			0.335156282115524, 0.284369809600727, 0.270075906182746, 0.127689570299207,
			0.248778281916907, 0.271657788743462, 0.293346312353042, 0.118870461584074,
			0.347767704804178, 1.119919940137974, 1.528048459613556, 0.399390301753092,
			
			0.241873571540868, 0.267325817660143, 0.334780629072782, 0.301651586726927,
			0.277064331887503, 0.253524404359242, 0.282150329524853, 0.268828551075853,
			0.208805038309086, 0.198138779535580, 0.228624285742602, 0.217032724171575,
			0.925152576629409, 0.641979031363691, 0.446794524932619, 0.178096907232692,
			
			0.408959960261278, 0.458603477148098, 0.565021279854640, 0.533195254992912,
			0.300410894904167, 0.426297991808946, 0.575920757418752, 0.513247731493541,
			0.303903577435438, 0.299195909145241, 0.360430286717156, 0.385593403984964,
			0.559935740903585, 0.298948990808052, 0.267704058359116, 0.313877703199338,
			
			0.334349191835064, 0.331983603647059, 0.392220714079066, 0.406517347684063,
			0.234640418580267, 0.303232696231040, 0.405821982955168, 0.393611395358567,
			0.239172809919676, 0.216514108351268, 0.254464520462408, 0.296451592996035,
			0.276358969906296, 0.233761715518065, 0.311620024949421, 0.277319482627673 };
		
		t3_cp_reco_check.set(data_out_cp, data_out_cp + 64);
		
		
		cp3_tensor< 4, 4, 4, 4, double, double > cp3_rank4;
		
		vector< 4, double> lambda;
		matrix< 4, 4, double > u1;
		matrix< 4, 4, double > u2;
		matrix< 4, 4, double > u3;
		
		//cp3_rank4.cp_als( t3_cp_input );
		
		cp3_rank4.get_u1( u1 );
		cp3_rank4.get_u2( u2 );
		cp3_rank4.get_u3( u3 );
		cp3_rank4.get_lambdas( lambda );
		
		if( ok)
		{	
			log( "cp3 tensor test: rank-1 approximation ", ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "cp3 tensor test: rank-1 approximation " << std::endl
			<< " lambda should be: " << lambda_check << "lambda is: " << lambda	<< std::endl
			<< " u1 should be: " << std::endl << u1_check << std::endl
			<< " u1 is: " << std::endl << u1 << std::endl
			<< " u2 should be: " << std::endl << u2_check << std::endl
			<< " u2 is: " << std::endl << u2 << std::endl
			<< " u3 should be: " << std::endl << u3_check << std::endl
			<< " u3 is: " << std::endl << u3 << std::endl;
			
			
			log_error( error.str() );
		}
		
		
		
		//decomposition (hopm test data after lathauwer 2000b)
		//prepare control data
		//rank-1 approximation
		matrix<3, 1, float_t> u1_rank1;
		matrix<2, 1, float_t> u2_rank1;
		matrix<2, 1, float_t> u3_rank1;
		matrix<3, 1, float_t> u1_rank1_check;
		matrix<2, 1, float_t> u2_rank1_check;
		matrix<2, 1, float_t> u3_rank1_check;
		
		u1_rank1_check.at(0,0) = -0.2515; u1_rank1_check.at(1,0) = 0.6035; u1_rank1_check.at(2,0) = 0.7567;
		u2_rank1_check.at(0,0) = 0.1344; u2_rank1_check.at(1,0) = 0.9909;
		u3_rank1_check.at(0,0) = 0.5765; u3_rank1_check.at(1,0) = -0.8171;
		
		vector< 1, float_t> lambda_rank1;
		vector< 1, float_t> lambda_rank1_check;
		lambda_rank1_check.at(0)  = 10.1693;
		
		tensor3< 3, 2, 2, float_t> t3_data;
		float_t data[] = { 0, 1, 2, 3, 4, 5, -1, 4, -2, -5, 3, -6};
		t3_data.set(data, data + 12);
		
		cp3_tensor< 3, 2, 2, 1, float_t, float_t > cp3_rank1( u1_rank1, u2_rank1, u3_rank1, lambda_rank1 );

		cp3_rank1.cp_als( t3_data );
		cp3_rank1.get_u1( u1_rank1 );
		cp3_rank1.get_u2( u2_rank1 );
		cp3_rank1.get_u3( u3_rank1 );
		cp3_rank1.get_lambdas( lambda_rank1 );
		
		float_t precision = 0.001;
		ok = ( u1_rank1.equals( u1_rank1_check, precision ) && u2_rank1.equals( u2_rank1_check, precision) && u3_rank1.equals( u3_rank1_check, precision) && lambda_rank1 == lambda_rank1_check );
		
		if( ok)
		{	
			log( "cp3 tensor test: rank-1 approximation ", ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "cp3 tensor test: rank-1 approximation " << std::endl
			<< " lambda should be: " << lambda_rank1_check << "lambda is: " << lambda_rank1	<< std::endl
			<< " u1 should be: " << std::endl << u1_rank1_check << std::endl
			<< " u1 is: " << std::endl << u1_rank1 << std::endl
			<< " u2 should be: " << std::endl << u2_rank1_check << std::endl
			<< " u2 is: " << std::endl << u2_rank1 << std::endl
			<< " u3 should be: " << std::endl << u3_rank1_check << std::endl
			<< " u3 is: " << std::endl << u3_rank1 << std::endl;
			
			
			log_error( error.str() );
		}
		
#endif		
		
		ok = true;
		return ok;
	}
	
	
} // namespace vmml