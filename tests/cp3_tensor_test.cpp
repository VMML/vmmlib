#include "cp3_tensor_test.hpp"
#include <vmmlib/cp3_tensor.hpp>
#include <sstream>

namespace vmml
{
	
	bool
	cp3_tensor_test::run()
	{
		bool ok = false;
		
		//CP decomposition
		typedef tensor3< 4, 4, 4, double > cp_t3_type;
		
		cp_t3_type t3_cp_input;
		double data_in_cp[] = { 
			0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858, 
			0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
			0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
			0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071 };
		t3_cp_input.set(data_in_cp, data_in_cp + 64);
		
		//FIXME: test for reconstruction
		//same test data s in tensor3.reconstruct_CP()
		cp_t3_type t3_cp_reco_check;
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
		//log( "cp3 tensor reconstruction ", false  );

		
		//rank-R approximation with R = I
		typedef matrix< 4, 4, double > cp1_u_type;
		typedef vector< 4, double> cp1_lambda_type;
		typedef cp3_tensor< 4, 4, 4, 4, double, double > cp1_decomp_type;
		
		cp1_lambda_type lambda_check; 
		lambda_check.at(0) = 2.669350624084473; lambda_check.at(1) = 1.478976130485535; 
		lambda_check.at(2) = 0.7834770679473877; lambda_check.at(3) = 0.6970056295394897;       
		double data_u1_cp[] = { 
			0.509520947933197, -0.05159009620547295, 0.004901818465441465, 0.3843765556812286,
			0.4130257368087769, -0.1193722859025002, -0.0818917378783226, 0.7114874124526978,
			0.3542648255825043, -0.101799413561821, -0.04321824759244919, 0.2759316265583038,
			0.6665542721748352, -0.9862684607505798, -0.9956916570663452, -0.5195209980010986 };
	    cp1_u_type u1_check; u1_check.set( data_u1_cp, data_u1_cp + 16);
		double data_u2_cp[] = { 
			-0.5690991282463074, 0.1410824805498123, -0.9016433954238892, 0.2951040863990784,
			-0.4794410169124603, -0.5260365605354309, -0.1053012236952782, -0.2450571358203888,
			-0.508577287197113, -0.8381729125976562, 0.4171290695667267, -0.5841856598854065,
			-0.4331415891647339, 0.02911213971674442, 0.04420609027147293, -0.7152535915374756 };
		cp1_u_type u2_check; u2_check.set( data_u2_cp, data_u2_cp + 16);
		double data_u3_cp[] = { 
			-0.4507713615894318, 0.9092601537704468, 0.1635313779115677, 0.3119933307170868,
			-0.3851732909679413, 0.3323416411876678, 0.8324259519577026, -0.2073864340782166,
			-0.6179288029670715, -0.1530936360359192, -0.2699224948883057, -0.8566842675209045,
			-0.5163435935974121, -0.1983867585659027, -0.4554848372936249, -0.3546026051044464	 };
		cp1_u_type u3_check; u3_check.set( data_u3_cp, data_u3_cp + 16);
		
		cp1_decomp_type cp3_rank4;
		cp1_lambda_type lambda;
		cp1_u_type u1;
		cp1_u_type u2;
		cp1_u_type u3;
		
		cp3_rank4.cp_als( t3_cp_input, 100 );
		cp3_rank4.get_u1( u1 );
		cp3_rank4.get_u2( u2 );
		cp3_rank4.get_u3( u3 );
		cp3_rank4.get_lambdas( lambda );
		
		double precision = 0.0001;
		ok = u1.equals( u1_check, precision );
		ok = u2.equals( u2_check, precision );
		ok = u3.equals( u3_check, precision );
		ok = lambda.equals( lambda_check, precision );

		if( ok)
		{	
			log( "HOPM: rank-R approximation (R = I) (ALS)", ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOPM: rank-R approximation (R = I) (ALS)" << std::setprecision(16) << std::endl
			<< " lambda should be: " << lambda_check << "lambda is: " << lambda	<< std::endl
			<< " u1 should be: " << std::endl << u1_check << std::endl
			<< " u1 is: " << std::endl << u1 << std::endl
			<< " u2 should be: " << std::endl << u2_check << std::endl
			<< " u2 is: " << std::endl << u2 << std::endl
			<< " u3 should be: " << std::endl << u3_check << std::endl
			<< " u3 is: " << std::endl << u3 << std::endl;
			
			log_error( error.str() );
		}
		
		//rank-R approximation (R > I)
		typedef matrix< 4, 2, double > cp2_u_type;
		typedef vector< 2, double> cp2_lambda_type;
		typedef cp3_tensor< 4, 4, 4, 2, double, double > cp2_decomp_type;

		cp2_lambda_type lambda_check2; 
		lambda_check2.at(0) = 2.905300140380859; lambda_check2.at(1) = 1.484431028366089; 
		double data_u1_cp2[] = { 
			0.5060855746269226, -0.001956983935087919,
			0.4604957103729248, 0.02855747938156128,
			0.359317272901535, -0.05445058643817902,
			0.634596049785614, -0.998106062412262 };
	    cp2_u_type u1_check2; u1_check2.set( data_u1_cp2, data_u1_cp2 + 8);
		double data_u2_cp2[] = { 
			-0.4963421821594238, -0.03891453891992569,
			-0.4994568228721619, -0.5732658505439758,
			-0.569895327091217, -0.8119437694549561,
			-0.4235640466213226, -0.1029524803161621 };
		cp2_u_type u2_check2; u2_check2.set( data_u2_cp2, data_u2_cp2 + 8);
		double data_u3_cp2[] = { 
			-0.3882103264331818, 0.9169266819953918,
			-0.466062068939209, 0.05367951840162277,
			-0.6417607665061951, -0.3428435325622559,
			-0.4692781269550323, -0.1970336884260178	 };
		cp2_u_type u3_check2; u3_check2.set( data_u3_cp2, data_u3_cp2 + 8);
		
		cp2_decomp_type cp3_rank2;
		cp2_lambda_type lambda_2;
		cp2_u_type u1_2;
		cp2_u_type u2_2;
		cp2_u_type u3_2;
		
		cp3_rank2.cp_als( t3_cp_input, 50 );
		cp3_rank2.get_u1( u1_2 );
		cp3_rank2.get_u2( u2_2 );
		cp3_rank2.get_u3( u3_2 );
		cp3_rank2.get_lambdas( lambda_2 );
		
		precision = 0.0001;
		ok = u1_2.equals( u1_check2, precision );
		ok = u2_2.equals( u2_check2, precision );
		ok = u3_2.equals( u3_check2, precision );
		ok = lambda_2.equals( lambda_check2, precision );
		
		if( ok)
		{	
			log( "HOPM: rank-R approximation (R < I) (ALS)", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "HOPM: rank-R approximation (R < I) (ALS)" << std::setprecision(16) << std::endl
			<< " lambda should be:\n" << lambda_check2 << "lambda is:\n" << lambda_2	<< std::endl
			<< " u1 should be: " << std::endl << u1_check2 << std::endl
			<< " u1 is: " << std::endl << u1_2 << std::endl
			<< " u2 should be: " << std::endl << u2_check2 << std::endl
			<< " u2 is: " << std::endl << u2_2 << std::endl
			<< " u3 should be: " << std::endl << u3_check2 << std::endl
			<< " u3 is: " << std::endl << u3_2 << std::endl;
			
			log_error( error.str() );
		}
		
		//rank-R approximation with R > I
		typedef matrix< 4, 6, double > cp3_u_type;
		typedef vector< 6, double> cp3_lambda_type;
		typedef cp3_tensor< 4, 4, 4, 6, double, double > cp3_decomp_type;

		cp3_lambda_type lambda_check3; 
		lambda_check3.at(0) = 2.949088573455811; lambda_check3.at(1) = 1.269379019737244; 
		lambda_check3.at(2) = 0.7361284494400024; lambda_check3.at(3) = 0.3052012622356415; 
		lambda_check3.at(4) = 0.6072210669517517; lambda_check3.at(5) = 0.3266942501068115; 
		double data_u1_cp3[] = { 
			0.4733894467353821, 0.1651687473058701, 0.0342925488948822, 0.1165953725576401, 0.5458840131759644, 0.003787092166021466,
			0.4347946643829346, 0.16663658618927, -0.3642348647117615, 0.841981053352356, 0.6078583598136902, -0.630190908908844,
			0.3151637613773346, 0.04783087223768234, 0.09616393595933914, -0.5232028365135193, 0.4610396325588226, -0.1384356766939163,
			0.6982319355010986, -0.9709087610244751, -0.9256940484046936, 0.0610918216407299, -0.3463544547557831, 0.7639899849891663 
		};
		
	    cp3_u_type u1_check3; u1_check3.set( data_u1_cp3, data_u1_cp3 + 24);
		double data_u2_cp3[] = { 
			-0.4672411382198334, 0.3213400542736053, -0.9265543818473816, -0.7766140103340149, 0.3166076242923737, 0.8067278861999512,
			-0.5258027911186218, -0.4592712223529816, -0.3573199212551117, -0.3767221570014954, -0.05090122297406197, -0.128976508975029,
			-0.6067699790000916, -0.8228066563606262, 0.1175554245710373, -0.02709124609827995, -0.02856527082622051, -0.4263260960578918,
			-0.3701991438865662, 0.09380617737770081, 0.0002545058378018439, -0.5041994452476501, 0.9467591047286987, -0.3883313834667206
			 };
		cp3_u_type u2_check3; u2_check3.set( data_u2_cp3, data_u2_cp3 + 24);
		double data_u3_cp3[] = { 
			-0.5580258369445801, 0.8586198091506958, 0.1629276275634766, 0.3001251816749573, -0.3994224965572357, 0.1571977436542511,
			-0.3904423415660858, 0.04710535332560539, 0.8503167033195496, 0.857235848903656, 0.3206111192703247, 0.3149990439414978,
			-0.5808628797531128, -0.3794513642787933, -0.4222103357315063, -0.3939476907253265, 0.6514447927474976, 0.8237122893333435,
			-0.4458253383636475, -0.3414232730865479, -0.2686156332492828, 0.1409851014614105, 0.5597229599952698, -0.4444800019264221
			 };
		cp3_u_type u3_check3; u3_check3.set( data_u3_cp3, data_u3_cp3 + 24);
		
		cp3_decomp_type cp3_rank6;
		cp3_lambda_type lambda_3;
		cp3_u_type u1_3;
		cp3_u_type u2_3;
		cp3_u_type u3_3;
		
		cp3_rank6.cp_als( t3_cp_input, 50 );
		cp3_rank6.get_u1( u1_3 );
		cp3_rank6.get_u2( u2_3 );
		cp3_rank6.get_u3( u3_3 );
		cp3_rank6.get_lambdas( lambda_3 );
		
		precision = 0.0001;
		ok = u1_3.equals( u1_check3, precision );
		ok = u2_3.equals( u2_check3, precision );
		ok = u3_3.equals( u3_check3, precision );
		ok = lambda_3.equals( lambda_check3, precision );
		if( ok)
		{	
			log( "HOPM: rank-R approximation (R > I) (ALS)", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "HOPM: rank-R approximation (R > I) (ALS)" << std::setprecision(16) << std::endl
			<< " lambda should be:\n" << lambda_check3 << "\nlambda is:\n" << lambda_3	<< std::endl
			<< " u1 should be: " << std::endl << u1_check3 << std::endl
			<< " u1 is: " << std::endl << u1_3 << std::endl
			<< " u2 should be: " << std::endl << u2_check3 << std::endl
			<< " u2 is: " << std::endl << u2_3 << std::endl
			<< " u3 should be: " << std::endl << u3_check3 << std::endl
			<< " u3 is: " << std::endl << u3_3 << std::endl;
			
			log_error( error.str() );
		}
		
		return ok;
	}
	
	
} // namespace vmml