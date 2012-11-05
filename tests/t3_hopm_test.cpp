#include "t3_hopm_test.hpp"
#include "vmmlib/t3_hopm.hpp"
#include <sstream>

namespace vmml
{
	
	bool
	t3_hopm_test::run()
	{
		bool ok = false;
		
		double precision = 0.001;

		//CP decomposition
		typedef tensor3< 4, 4, 4, double > cp_t3_type;
		
		cp_t3_type t3_cp_input;
		double data_in_cp[] = { 
			0.3780, 0.3150, 0.3386, 0.2047, 0.2913, 0.3071, 0.2835, 0.1024, 0.2362, 0.2835, 0.2677, 0.1024, 0.3543, 1.1181, 1.5354, 0.3858, 
			0.2520, 0.2283, 0.3228, 0.2835, 0.2677, 0.2598, 0.2992, 0.2126, 0.2441, 0.2205, 0.2441, 0.2913, 0.9213, 0.6457, 0.4331, 0.1890,
			0.4409, 0.4409, 0.5591, 0.5039, 0.2362, 0.4409, 0.5984, 0.6142, 0.2520, 0.2835, 0.3465, 0.3543, 0.5748, 0.2835, 0.2992, 0.2835,
			0.3386, 0.3150, 0.4488, 0.4173, 0.2756, 0.3150, 0.3465, 0.3386, 0.2835, 0.2677, 0.2362, 0.2913, 0.2598, 0.2520, 0.2756, 0.3071 };
		t3_cp_input.set(data_in_cp, data_in_cp + 64);
		
		//std::cout << "cpinput\n" << t3_cp_input	<< std::endl;
		
		
		//rank-R approximation with R = I
		typedef matrix< 4, 4, double > cp1_u_type;
		typedef vector< 4, double> cp1_lambda_type;
		
		cp1_lambda_type lambda_check;
		lambda_check.at(0) = 3.0476; lambda_check.at(1) = 1.4161; 
		lambda_check.at(2) = 0.7849; lambda_check.at(3) = 0.6352;
                
		double data_u1_cp[] = { 
                0.4671,   -0.0204,   -0.0469,   -0.3844,
                0.4149,   -0.0260,    0.0611,   -0.5919,
                0.3321,    0.0548,    0.0140,   -0.2944,
                0.7067,    0.9979,    0.9969,    0.6444,
		};
	    cp1_u_type u1_check; u1_check.set( data_u1_cp, data_u1_cp + 16);
		double data_u2_cp[] = { 
                -0.5169,   0.1747,    -0.9317,   -0.3306,
                -0.5142,    -0.5134,    -0.2352,   -0.0583,
                -0.5675,    -0.8373,   0.2660,    0.2163,
                -0.3825,   0.0694,   0.0771,    0.9168,
                };
		cp1_u_type u2_check; u2_check.set( data_u2_cp, data_u2_cp + 16);
		double data_u3_cp[] = { 
		-0.4495,    -0.8320,    -0.1279,    0.2836,
                -0.3973,    -0.0952,    -0.6569,   -0.2261,
                -0.6338,   0.4455,   0.5090,   -0.8426,
                -0.4882,   0.3167,   0.5413,   -0.3981,
                };
		cp1_u_type u3_check; u3_check.set( data_u3_cp, data_u3_cp + 16);
		
		cp1_lambda_type lambda;
		cp1_u_type u1;
		cp1_u_type u2;
		cp1_u_type u3;
		typedef t3_hopm< 4, 4, 4, 4, double > hopm_type;
		hopm_type::als( t3_cp_input, u1, u2, u3, lambda, hopm_type::init_hosvd(), 20, -1 );
				
		ok = u1.equals( u1_check, precision );
		ok = u2.equals( u2_check, precision ) && ok;
		ok = u3.equals( u3_check, precision ) && ok;
		ok = lambda.equals( lambda_check, precision ) && ok;
		
		if( ok )
		{	
			log( "HOPM/CP-ALS: rank-R approximation (R = I)", ok  );
		} else
		{
			std::stringstream error;
			error 
			<< "HOPM/CP-ALS: rank-R approximation (R = I)" << std::setprecision(16) << std::endl
			<< " lambda should be: " << lambda_check << "lambda is: " << lambda	<< std::endl
			<< " u1 should be: " << std::endl << u1_check << std::endl
			<< " u1 is: " << std::endl << u1 << std::endl
			<< " u2 should be: " << std::endl << u2_check << std::endl
			<< " u2 is: " << std::endl << u2 << std::endl
			<< " u3 should be: " << std::endl << u3_check << std::endl
			<< " u3 is: " << std::endl << u3 << std::endl;
			
			log_error( error.str() );
		}
		
		//test with uint8 input values
		
		
		
		//rank-R approximation (R < I)
		typedef matrix< 4, 2, double > cp2_u_type;
		typedef vector< 2, double> cp2_lambda_type;
		
		cp2_lambda_type lambda_check2; 
    
		lambda_check2.at(0) = 2.8700; lambda_check2.at(1) = 1.5134;

		double data_u1_cp2[] = {
		0.5130,    0.0172,
                0.4671,   -0.0112,
                0.3632,    0.0641,
                0.6218,    0.9977,
                };
	    cp2_u_type u1_check2; u1_check2.set( data_u1_cp2, data_u1_cp2 + 8);
		double data_u2_cp2[] = { 
                -0.4996,    -0.0615,
                -0.4958,    -0.5782,
                -0.5648,    -0.8075,
                -0.4309,    -0.0994, 
                };
		cp2_u_type u2_check2; u2_check2.set( data_u2_cp2, data_u2_cp2 + 8);
		double data_u3_cp2[] = { 
		-0.3748,    -0.9273,
                -0.4621,    -0.0808,
                -0.6483,   0.3165,
                -0.4750,   0.1828,
                };
		cp2_u_type u3_check2; u3_check2.set( data_u3_cp2, data_u3_cp2 + 8);
		
		cp2_lambda_type lambda_2;
		cp2_u_type u1_2;
		cp2_u_type u2_2;
		cp2_u_type u3_2;
		
		typedef t3_hopm< 2, 4, 4, 4, double > hopm2_type;
		hopm2_type::als( t3_cp_input, u1_2, u2_2, u3_2, lambda_2, hopm2_type::init_hosvd(), 20, -1 );
		
		precision = 0.0001;
		ok = u1_2.equals( u1_check2, precision );
		ok = u2_2.equals( u2_check2, precision ) && ok;
		ok = u3_2.equals( u3_check2, precision ) && ok;
		ok = lambda_2.equals( lambda_check2, precision ) && ok;
		
		if( ok)
		{	
			log( "HOPM/CP-ALS: rank-R approximation (R < I)", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "HOPM/CP-ALS: rank-R approximation (R < I)" << std::setprecision(16) << std::endl
			<< " lambda should be:\n" << lambda_check2 << "lambda is:\n" << lambda_2	<< std::endl
			<< " u1 should be: " << std::endl << u1_check2 << std::endl
			<< " u1 is: " << std::endl << u1_2 << std::endl
			<< " u2 should be: " << std::endl << u2_check2 << std::endl
			<< " u2 is: " << std::endl << u2_2 << std::endl
			<< " u3 should be: " << std::endl << u3_check2 << std::endl
			<< " u3 is: " << std::endl << u3_2 << std::endl;
			
			log_error( error.str() );
		}
		
#if 0
		typedef tensor3< 4, 4, 4, unsigned char > cp_t3_uint8_type;
		
		cp_t3_uint8_type t3_cp_input2;
		t3_cp_input2.fill_increasing_values();		
		std::cout << "cpinput2\n" << t3_cp_input2	<< std::endl;
		
		cp2_lambda_type lambda_10;
		cp2_u_type u1_10;
		cp2_u_type u2_10;
		cp2_u_type u3_10;
		
		
		//should be converted to double
		hopm2_type::als( t3_cp_input2, u1_10, u2_10, u3_10, lambda_10, hopm2_type::init_hosvd(), 50 );
		
		
#endif
		
		
		
		
		
		//rank-R approximation with R > I
		typedef matrix< 4, 6, double > cp3_u_type;
		typedef vector< 6, double> cp3_lambda_type;
		
		cp3_lambda_type lambda_check3; 
		lambda_check3.at(0) = 2.949088573455811; lambda_check3.at(1) = 1.269379019737244; 
		lambda_check3.at(2) = 0.7361284494400024; lambda_check3.at(5) = 0.3052012622356415; 
		lambda_check3.at(3) = 0.6072210669517517; lambda_check3.at(4) = 0.3266942501068115; 
		
		double data_u1_cp3[] = { 
			0.4733901800329614, 0.1651671976476345, 0.03429142531960762, 0.5458835313216009, 0.003786561168133651, 0.1165964383253547, 
			0.4347951789645572, 0.166634467441311, -0.3642364139173337, 0.6078584294591635, -0.6301898404211475, 0.8419814543957256, 
			0.3151639631872558, 0.04782939064330856, 0.09616514805376049, 0.4610411401757213, -0.1384341324514788, -0.5232021638435955, 
			0.6982310263126678, -0.9709094192933075, -0.9256933602581305, -0.3463529515288832, 0.7639911111766255, 0.06109007098440862
		};
		cp3_u_type u1_check3; u1_check3.set( data_u1_cp3, data_u1_cp3 + 24);
		
		double data_u2_cp3[] = { 
			-0.4672415036543002, 0.3213422331793286, -0.9265545290349607, 0.3166075143499547, 0.8067300709812761, -0.776617286082291, 
			-0.5258024640955455, -0.4592700321284731, -0.3573195170303192, -0.05089615712113303, -0.1289753777336823, -0.3767214814465603, 
			-0.6067694317031729, -0.8228066042210956, 0.1175559547231347, -0.02855775384536048, -0.4263246085807784, -0.02708859484252683, 
			-0.3701999497089777, 0.09380564382905007, 0.00025491793523442, 0.946759588144082, -0.3883288202431613, -0.5041951253278271
		};
		cp3_u_type u2_check3; u2_check3.set( data_u2_cp3, data_u2_cp3 + 24);
		
		double data_u3_cp3[] = { 
			-0.5580260829866757, 0.858622294938589, 0.1629322461162687, -0.399424471555235, 0.157199202283234, 0.3001292203805404, 
			-0.3904421270404547, 0.04710793275822141, 0.8503159915200578, 0.3206114044275722, 0.3149977299563115, 0.857234067975567,
			-0.5808629004086727, -0.3794478834543623, -0.4222111559057403, 0.6514441744236624, 0.8237143771440413, -0.3939488133263827, 
			-0.4458252203134067, -0.3414204168582224, -0.2686141053063062, 0.5597221690117998, -0.4444766200938916, 0.1409841702268896
			 };
		cp3_u_type u3_check3; u3_check3.set( data_u3_cp3, data_u3_cp3 + 24);
		
		cp3_lambda_type lambda_3;
		cp3_u_type u1_3;
		cp3_u_type u2_3;
		cp3_u_type u3_3;
		
		
		//test norm ktensor
		double check_norm = 3.309332639286851;
		double nrm = hopm2_type::norm_ktensor( u1_2, u2_2, u3_2, lambda_2 );
		
		precision = 0.00000001;
		ok = (check_norm - nrm < precision );
		
		if( ok)
		{	
			log( "norm ktensor", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "norm ktensor" << std::setprecision(16) << std::endl
			<< "should be: " << check_norm << std::endl
			<< "is: " << nrm
			<< std::endl;
			
			log_error( error.str() );
		}
		
		
#if 0 //test hosvd with SVD
		t3_hopm< 6, 4, 4, 4, double >::als( t3_cp_input, u1_3, u2_3, u3_3, lambda_3, init_hosvd_e, 50 );
		
		ok = u1_3.equals( u1_check3, precision );
		ok = u2_3.equals( u2_check3, precision ) && ok;
		ok = u3_3.equals( u3_check3, precision ) && ok;
		ok = lambda_3.equals( lambda_check3, precision ) && ok;
		if( ok)
		{	
			log( "HOPM/CP-ALS: rank-R approximation (R > I) - init with DCT", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "HOPM/CP-ALS: rank-R approximation (R > I) - init with DCT" << std::setprecision(16) << std::endl
			<< " lambda should be:\n" << lambda_check3 << "\nlambda is:\n" << lambda_3	<< std::endl
			<< " u1 should be: " << std::endl << u1_check3 << std::endl
			<< " u1 is: " << std::endl << u1_3 << std::endl
			<< " u2 should be: " << std::endl << u2_check3 << std::endl
			<< " u2 is: " << std::endl << u2_3 << std::endl
			<< " u3 should be: " << std::endl << u3_check3 << std::endl
			<< " u3 is: " << std::endl << u3_3 << std::endl;
			
			log_error( error.str() );
		}
#endif
		
//FIXME: check test on linux
#if 0
		//rank-R approximation with R > I (init with DCT)
#define D 5
		typedef matrix< 4, D, double > cp4_u_type;
		typedef vector< D, double> cp4_lambda_type;
		
		cp4_lambda_type lambda_check4; 
		lambda_check4.at(0) = 2.498177053837107; lambda_check4.at(1) = 1.562078744247265; 
		lambda_check4.at(2) = 1.489490792770405; lambda_check4.at(3) = 1.344552972087995; 
		lambda_check4.at(4) = 1.00387849517148; 
		
		double data_u1_cp4[] = { 
			0.4934237561073131, 0.01845047663911101, 0.1782548118024551, -0.09798657543290282, 0.5487301760825246,
			0.3947501379000385, 0.03515044879384801, 0.03535751132544677, -0.1432457637568596, 0.7087410243000387,
			0.3590339322332908, -0.01781538889190339, 0.08821210931538905, -0.107110463271061, 0.3406927700691575,
			0.6868769620841511, 0.9990528703627306, -0.9793843434690348, -0.9789824466549462, -0.2837424725929791 };
		cp4_u_type u1_check4; u1_check4.set( data_u1_cp4, data_u1_cp4 + 24);
		
		double data_u2_cp4[] = { 
			0.6070441220297691, 0.8933125585845675, 0.8883053491235372, 0.1060204947494846, 0.0760074664421849,
			0.4713983671605251, 0.4218304924097846, 0.3679361948033744, -0.5453786675836705, -0.3260460774679831,
			0.476310434351972, 0.14017888770209, -0.07276909373517129, -0.8311900202676862, -0.6431636455439254, 
			0.4270941154768272, -0.06634446315267831, -0.2650306062874359, 0.0210929809448136, -0.688663448618721 };
		cp4_u_type u2_check4; u2_check4.set( data_u2_cp4, data_u2_cp4 + 24);
		
		double data_u3_cp4[] = { 
			0.4879274471565801, 0.03925648729505715, 0.05178189302040599, 0.9909298694606465, 0.1945927925001754,
			0.5390409774503466, -0.299246573525204, -0.6105174583715575, 0.1122841549455829, 0.00804843394886336,
			0.3757757219193364, 0.7753824955374058, 0.7175999088679488, -0.07368107715613589, -0.9593409048890027,
			0.5745904957068463, -0.5547002812806773, -0.3311154471954886, -0.004621820850501224, -0.204288756481369 };
		cp4_u_type u3_check4; u3_check4.set( data_u3_cp4, data_u3_cp4 + 24);
		
		cp4_lambda_type lambda_4;
		cp4_u_type u1_4;
		cp4_u_type u2_4;
		cp4_u_type u3_4;
		
		typedef t3_hopm< D, 4, 4, 4, double > hopm4_type;
		hopm4_type::als( t3_cp_input, u1_4, u2_4, u3_4, lambda_4, hopm4_type::init_dct(), 100 );
		
		ok = u1_4.equals( u1_check4, precision );
		ok = u2_4.equals( u2_check4, precision ) && ok;
		ok = u3_4.equals( u3_check4, precision ) && ok;
		ok = lambda_4.equals( lambda_check4, precision ) && ok;
		if( ok)
		{	
			log( "HOPM/CP-ALS with init DCT: rank-R approximation (R > I)", ok  );
		} 
		else
		{
			std::stringstream error;
			error 
			<< "HOPM/CP-ALS with init DCT: rank-R approximation (R > I)" << std::setprecision(16) << std::endl
			<< " lambda should be:\n" << lambda_check4 
			<< "\nlambda is:\n" << lambda_4	<< std::endl
			<< " u1 should be: " << std::endl << u1_check4 << std::endl
			<< " u1 is: " << std::endl << u1_4 << std::endl
			<< " u2 should be: " << std::endl << u2_check4 << std::endl
			<< " u2 is: " << std::endl << u2_4 << std::endl
			<< " u3 should be: " << std::endl << u3_check4 << std::endl
			<< " u3 is: " << std::endl << u3_4 << std::endl;
			
			log_error( error.str() );
		}
		
#endif
		
		
		
		return ok;
	}
	
} //end vmml namespace

