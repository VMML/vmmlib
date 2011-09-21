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
		

		//rank-R approximation with R = I
		typedef matrix< 4, 4, double > cp1_u_type;
		typedef vector< 4, double> cp1_lambda_type;
		
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
		
		cp1_lambda_type lambda;
		cp1_u_type u1;
		cp1_u_type u2;
		cp1_u_type u3;
		t3_hopm< 4, 4, 4, 4, double >::als( t3_cp_input, u1, u2, u3, lambda, init_hosvd_e, 100 );
				
		ok = u1.equals( u1_check, precision );
		ok = u2.equals( u2_check, precision ) && ok;
		ok = u3.equals( u3_check, precision ) && ok;
		ok = lambda.equals( lambda_check, precision ) && ok;
		
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
		
		cp2_lambda_type lambda_2;
		cp2_u_type u1_2;
		cp2_u_type u2_2;
		cp2_u_type u3_2;
		
		t3_hopm< 2, 4, 4, 4, double >::als( t3_cp_input, u1_2, u2_2, u3_2, lambda_2, init_hosvd_e, 50 );
		
		precision = 0.0001;
		ok = u1_2.equals( u1_check2, precision );
		ok = u2_2.equals( u2_check2, precision ) && ok;
		ok = u3_2.equals( u3_check2, precision ) && ok;
		ok = lambda_2.equals( lambda_check2, precision ) && ok;
		
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
		
		t3_hopm< 6, 4, 4, 4, double >::als( t3_cp_input, u1_3, u2_3, u3_3, lambda_3, init_hosvd_e, 50 );
		
		ok = u1_3.equals( u1_check3, precision );
		ok = u2_3.equals( u2_check3, precision ) && ok;
		ok = u3_3.equals( u3_check3, precision ) && ok;
		ok = lambda_3.equals( lambda_check3, precision ) && ok;
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
	
} //end vmml namespace