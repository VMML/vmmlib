#include "lapack_sym_eigs_test.hpp"
#include <vmmlib/matrix.hpp>
#include <vmmlib/lapack_sym_eigs.hpp>


namespace vmml
{
	
	bool
	lapack_sym_eigs_test::run()
	{
		bool ok = true;
		
		matrix< 4, 4, double > A;
		
		double data[] = { -2.0000, -0.6714, 0.8698, 0.5792,
			-0.6714, -1.1242, -0.0365, -0.5731,
			0.8698, -0.0365, -0.4660, -0.8542,
			0.5792, -0.5731, -0.8542, 0.1188};
		A.set( data, data+16 );
		
		
		//compute all eigenvalues
		matrix< 4, 4, double > eigvectors;
		vector< 4, double > eigvalues;
		lapack_sym_eigs< 4, double >  eig;
		eig.compute_all( A, eigvectors, eigvalues);
		
		double data_eigvalues[] = {-2.776873597454109, -1.505827874899289, -0.03704459256313505, 0.8483149985370052};
		double data_eigvectors[] = {
			-0.846872951394099, -0.05486702445763353, 0.522012007849923, 0.08543580914755927,
			-0.2525272340154449, -0.8141033734762484, -0.45228815766001, -0.2624902204189269,
			0.4028336856489185, -0.3565378545075824, 0.6942974502977717, -0.4780761764820428,
			0.2382461373958975, -0.4550890739458455, 0.2022126025029584, 0.8338121947897114};
		matrix< 4, 4, double > eigvectors_check;
		vector< 4, double > eigvalues_check;
		eigvalues_check = data_eigvalues;
		eigvectors_check.set( data_eigvectors, data_eigvectors + 16);
				
		double precision = 1.0e-10;
		ok = eigvalues.equals( eigvalues_check, precision );
		ok = ok && eigvectors.equals( eigvectors_check, precision );
		
		if ( ok ) {
			log( "symmetric eigenvalue decomposition using lapack", ok );
		} else {
			std::stringstream error;
			error 
			<< "symmetric eigenvalue decomposition using lapack: " << std::endl
			<< "eigenvalues should be: " << std::endl << eigvalues_check << std::endl
			<< "are: " << std::endl << eigvalues << std::endl	
			<< "eigenvectors should be: " << std::endl << eigvectors_check << std::endl
			<< "are: " << std::endl << eigvectors << std::endl;
			
			log_error( error.str() );
		}
		//end compute all eigenvalues
		
		//compute x largest magnitude eigenvalues
		matrix< 3, 4, double > eigxvectors;
		vector< 3, double > eigxvalues;
		lapack_sym_eigs< 4, double >  eigs;
		eigs.compute_x( A, eigxvectors, eigxvalues);
		
#if 0
		std::cout << std::setprecision(16) << "eigvalues check:\n" << eigxvalues_check
		<< std::endl << "eigvalues is\n" << eigxvalues << std::endl << std::endl;
		
		std::cout << std::setprecision(16) << "eigvectors check:\n" << eigxvectors_check
		<< std::endl << "eigvectors is\n" << eigxvectors << std::endl;
#endif
		
#if 0		
		matrix< 3, 4, double > eigxvectors_check;
		vector< 3, double > eigxvalues_check;

		ok = eigxvalues.equals( eigvalues_check, precision );
		ok = ok && eigxvectors.equals( eigvectors_check, precision );
		
		if ( ok ) {
			log( "symmetric eigenvalue decomposition (x largest eigenvalues) using lapack", ok );
		} else {
			std::stringstream error;
			error 
			<< "symmetric eigenvalue decomposition (x largest eigenvalues) using lapack: " << std::endl
			<< "eigenvalues should be: " << std::endl << eigvalues_check << std::endl
			<< "are: " << std::endl << eigxvalues << std::endl	
			<< "eigenvectors should be: " << std::endl << eigvectors_check << std::endl
			<< "are: " << std::endl << eigxvectors << std::endl;
			
			log_error( error.str() );
		}
		//end compute x largest eigenvalues
#endif
		
		return true;
	}
	
	
	
} // namespace vmml

