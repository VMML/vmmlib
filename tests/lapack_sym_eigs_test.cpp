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
		matrix< 4, 3, double > eigxvectors;
		vector< 3, double > eigxvalues;
		lapack_sym_eigs< 4, double >  eigs;
		eigs.compute_x( A, eigxvectors, eigxvalues);
		
		matrix< 4, 3, double > eigxvectors_check;
		vector< 3, double > eigxvalues_check;
		vector< 4, double > eigvector_check;
		double first_eigvalue = -2.776873597454109;
		double data_eigxvalues[] = {-2.776873597454109, -1.505827874899289, 0.8483149985370052 };
		double data_eigxvectors[] = {
			-0.846872951394099, -0.05486702445763353, 0.08543580914755927, 
			-0.2525272340154449, -0.8141033734762484, -0.2624902204189269,
			0.4028336856489185, -0.3565378545075824, -0.4780761764820428,
			0.2382461373958975, -0.4550890739458455, 0.8338121947897114 };
		eigxvalues_check = data_eigxvalues;
		eigxvectors_check.set( data_eigxvectors, data_eigxvectors + 12);
		eigxvectors_check.get_column( 0, eigvector_check );
	
		//get first eigvalue and eigvector
		vector< 4, double > eigvector;
		double eigvalue;
		eigs.compute_1st( A, eigvector, eigvalue);
				
		ok = eigxvalues.equals( eigxvalues_check, precision );
		ok = ok && eigxvectors.equals( eigxvectors_check, precision );
		ok = ok && (eigvector.equals( eigvector_check, precision ));
		//ok = ok && (fabs(eigvalue - first_eigvalue) < precision);
		
		if ( ok ) {
			log( "symmetric eigenvalue decomposition (x largest eigenvalues) using lapack", ok );
		} else {
			std::stringstream error;
			error 
			<< std::setprecision(10) 
			<< "symmetric eigenvalue decomposition (x largest eigenvalues) using lapack: " << std::endl
			<< "eigenvalues should be: " << std::endl << eigxvalues_check << std::endl
			<< "are: " << std::endl << eigxvalues << std::endl	
			<< "eigenvectors should be: " << std::endl << eigxvectors_check
			<< "are: " << std::endl << eigxvectors  << std::endl
			<< "first eigenvalue should be: " << first_eigvalue << ", is: " << eigvalue  << std::endl
			<< "first eigenvector should be:\n" << eigvector_check << "\n is:\n" << eigvector
			<< std::endl;
			
			log_error( error.str() );
		}
		//end compute x largest eigenvalues
		
		return true;
	}
	
	
	
} // namespace vmml

