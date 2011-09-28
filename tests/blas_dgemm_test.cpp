#include "blas_dgemm_test.hpp"

#include <vmmlib/blas_dgemm.hpp>

namespace vmml
{
	
	bool
	blas_dgemm_test::run()
	{
		bool ok = false;
		
		matrix< 3, 6, double > A;
		matrix< 3, 3, double > C;
		matrix< 3, 3, double > C_check;
		
		double AData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
		A = AData;
		
		blas_dgemm< 3, 6, 3, double > blas_cov;
		blas_cov.compute( A, C );
		
		double CData[] = { 91, 217, 343, 217, 559, 901, 343, 901, 1459 };
		C_check = CData;
		
		ok = C == C_check;
		
		log( "symmetric matrix-matrix multiplication (input left matrix) (MxK) x (KxN) = (MxN), while M=N", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix (left matrix)\n" << A << "\n"
            << "covariance matrix should be\n" << C_check << "\n"
            << "covariance matrix is\n" << C << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		
		//A*B = D (MxK, KxN, MxN) (input right matrix)
		matrix< 6, 3, double > E;
		matrix< 3, 3, double > F;
		matrix< 3, 3, double > F_check;
		
		double EData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
		E = EData;
		
		blas_dgemm< 3, 6, 3, double > blas_cov3;
		blas_cov3.compute_t( E, F );
		
		double FData[] = { 591, 642, 693, 642, 699, 756, 693, 756, 819 };
		F_check = FData;
		
		ok = F == F_check;
		
		log( "symmetric matrix-matrix multiplication (input right matrix) (MxK) x (KxN) = (MxN), while M=N", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix (right matrix) \n" << E << "\n"
            << "covariance matrix should be\n" << F_check << "\n"
            << "covariance matrix is\n" << F << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}

		//A*B = D (MxK, KxN, MxN)
		matrix< 6, 2, double > B;
		matrix< 3, 2, double > D;
		matrix< 3, 2, double > D_check;
		
		double BData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		B = BData;

		blas_dgemm< 3, 6, 2, double > blas_v_dgemm;
		blas_v_dgemm.compute( A, B, D );

		double DData[] = { 161, 182, 377, 434, 593, 686 };
		D_check = DData;
		
		ok = D == D_check;
		
		log( "matrix-matrix multiplication (MxK) x (KxN) = (MxN)", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix A\n" << A << "\n"
            << "input matrix B\n" << B << "\n"
            << "matrix C should be\n" << D_check << "\n"
            << "matrix C is\n" << D << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		
		//A*B = D (MxK, KxN, MxN) (input B is given as transpose)
		matrix< 2, 6, double > B_t;
		matrix< 3, 2, double > G;
		matrix< 3, 2, double > G_check;
		
		double BtData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		B_t = BtData;
		
		blas_dgemm< 3, 6, 2, double > blas_dgemm5;
		blas_dgemm5.compute_bt( A, B_t, G );
		
		double GData[] = { 91, 217, 217, 559, 343, 901 };
		G_check = GData;
		
		ok = G == G_check;
		
		log( "matrix-matrix multiplication (input B is transpose) (MxK) x (NxK) = (MxN)", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix A\n" << A << "\n"
            << "input matrix B_t\n" << B_t << "\n"
            << "matrix C should be\n" << G_check << "\n"
            << "matrix C is\n" << G << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		
		//A*B = H (MxK, KxN, MxN) (input A+B are given as transpose)
		matrix< 2, 5, double > At;
		matrix< 6, 2, double > Bt;
		matrix< 5, 6, double > H;
		matrix< 5, 6, double > H_check;
		
		Bt = BtData;
		At = BtData;
		
		blas_dgemm< 5, 2, 6, double > blas_dgemm6;
		blas_dgemm6.compute_t( At, Bt, H );

		double HData[] = { 
			13, 27, 41, 55, 69, 83,
			16, 34, 52, 70, 88, 106,
			19, 41, 63, 85, 107, 129,
			22, 48, 74, 100, 126, 152,
			25, 55, 85, 115, 145, 175 };
		H_check = HData;
		
		ok = H == H_check;
				
		log( "matrix-matrix multiplication (input A+B are transpose) (MxK) x (NxK) = (MxN)", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input matrix A\n" << A << "\n"
            << "input matrix B_t\n" << B << "\n"
            << "matrix H should be\n" << H_check << "\n"
            << "matrix H is\n" << H << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
	
		//vector vector outer product
		vector< 5, double > vA;
		vector< 6, double > vB;
		vA = BtData;
		vB = BtData;
		blas_dgemm< 5, 1, 6, double > blas_dgemm7;
		blas_dgemm7.compute_vv_outer( vA, vB, H );
		
		double H2Data[] = {
			1, 2, 3, 4, 5, 6,
			2, 4, 6, 8, 10, 12,
			3, 6, 9, 12, 15, 18,
			4, 8, 12, 16, 20, 24,
			5, 10, 15, 20, 25, 30};
		H_check = H2Data;
		
	
		ok = (H == H_check) ;
		
		log( "vector vector outer product with dgemm (M)' x (N) = (MxN)", ok );
		if ( ! ok )
		{
			std::stringstream ss;
			ss
            << "input vector A\n" << vA << "\n"
            << "input vector B_t\n" << vB << "\n"
            << "matrix H should be\n" << H_check << "\n"
            << "matrix H is\n" << H << "\n"
            << std::endl;
			log_error( ss.str() );
            
		}
		

		return ok;
	}
	
	
	
} // namespace vmml

