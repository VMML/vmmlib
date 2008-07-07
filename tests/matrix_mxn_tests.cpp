#include "matrix_mxn_tests.hpp"

#include <sstream>

namespace vmml
{

bool
matrix_mxn_tests::run()
{
    bool ok = true;
    
    matrix_mxn< 2, 3 > m0;
    double data[] = { 1, 2, 3, 4, 5, 6 };
    
    m0.copyFrom1DimCArray( data );
    
    // tests copyFrom1DimCArray function with row_by_row data
	ok = true;
	{
		size_t tmp = 1;
		for( size_t rowIndex = 0; ok && rowIndex < 2; ++rowIndex )
		{
			for( size_t colIndex = 0; ok && colIndex < 3; ++colIndex, ++tmp )
			{
				ok = m0[ rowIndex ][ colIndex ] == tmp;
			}
		}
		log( "copyFrom1DimCArray( .. , row_by_row )", ok  );
		if ( ! ok )
		{
			std::stringstream error;
			error << m0 << std::endl;
			log_error( error.str() );
		}
	
	}
	
    matrix_mxn< 3, 2 > m1;
    m1.copyFrom1DimCArray( data, false );

    // tests copyFrom1DimCArray function with col_by_col data
	ok = true;
	{
		for( size_t tmp = 1, colIndex = 0; ok && colIndex < 2; ++colIndex )
		{
			for( size_t rowIndex = 0; ok && rowIndex < 3; ++rowIndex, ++tmp )
			{
				ok = m1[ rowIndex ][ colIndex ] == tmp;
			}
		}
		log( "copyFrom1DimCArray( .. , col_by_col )", ok );
		if ( ! ok )
		{
			std::stringstream error;
			error << m1 << std::endl;
			log_error( error.str() );
		}
	}

    // test operator== / operator !=
	ok = true;
	{
		matrix_mxn< 2, 3 > m0_copy;
		m0.copyFrom1DimCArray( data );
		m0_copy.copyFrom1DimCArray( data );
		ok = m0 == m0_copy;
		if ( ok )
		{
			ok = ! ( m0 != m0_copy );
		}
		log( "operator== / operator!=", ok );
		if ( ! ok )
		{
			std::stringstream error;
			error << m0 << m1 << std::endl;
			log_error( error.str() );
		}
	}
	
    // test transpose functionality 
	ok = true;
    {
        matrix_mxn< 3, 2 > m0t = m0.getTransposed();
        m1.copyFrom1DimCArray( data, false );
        
        ok = m1 == m0t;
		log( "getTransposed() / transposeTo()", ok );
        if ( !ok )
        {
			std::stringstream error;
			error << m1 << m0t << std::endl;
			log_error( error.str() );
		}
    }
    
    
    // test multiplication
    if ( ok )
    {
        matrix_mxn< 2, 3 > mul0;
        double mul0data[] = { 1, 0, 2, -1, 3, 1 };
        mul0 = mul0data;
        
        matrix_mxn< 3, 2 > mul1;
        double mul1data[] = { 3, 1, 2, 1, 1, 0 };
        mul1 = mul1data;
        
        matrix_mxn< 2, 2 > result;
        result.mul( mul0, mul1 );
        
        matrix_mxn< 2, 2 > correct_result;
        double correct_result_data[] = { 5, 1, 4, 2 };
        correct_result = correct_result_data;
        ok = result == correct_result;
    
        if ( ok )
        {
            result = mul0 * mul1;
            ok = result == correct_result;
        }

		log( "matrix multiplication ( mul(), operator* )", ok );

        if ( ! ok )
        {
			std::stringstream error;
			error
				<< "result = M0 * M1 \n"
                << "M0 "        << mul0 
                << "M1 "        << mul1 
                << "result "    << result
                << std::endl;
			log_error( error.str() );
        }
    }

	// matrix addition
	
	if ( ok )
	{
		matrix_mxn< 2, 3 > MplusM = m0 + m0;
	
	}

    #if 0
    // matrix inversion for 2x2 using functor
    
    if ( ok )
    {
        matrix_mxn< 2, 2 > M, M_inverse, M_inverse_correct;
        double Mdata[] = { 1, 3, 4, 2 };
        M = Mdata;

        double M_inverse_correct_data[] = { -0.2, 0.3, 0.4, -0.1 };
        M_inverse_correct = M_inverse_correct_data;
            
        invert_2x2_matrix_functor<> computeInverse;
        computeInverse( M, M_inverse );
        
        ok = M_inverse == M_inverse_correct;
        if ( ! ok )
        {
            log_error( "matrix inverse computation failed, possibly due to precision errors.", true );
            ok = M_inverse.isEqualTo( M_inverse_correct, 1e-15 );
			log( "matrix inversion with reduced precision (tolerance: 1e-15).", ok );
        }
		log( "matrix inversion for 2x2 matrices (functor)", ok );
        if ( ! ok )
        {  
			std::stringstream error;
			error
                << "matrix M " << M 
                << "inverse (computed)" << M_inverse 
                << "inverse (correct)" << M_inverse_correct 
                << std::endl;
			log_error( error.str() );
        }
    }
    #endif
    
    // matrix inversion for 2x2 using matrix_mxm

    if ( ok )
    {
        matrix_mxm< 2 > M, M_inverse, M_inverse_correct;
        double Mdata[] = { 1, 3, 4, 2 };
        M.copyFrom1DimCArray( Mdata );

        double M_inverse_correct_data[] = { -0.2, 0.3, 0.4, -0.1 };
        M_inverse_correct.copyFrom1DimCArray( M_inverse_correct_data );
            
        M.computeInverse( M_inverse );
        ok = M_inverse == M_inverse_correct;
        if ( ! ok )
        {
            log_error( "matrix inverse computation failed, possibly due to precision errors.", true );
            ok = M_inverse.isEqualTo( M_inverse_correct, 1e-15 );
			log( "matrix inversion with reduced precision (tolerance: 1e-15).", ok );
        }
		log( "matrix inversion for 2x2 matrices (specialized template)", ok );
        if ( ! ok )
        {  
			std::stringstream error;
			error
                << "matrix M " << M 
                << "inverse (computed)" << M_inverse 
                << "inverse (correct)" << M_inverse_correct 
                << std::endl;
			log_error( error.str() );
        }
    }
    
    return ok;
}




} // namespace vmml

