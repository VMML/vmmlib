#include "matrix_test.hpp"

#include <sstream>

namespace vmml
{

bool
matrix_test::run()
{
    bool ok = true;
    
    matrix< 2, 3 > m0;
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
				ok = m0.at( rowIndex, colIndex ) == tmp;
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
	
    matrix< 3, 2 > m1;
    m1.copyFrom1DimCArray( data, false );

    // tests copyFrom1DimCArray function with col_by_col data
	ok = true;
	{
		for( size_t tmp = 1, colIndex = 0; ok && colIndex < 2; ++colIndex )
		{
			for( size_t rowIndex = 0; ok && rowIndex < 3; ++rowIndex, ++tmp )
			{
				ok = m1.at( rowIndex, colIndex ) == tmp;
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
		matrix< 2, 3 > m0_copy;
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

    // test operator[]
	ok = true;
	{
		m0 = data;
		ok = m0[ 1 ][ 1 ] == 5;
		if ( ok )
			ok = m0[ 1 ][ 1 ] == m0.at( 1, 1 );
		if ( ok ) 
		{
			m0[ 1 ][ 2 ] = 23;
			ok = m0.at( 1, 2 ) == 23;
		}
		
		log( "operator[]", ok );
		if ( ! ok )
		{
			std::stringstream error;
			error 
				<< " m0 " << m0 
				<< " m0[ 1 ][ 1 ] " << m0[ 1 ][ 1 ] 
				<< " m0.at( 1, 1 ) " << m0.at( 1, 1 )
				<< std::endl;
			log_error( error.str() );
		}
	}

    // test getRow/setRow/getColumn/setColumn
	ok = true;
	{
		matrix< 2, 3 > M;
		double Mdata[] = { 1, 2, 3, 4, 5, 6 };
		M = Mdata;
		matrix< 1, 3 > M_row = M.getRow( 1 );
		matrix< 2, 1 > M_column = M.getColumn( 2 );
		
		for( size_t column = 0; ok && column < 3; ++column )
		{
			ok = M.at( 1, column ) == M_row.at( 0, column );
		}

		for( size_t row = 0; ok && row < 2; ++row )
		{
			ok = M.at( row, 2 ) == M_column.at( row, 0 );
		}
		
		double Mdata_row[] = { 3, 2, 5, 4, 5, 6 };
		matrix< 2, 3 > Mr;
		Mr = Mdata_row;
		
		M = Mdata;
		M_row.at( 0, 0 ) = 3;
		M_row.at( 0, 1 ) = 2;
		M_row.at( 0, 2 ) = 5;
		M.setRow( 0, M_row );
		for( size_t column = 0; ok && column < 3; ++column )
		{
			ok = M == Mr;
		}
		
		double Mdata_column[] = { 1, 5, 3, 4, 2, 6 };
		matrix< 2, 3 > Mc;
		Mc = Mdata_column;
		M = Mdata;
		M_column.at( 0, 0 ) = 5;
		M_column.at( 1, 0 ) = 2;
		M.setColumn( 1, M_column );
		for( size_t row = 0; ok && row < 2; ++row )
		{
			ok = M == Mc;
		}

		log( "getRow/setRow/getColumn/setColumn", ok );
		if ( ! ok )
		{
			std::stringstream error;
			error 
				<< "M " << M 
				<< "M_row " << M_row 
				<< "Mr " << Mr 
				<< "M_column " << M_column 
				<< "Mc " << Mc 
				<< std::endl;
			log_error( error.str() );
		}
	}

	
    // test transpose functionality 
	ok = true;
    {
		m0 = data;
        matrix< 3, 2 > m0t = m0.getTransposed();
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
        matrix< 2, 3 > mul0;
        double mul0data[] = { 1, 0, 2, -1, 3, 1 };
        mul0 = mul0data;
        
        matrix< 3, 2 > mul1;
        double mul1data[] = { 3, 1, 2, 1, 1, 0 };
        mul1 = mul1data;
        
        matrix< 2, 2 > result;
        result.multiply( mul0, mul1 );
        
        matrix< 2, 2 > correct_result;
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


	#ifdef VMMLIB_SAFE_ACCESSORS
	ok = true;
	{
		matrix< 3, 2 > m;
		try
		{
			m.at( 3, 2 );
			ok = false;
		}
		catch(...)
		{}
		try
		{
			if ( ok ) m.at( 1, 1 );
		}
		catch(...)
		{
			ok = false;
		}
		try
		{
			if ( ok ) 
			{
				m[ 3 ][ 2 ];
				ok = false;
			}
		}
		catch(...)
		{}
		log( "safe accessors (debug option)", ok );
	}
	#endif

	// matrix addition
	ok = true;
	{
		// FIXME
		matrix< 2, 3 > MplusM = m0 + m0;
	
	}

	// getSubMatrix
	ok = true;
	{
		matrix< 1, 2 > m;
		m = m0.getSubMatrix< 1, 2 >( 1, 1 );
		ok = m.at( 0, 0 ) == m0.at( 1, 1 );
		if ( ok ) ok = m.at( 0, 1 ) == m0.at( 1, 2 );
		log( "getSubMatrix()", ok );
        if ( ! ok )
        {
			std::stringstream error;
			error
                << "m0 " << m0 
                << "m " << m 
				<< "m0.at( 1, 1 ): " << m0.at( 1, 1 ) << "\n"
                << "m.at( 0, 0 ): " << m.at( 0, 0 ) << "\n"
				<< "m0.at( 1, 2 ): " << m0.at( 1, 2 )<< "\n"
                << "m.at( 0, 1 ): " << m.at( 0, 1 )<< "\n"
                << std::endl;
			log_error( error.str() );
        }
	}


    // matrix inversion for 2x2
	ok = true;
    {
        matrix< 2, 2 > M, M_inverse, M_inverse_correct;
        double Mdata[] = { 1, 3, 4, 2 };
        M.copyFrom1DimCArray( Mdata );

        double M_inverse_correct_data[] = { -0.2, 0.3, 0.4, -0.1 };
        M_inverse_correct.copyFrom1DimCArray( M_inverse_correct_data );

        M.computeInverse( M_inverse );
        
        ok = M_inverse == M_inverse_correct;
		log( "matrix inversion for 2x2 matrices, maximum precision", ok, true );
        if ( ! ok )
        {
            ok = M_inverse.isEqualTo( M_inverse_correct, 1e-15 );
			log( "matrix inversion 2x2 with reduced precision (tolerance: 1e-15).", ok );
        }
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


    // matrix inversion for 3x3 
    if ( ok )
    {
        matrix< 3, 3 > M, M_inverse, M_inverse_correct;
        double Mdata[] = { 8, 1, 6, 3, 5, 7, 4, 9, 2 };
        M.copyFrom1DimCArray( Mdata );

        double M_inverse_correct_data[] = { 0.1472, -0.1444, 0.0639, 
            -0.0611, 0.0222, 0.1056, -0.0194, 0.1889, -0.1028 };

        M_inverse_correct.copyFrom1DimCArray( M_inverse_correct_data );

        M.computeInverse( M_inverse );
        
        ok = M_inverse == M_inverse_correct;
		log( "matrix inversion for 3x3 matrices, maximum precision", ok, true );
        if ( ! ok )
        {
            ok = M_inverse.isEqualTo( M_inverse_correct, 1e-15 );
			log( "matrix inversion 3x3 with reduced precision (tolerance: 1e-15).", ok );
        }
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


    // matrix inversion for 4x4 
    if ( ok )
    {
        matrix< 4, 4 > M, M_inverse, M_inverse_correct;
       double Mdata[] = { 17., 24., 1., 8., 23., 5., 7., 14.,
                 4., 6., 13., 20., 10., 12., 19., 21. };
        M.copyFrom1DimCArray( Mdata );

        double M_inverse_correct_data[] = { -5.780346820809248e-03, 4.962205424633170e-02, -4.811027123165852e-02, 1.493997332147622e-02, 
                4.277456647398844e-02, -3.797243219208537e-02, -1.013783903957314e-02, 1.867496665184526e-02, 
                -3.930635838150288e-02, -1.333926189417519e-02, -1.333036905291240e-01, 1.508225878168074e-01, 
                1.387283236994219e-02, 1.013783903957314e-02, 1.493108048021343e-01, -1.066251667407737e-01 };

        M_inverse_correct.copyFrom1DimCArray( M_inverse_correct_data );

        M.computeInverse( M_inverse );
        
        ok = M_inverse == M_inverse_correct;
		log( "matrix inversion for 4x4 matrices, maximum precision", ok, true );
        if ( ! ok )
        {
            ok = M_inverse.isEqualTo( M_inverse_correct, 1e-15 );
			log( "matrix inversion 4x4 with reduced precision (tolerance: 1e-15).", ok );
        }
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

