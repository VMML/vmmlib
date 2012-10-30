/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * t3_padder is a cubic tensor3 data structure, but not all of its content is filled
 * blocks of data can be accessed (read) on demand: if the data is available it is read, if there is no data available it is zero-padded and returned
 * 
 */

#ifndef __VMML__T3_PADDER__HPP__
#define __VMML__T3_PADDER__HPP__

#include <vmmlib/tensor3.hpp>

namespace vmml
{
	// I = cubic virtual data structure, I1-I2: data dimensions where data is available
	template< size_t I, size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_padder
	{ //tensor3 block reader, 
	public:    
		
		typedef tensor3< I1, I2, I3, T > t3_type;
		typedef tensor3< I, I, I, T > t3_pow2_type;


		t3_padder( const std::string& dir_, const std::string& filename_ );  
		
		~t3_padder();
		
		template< size_t B >
		void get_tensor3_block( tensor3< B, B, B, T >& t3_block_, size_t j1_, size_t j2_, size_t j3_ );
		
		
	private:
		std::string _path;
		
	}; // class t3_padder
	
	
#define VMML_TEMPLATE_STRING       template< size_t I, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME    t3_padder< I, I1, I2, I3, T >
	
	
	VMML_TEMPLATE_STRING
	VMML_TEMPLATE_CLASSNAME::t3_padder( const std::string& dir_, const std::string& filename_ )
	{
		_path = dir_;
		_path.append( filename_ );
		
		std::ifstream infile;
		infile.open( _path.c_str(), std::ios::in); 
		
		if( ! infile.is_open() )
		{
			VMMLIB_ERROR( "t3_padder - no input file open", VMMLIB_HERE );
		}
			
	}
	

	
	VMML_TEMPLATE_STRING
	VMML_TEMPLATE_CLASSNAME::~t3_padder( )
	{
	}
	
	VMML_TEMPLATE_STRING
	template< size_t B >
	void
	VMML_TEMPLATE_CLASSNAME::get_tensor3_block( tensor3< B, B, B, T >& t3_block_, size_t j1_, size_t j2_, size_t j3_ )
	{
		// if block is within data area -> return respective block
		if( j1_ < I1 && j2_ < I2 && j3_ < I3 )
		{
			
		}
		// if block is outside data area -> return zero-padded block
		else if( j1_ >= I1 && j2_ >= I2 && j3_ >= I3 )
		{
			t3_block_.zero();
		}
		// if block is partly outside data area -> return data block zero-padded
		else 
		{
			t3_block_.zero();
		}

	}
	
	
	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
} // namespace vmml

#endif

