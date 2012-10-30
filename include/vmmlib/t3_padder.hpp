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
	template< size_t B, size_t I, size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_padder
	{ //tensor3 block reader, 
	public:    
		
		typedef tensor3< I1, I2, I3, T > t3_t;
		typedef tensor3< I, I, I, T > t3_pow2_t;
		typedef tensor3< B, B, B, T > t3_block_t;
		typedef tensor3< I1, I2, B, T > t3_cached_t;


		t3_padder( const std::string& dir_, const std::string& filename_ );  
		
		~t3_padder();
		
		void get_data_block( t3_block_t& t3_block_, size_t j1_, size_t j2_, size_t j3_ );
		void read_from_raw( const size_t j3_ );
		
		
	private:
		std::string _path;
		t3_cached_t _cached_data;
		size_t _cache_idx;
		
	}; // class t3_padder
	
	
#define VMML_TEMPLATE_STRING       template< size_t B, size_t I, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME    t3_padder< B, I, I1, I2, I3, T >
	
	
	VMML_TEMPLATE_STRING
	VMML_TEMPLATE_CLASSNAME::t3_padder( const std::string& dir_, const std::string& filename_ )
	: _cache_idx( 0 )
	{
		_path = dir_;
		_path.append( filename_ );
		
		std::ifstream infile;
		infile.open( _path.c_str(), std::ios::in); 
		
		if( ! infile.is_open() )
		{
			VMMLIB_ERROR( "t3_padder - no input file open", VMMLIB_HERE );
		}
		
		read_from_raw( _cache_idx );
	}
	

	
	VMML_TEMPLATE_STRING
	VMML_TEMPLATE_CLASSNAME::~t3_padder( )
	{
	}
	
	VMML_TEMPLATE_STRING
	void
	VMML_TEMPLATE_CLASSNAME::get_data_block( t3_block_t& t3_block_, size_t j1_, size_t j2_, size_t j3_ )
	{
		//check if data at j3_ index was already cached
		
		if ( _cache_idx == j3_ )
		{
			//do nothing
		}
		else
		{
			//if no, read all slices for given J3 range, j3_ until j3_+B
			read_from_raw( j3_ );
		}
		
		// if block is within data area -> return respective block
		if( j1_ < I1 && j2_ < I2 && j3_ < I3 && (j1_ + B) < I1 && (j2_ + B) < I2 && (j3_ + B) < I3 )
		{
			_cached_data.get_sub_tensor3( t3_block_, j1_, j2_ );
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
			//FIXME: check sizes
			//Maybe do it withouth templates
			const size_t k1 = 1;
			const size_t k2 = 1;
			const size_t k3 = 1;
			//std::cout << "subvolume: (" << k1 << "," << k2 << "," << k3 << ")" << std::endl;
			tensor3< k1, k2, k3, T > t3_sub_block;
			_cached_data.get_sub_tensor3( t3_sub_block, j1_, j2_ );
			t3_block_.set_sub_tensor3( t3_sub_block, B-k1, B-k2, B-k3 );
		}

	}
	
	VMML_TEMPLATE_STRING
	void
	VMML_TEMPLATE_CLASSNAME::read_from_raw( const size_t j3_ )
	{
		//FIXME: replace with actual reading
		_cached_data.fill_increasing_values();
	}
	
	
	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
} // namespace vmml

#endif

