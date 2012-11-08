/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 *
 * t3_virtual_padder is a cubic tensor3 data structure, but not all of its content is filled
 * blocks of data can be accessed (read) on demand: if the data is available it is read, if there is no data available it is zero-padded and returned
 * 
 */

#ifndef __VMML__T3_VIRTUAL_PADDER__HPP__
#define __VMML__T3_VIRTUAL_PADDER__HPP__

#include <vmmlib/tensor3.hpp>

namespace vmml
{
	// I = cubic virtual data structure, I1-I2: data dimensions where data is available
	template< size_t B, size_t I, size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_virtual_padder
	{ //tensor3 block reader, 
	public:    
		
		typedef tensor3< I1, I2, I3, T > t3_t;
		typedef tensor3< I, I, I, T > t3_pow2_t;
		typedef tensor3< B, B, B, T > t3_block_t;
		typedef tensor3< I1, I2, B, T > t3_cached_t;
		
		typedef typename vmml::tensor3_iterator< t3_cached_t >	         iterator;


		t3_virtual_padder( const std::string& dir_, const std::string& filename_ );  
		
		~t3_virtual_padder();
		
		void get_data_block( t3_block_t& t3_block_, size_t j1_, size_t j2_, size_t j3_ );
		void read_from_raw();

		
	private:
		std::string _path;
		t3_cached_t _cached_data;
		size_t _cache_idx;
		
	}; // class t3_virtual_padder
	
	
#define VMML_TEMPLATE_STRING       template< size_t B, size_t I, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME    t3_virtual_padder< B, I, I1, I2, I3, T >
	
	
	VMML_TEMPLATE_STRING
	VMML_TEMPLATE_CLASSNAME::t3_virtual_padder( const std::string& dir_, const std::string& filename_ )
	: _cache_idx( 0 )
	{
		_path = dir_;
		_path.append( filename_ );
		
		std::ifstream infile;
		infile.open( _path.c_str(), std::ios::in); 
		
		if( ! infile.is_open() )
		{
			VMMLIB_ERROR( "t3_virtual_padder - no input file open", VMMLIB_HERE );
		}
		
		read_from_raw();
	}
	

	
	VMML_TEMPLATE_STRING
	VMML_TEMPLATE_CLASSNAME::~t3_virtual_padder( )
	{
	}
	
	VMML_TEMPLATE_STRING
	void
	VMML_TEMPLATE_CLASSNAME::get_data_block( t3_block_t& t3_block_, size_t j1_, size_t j2_, size_t j3_ )
	{
		
		// if block is outside data area -> return zero-padded block
		if( j1_ >= I1 && j2_ >= I2 && j3_ >= I3 )
		{
			t3_block_.zero();
		}
		else {
			//check if data at j3_ index was already cached
			if ( _cache_idx == j3_ )
			{
				//do nothing
			}
			else
			{
				//if no, read all slices for given J3 range, j3_ until j3_+B
				_cache_idx = j3_;
				read_from_raw();
			}
			
			// if block is within data area -> return respective block
			if( j1_ < I1 && j2_ < I2 && j3_ < I3 && (j1_ + B) < I1 && (j2_ + B) < I2 && (j3_ + B) < I3 )
			{
				_cached_data.get_sub_tensor3( t3_block_, j1_, j2_ );
			}
			// if block is partly outside data area -> return data block zero-padded
			else 
			{
				t3_block_.zero();
				const size_t len1 = I1 - j1_ ;
				const size_t len2 = I2 - j2_ ;
				const size_t len3 = I3 - j3_ ;
				
				size_t block_data_len = len1 * len2 * len3 * sizeof(T);
				char* data = new char[ block_data_len ];
				
				//std::cout << "cached data:\n" << _cached_data << std::endl;
				_cached_data.get_sub_tensor3( data, j1_, j1_+len1-1, j2_, j2_+len2-1, 0, len3-1 );
				
				//assume that block size is small
				
				T* t_ptr = (T*)&(data[0]);
				for ( size_t k3 = 0; k3 < len3; ++k3)
				{
					for ( size_t k1 = 0; k1 < len1; ++k1)
					{
						for ( size_t k2 = 0; k2 < len2; ++k2)
						{
							t3_block_.at( k1, k2, k3 ) = *t_ptr;
							//std::cout << "ptr_val = " << int(*t_ptr) << std::endl;
							++t_ptr;
						}
					}
				}
				delete[] data;
			}
		}

	}
	
	VMML_TEMPLATE_STRING
	void
	VMML_TEMPLATE_CLASSNAME::read_from_raw()
	{
		_cached_data.zero();
		
		size_t max_file_len = 2147483648u - sizeof(T) ;
		size_t len3 = (_cache_idx < (I3 - B)) ? B : (I3 - _cache_idx);
		size_t cache_data_len = I1 * I2 * len3 * sizeof(T);
		size_t len_read = 0;
		size_t start_idx = _cache_idx * I1 * I2 * sizeof(T);
		char* data = new char[ cache_data_len ];
		std::ifstream infile;
		infile.open( _path.c_str(), std::ios::in);
		
		if( infile.is_open())
		{
			infile.seekg( start_idx, std::ios::beg );
			if( infile.tellg() != start_idx )
			{
				std::cout << "Can't proceed to the offset: " << start_idx << " to read file: " << _path.c_str() << std::endl;
				infile.close();
			}
			
			infile.read( &data[0], cache_data_len );
			
			iterator  it = _cached_data.begin(),
			it_end = _cached_data.end();
			
			while ( cache_data_len > 0 )
			{
				len_read = (cache_data_len % max_file_len ) > 0 ? cache_data_len % max_file_len : cache_data_len;
				cache_data_len -= len_read;
				
				T* T_ptr = (T*)&(data[0]);
				for( ; (it != it_end) && (len_read > 0); ++it, len_read -= sizeof(T) )
				{
					*it = *T_ptr; ++T_ptr;
				}
			}
			
			delete[] data;
			infile.close();
		} else {
			std::cout << "no file open" << std::endl;
		}
		//std::cout << "read file at z index: " << _cache_idx << "\n" << _cached_data << std::endl;
		
	}
	
	
	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
} // namespace vmml

#endif

