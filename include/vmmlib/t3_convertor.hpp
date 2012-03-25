#ifndef __VMML__T3_CONVERTOR__HPP__
#define __VMML__T3_CONVERTOR__HPP__

#include "tensor3.hpp"



namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_convertor
	{
		
	public:	
		
				
		typedef tensor3< I1, I2, I3, T > t3_t;
		
		template< typename T_convert >
		static void convert_raw( const std::string& dir_, const std::string& in_filename_, const std::string& out_filename_ );
				
		static void export_to( std::vector< T >& data_ ) ;
		static void import_from( const std::vector< T >& data_ ) ;
		static void write_to_raw( const t3_t& data_, const std::string& dir_, const std::string& filename_ );
		static void read_from_raw( t3_t& data_, const std::string& dir_, const std::string& filename_ ) ;
		static void write_datfile( const std::string& dir_, const std::string& filename_ );
		static void write_to_csv( const t3_t& data_, const std::string& dir_, const std::string& filename_ );
		//static void remove_normals_from_raw( const std::string& dir_, const std::string& filename_ ) ;
		//static void remove_uct_cylinder( const size_t radius_offset_, int seed_ = 0 ) ;
		
		
	protected:
		
	}; //end t3_convertor
	
	
	
#define VMML_TEMPLATE_STRING        template< size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME     t3_convertor< I1, I2, I3, T >
	
	
	VMML_TEMPLATE_STRING
	template< typename T_convert >
	void
	VMML_TEMPLATE_CLASSNAME::convert_raw( const std::string& dir_, const std::string& in_filename_, const std::string& out_filename_ )
	{
		int dir_length = dir_.size() -1;
		int last_separator = dir_.find_last_of( "/");
		std::string path_in = dir_;
		std::string path_out = dir_;
		if (last_separator < dir_length ) {
			path_in.append( "/" );
			path_out.append( "/" );
		}
		path_in.append( in_filename_ );
		path_out.append( out_filename_ );
		
		//check for format
		if( in_filename_.find( "raw", in_filename_.size() -3) == (-1)) {
			path_in.append( ".");
			path_in.append( "raw" );
		}
		std::string path_in_raw = path_in;

		//check for format
		if( out_filename_.find( "raw", out_filename_.size() -3) == (-1)) {
			path_out.append( ".");
			path_out.append( "raw" );
		}
		std::string path_out_raw = path_out;
		
		std::ofstream outfile;	
		outfile.open( path_out_raw.c_str() );
		
		std::ifstream infile;
		infile.open( path_in_raw.c_str(), std::ios::in); 
		
		if( infile.is_open() && outfile.is_open() )
		{
			T* in_value;
			T_convert out_value;
			size_t len_in = sizeof(T);
			size_t len_out = sizeof(T_convert);
			char* data = new char[ len_in ];
	
			for( size_t i3 = 0; i3 < I3; ++i3 )
			{
				for ( size_t i1 = 0; i1 < I1; ++i1 )
				{
					for( size_t i2 = 0; i2 < I2; ++i2 )
					{
						infile.read( data, len_in );
						in_value = (T*)&(data[0]);
						out_value = static_cast< T_convert> (*in_value);
						outfile.write( (char*)&(out_value), len_out );
					}
				}
			}	
			
			infile.close();
			outfile.close();
		} else {
			infile.close();
			outfile.close();
			std::cout << "no file open" << std::endl;
		}
	}
	
	
	
	VMML_TEMPLATE_STRING
	void
	VMML_TEMPLATE_CLASSNAME::read_from_raw( t3_t& data_, const std::string& dir_, const std::string& filename_ ) 
	{	
#if FIXME
		int dir_length = dir_.size() -1;
		int last_separator = dir_.find_last_of( "/");
		std::string path = dir_;
		if (last_separator < dir_length ) {
			path.append( "/" );
		}
		path.append( filename_ );
		
		size_t max_file_len = 2147483648u - sizeof(T) ;
		size_t len_data = sizeof(T) * SIZE;
		size_t len_read = 0;
		char* data = new char[ len_data ];
		std::ifstream infile;
		infile.open( path.c_str(), std::ios::in); 
		
		if( infile.is_open())
		{
			iterator  it = begin(),
			it_end = end();
						
			while ( len_data > 0 ) 
			{
				len_read = (len_data % max_file_len ) > 0 ? len_data % max_file_len : len_data;
				len_data -= len_read;
				infile.read( data, len_read );
				
				T* T_ptr = (T*)&(data[0]);
				for( ; (it != it_end) && (len_read > 0); ++it, len_read -= sizeof(T) )
				{
					*it = *T_ptr; ++T_ptr;
				}
			}
			
			delete[] data;
			infile.close();
		} else {
			infile.close();
			std::cout << "no file open" << std::endl;
		}
#endif
		
	}	
	
	
	
	
	VMML_TEMPLATE_STRING
	void
	VMML_TEMPLATE_CLASSNAME::write_to_raw( const t3_t& data_, const std::string& dir_, const std::string& filename_ ) 
	{		
		int dir_length = dir_.size() -1;
		int last_separator = dir_.find_last_of( "/");
		std::string path = dir_;
		if (last_separator < dir_length ) {
			path.append( "/" );
		}
		path.append( filename_ );
		//check for format
		if( filename_.find( "raw", filename_.size() -3) == (-1)) {
			path.append( ".");
			path.append( "raw" );
		}
		std::string path_raw = path;
		
		std::ofstream outfile;	
		outfile.open( path_raw.c_str() );
		if( outfile.is_open() ) {
			size_t len_slice = sizeof(T) * I1 * I2;
			for( size_t index = 0; index < I3; ++index )
			{
				outfile.write( (char*)&(data_.get_frontal_slice_fwd( index )), len_slice );
			}
			outfile.close();
		} else {
			outfile.close();
			std::cout << "no file open" << std::endl;
		}
	}	
	
	
	
	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
	
}//end vmml namespace

#endif
