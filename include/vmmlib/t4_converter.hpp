#ifndef __VMML__T4_CONVERTER__HPP__
#define __VMML__T4_CONVERTER__HPP__

#include "tensor4.hpp"

namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, size_t I4, typename T = float >
	class t4_converter
	{
		
	public:	
		
		
		typedef tensor3< I1, I2, I3, T > t3_t;
		
		template< typename T_convert >
		static void convert_raw( const std::string& dir_, const std::string& in_filename_, const std::string& out_filename_ );
		
		
		static void export_to( std::vector< T >& data_ ) ;
		static void import_from( const std::vector< T >& data_ ) ;
		
		static void write_to_raw( const t3_t& data_, const std::string& dir_, const std::string& filename_ ); //TODO: DK
		static void read_from_raw( t3_t& data_, const std::string& dir_, const std::string& filename_ ) ; //TODO: DK
		
		static void write_datfile( const std::string& dir_, const std::string& filename_ );
		static void write_to_csv( const t3_t& data_, const std::string& dir_, const std::string& filename_ ); //TODO: DK
		
		
	protected:
		
		static void concat_path( const std::string& dir_, const std::string& filename_, std::string& path_ );
		
	}; //end t4_converter
	
	
	
#define VMML_TEMPLATE_STRING        template< size_t I1, size_t I2, size_t I3, size_t I4, typename T >
#define VMML_TEMPLATE_CLASSNAME     t4_converter< I1, I2, I3, I4, T >
	
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::concat_path( const std::string& dir_, const std::string& filename_, std::string& path_ )
{
	int dir_length = dir_.size() -1;
	int last_separator = dir_.find_last_of( "/");
	path_ = dir_;
	if (last_separator < dir_length ) {
		path_.append( "/" );
	}
	path_.append( filename_ );
	
	//check for format
	if( filename_.find( "raw", filename_.size() -3) == (-1)) {
		path_.append( ".");
		path_.append( "raw" );
	}
}		


VMML_TEMPLATE_STRING
template< typename T_convert >
void
VMML_TEMPLATE_CLASSNAME::convert_raw( const std::string& dir_, const std::string& in_filename_, const std::string& out_filename_ )
{
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::read_from_raw( t3_t& data_, const std::string& dir_, const std::string& filename_ ) 
{	
}	




VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::write_to_raw( const t3_t& data_, const std::string& dir_, const std::string& filename_ ) 
{	
	// TODO: adapt for tensor4
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
