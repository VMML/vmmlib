#ifndef __VMML__t3_converter__HPP__
#define __VMML__t3_converter__HPP__

#include "tensor3.hpp"

//TODO: make open file methods and move other write/read methods

namespace vmml
{
	
	template< size_t I1, size_t I2, size_t I3, typename T = float >
	class t3_converter
	{
		
	public:	
		
				
		typedef tensor3< I1, I2, I3, T > t3_t;
		
		template< typename T_convert >
		static void convert_raw( const std::string& dir_, const std::string& in_filename_, const std::string& out_filename_ );

		//header size as bytes
		static void remove_uct_cylinder( const std::string& dir_, 
										const std::string& in_filename_, 
										const std::string& out_filename_,
										const double& sigma_, 
										const size_t header_size_,
										const size_t radius_offset_, 
										int seed_ = 0);

		static void export_to( std::vector< T >& data_ ) ;
		static void import_from( const std::vector< T >& data_ ) ;
		static void write_to_raw( const t3_t& data_, const std::string& dir_, const std::string& filename_ );
		static void read_from_raw( t3_t& data_, const std::string& dir_, const std::string& filename_ ) ;
		static void write_datfile( const std::string& dir_, const std::string& filename_ );
		static void write_to_csv( const t3_t& data_, const std::string& dir_, const std::string& filename_ );

		template< typename TT  >
		static void quantize_to( const std::string& dir_, 
						 const std::string& in_filename_, const std::string& out_filename_,
						 const T& min_value_, const T& max_value_ );

		
	protected:
		
		static void concat_path( const std::string& dir_, const std::string& filename_, std::string& path_ );
		
	}; //end t3_converter
	
	
	
#define VMML_TEMPLATE_STRING        template< size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME     t3_converter< I1, I2, I3, T >

	
	VMML_TEMPLATE_STRING
	template< typename TT  >
	void
	VMML_TEMPLATE_CLASSNAME::quantize_to( const std::string& dir_, 
										 const std::string& in_filename_, const std::string& out_filename_,
										 const T& min_value_, const T& max_value_ )
	{
		std::string path_in_raw = "";
		std::string path_out_raw = "";
		concat_path( dir_, in_filename_, path_in_raw );
		concat_path( dir_, out_filename_, path_out_raw );

		std::ofstream outfile;	
		outfile.open( path_out_raw.c_str() );
		
		std::ifstream infile;
		infile.open( path_in_raw.c_str(), std::ios::in); 
		
		if( infile.is_open() && outfile.is_open() )
		{
			double max_tt_range = double(std::numeric_limits< TT >::max());
			double min_tt_range = double(std::numeric_limits< TT >::min());
			double tt_range = max_tt_range - min_tt_range;
			double t_range = max_value_ - min_value_;
			
			//std::cout << "tt min= " << min_tt_range << ", tt max= " << max_tt_range << ", t min= " << min_value_ << ", t max= " << max_value_ << std::endl;
			//std::cout << "tt range=" << tt_range << ", t range= " << t_range << std::endl;

			T* in_value;
			TT out_value;
			size_t len_in = sizeof(T);
			size_t len_out = sizeof(TT);
			char* data = new char[ len_in ];
			
			for( size_t i3 = 0; i3 < I3; ++i3 )
			{
				for ( size_t i1 = 0; i1 < I1; ++i1 )
				{
					for( size_t i2 = 0; i2 < I2; ++i2 )
					{
						//read value
						infile.read( data, len_in );
						in_value = (T*)&(data[0]);
						
						//Quantize value
						if (std::numeric_limits<TT>::is_signed ) {
							out_value = TT( std::min( std::max( min_tt_range, double(( *in_value * tt_range / t_range ) + 0.5)), max_tt_range ));
						} else {
							out_value = TT( std::min( std::max( min_tt_range, double(((*in_value - min_value_ ) * tt_range / t_range) + 0.5)), max_tt_range ));
						}
						
						//write_value
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
		std::string path_in_raw = "";
		std::string path_out_raw = "";
		concat_path( dir_, in_filename_, path_in_raw );
		concat_path( dir_, out_filename_, path_out_raw );
		
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
	VMML_TEMPLATE_CLASSNAME::remove_uct_cylinder( const std::string& dir_, 
												 const std::string& in_filename_, 
												 const std::string& out_filename_, 
												 const double& sigma_, 
												 const size_t header_size_,
												 const size_t radius_offset_, 
												 int seed_ )
	{
		std::string path_in_raw = "";
		std::string path_out_raw = "";
		concat_path( dir_, in_filename_, path_in_raw );
		concat_path( dir_, out_filename_, path_out_raw );
		
		std::ofstream outfile;	
		outfile.open( path_out_raw.c_str() );
		
		std::ifstream infile;
		infile.open( path_in_raw.c_str(), std::ios::in); 
		
		//for noise adding in outer area
		if ( seed_ >= 0 )
			srand( seed_ );
		
		double length = 0;
		double radius = (I1-1.0)/2.0 - radius_offset_;
		radius *= radius;
		double k1 = 0;
		double k2 = 0;
		double fill_value = 0;
			
		if( infile.is_open() && outfile.is_open() )
		{
			T* in_value;
			T out_value;
			size_t len_val = sizeof(T);
			char* data = new char[ len_val ];
			
			//skip header
			infile.read( data, header_size_ );
			
			//Read/write data
			for( size_t i3 = 0; i3 < I3; ++i3 )
			{
				for ( size_t i1 = 0; i1 < I1; ++i1 )
				{
					k1 = i1 - (I1-1.0)/2.0;
					for( size_t i2 = 0; i2 < I2; ++i2 )
					{
						infile.read( data, len_val );
						in_value = (T*)&(data[0]);
						fill_value = static_cast< T > (*in_value);
						
						//check if value is outside cylinder
						k2 = i2 - (I2-1.0)/2.0;
						length = k1*k1 + k2*k2 ;
						if ( length >= radius )
						{
							fill_value = rand();
							fill_value /= RAND_MAX;
							fill_value *= sigma_;
						}
						
						out_value = static_cast< T > ( fill_value );
						outfile.write( (char*)&(out_value), len_val );
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
