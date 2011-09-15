/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 * @author Jonas Boesch
 *
 * The tucker3 tensor class is consists of the same components (core tensor, basis matrices u1-u3) as the tucker3 model described in:
 * - Tucker, 1966: Some mathematical notes on three-mode factor analysis, Psychometrika.
 * - De Lathauwer, De Moor, Vandewalle, 2000a: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - De Lathauwer, De Moor, Vandewalle, 2000b: On the Best rank-1 and Rank-(R_1, R_2, ..., R_N) Approximation and Applications of Higher-Order Tensors, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * 
 */

#ifndef __VMML__TUCKER3_TENSOR__HPP__
#define __VMML__TUCKER3_TENSOR__HPP__

#define CODE_ALL_U_MIN_MAX 0
#define CORE_RANGE 127


#include <vmmlib/t3_hosvd.hpp>
#include <vmmlib/t3_hooi.hpp>

namespace vmml
{
	
	template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T_value = float, typename T_coeff = double >
	class tucker3_tensor
	{
public:    
	
	typedef float T_internal;	
		
	typedef tucker3_tensor< R1, R2, R3, I1, I2, I3, T_value, T_coeff > tucker3_type;
		
	typedef tensor3< I1, I2, I3, T_value > t3_type;
	typedef typename t3_type::iterator t3_iterator;
	typedef typename t3_type::const_iterator t3_const_iterator;
		
	typedef tensor3< I1, I2, I3, T_coeff > t3_coeff_type;
	typedef typename t3_coeff_type::iterator t3_coeff_iterator;
	typedef typename t3_coeff_type::const_iterator t3_coeff_const_iterator;
		
	typedef tensor3< R1, R2, R3, T_coeff > t3_core_type;
	typedef typename t3_core_type::iterator t3_core_iterator;
	typedef typename t3_core_type::const_iterator t3_core_const_iterator;
		
	typedef matrix< R1, R2, T_coeff >        front_core_slice_type; //fwd: forward cylcling (after kiers et al., 2000)

	typedef matrix< I1, R1, T_coeff > u1_type;
	typedef typename u1_type::iterator u1_iterator;
	typedef typename u1_type::const_iterator u1_const_iterator;

	typedef matrix< I2, R2, T_coeff > u2_type;
	typedef typename u2_type::iterator u2_iterator;
	typedef typename u2_type::const_iterator u2_const_iterator;
	
	typedef matrix< I3, R3, T_coeff > u3_type;
	typedef typename u3_type::iterator u3_iterator;
	typedef typename u3_type::const_iterator u3_const_iterator;
	
	typedef tensor3< I1, I2, I3, T_internal > t3_comp_type;
	typedef typename t3_comp_type::iterator t3_comp_iterator;
	typedef typename t3_comp_type::const_iterator t3_comp_const_iterator;
		
	typedef tensor3< R1, R2, R3, T_internal > t3_core_comp_type;
	typedef matrix< I1, R1, T_internal > u1_comp_type;
	typedef matrix< I2, R2, T_internal > u2_comp_type;
	typedef matrix< I3, R3, T_internal > u3_comp_type;
	
	static const size_t SIZE = R1*R2*R3 + I1*R1 + I2*R2 + I3*R3;
		
	tucker3_tensor();
	tucker3_tensor( t3_core_type& core );
	tucker3_tensor( t3_core_type& core, u1_type& U1, u2_type& U2, u3_type& U3 );
	tucker3_tensor( const tucker3_type& other );
	~tucker3_tensor();
		
	void enable_quantify_coeff() { _is_quantify_coeff = true; };
	void disable_quantify_coeff() { _is_quantify_coeff = false; } ;
	void enable_quantify_hot() { _is_quantify_hot = true; _is_quantify_coeff = true; _is_quantify_log = false; _is_quantify_linear = false;};
	void disable_quantify_hot() { _is_quantify_hot = false; _is_quantify_coeff = false;} ;
	void enable_quantify_linear() { _is_quantify_linear = true; _is_quantify_coeff = true; _is_quantify_hot = false;};
	void disable_quantify_linear() { _is_quantify_linear = false; _is_quantify_coeff = false;} ;
	void enable_quantify_log() { _is_quantify_log = true; _is_quantify_coeff = true; _is_quantify_hot = false;};
	void disable_quantify_log() { _is_quantify_log = false; _is_quantify_coeff = false;} ;
		
	tensor3< R1, R2, R3, char> get_core_signs() { return _signs; };
	void set_core_signs(	const tensor3< R1, R2, R3, char> signs_ ) { _signs = signs_; } ;
		
	void set_core( t3_core_type& core )  { _core = t3_core_type( core ); _core_comp.cast_from( core ); } ;
	void set_u1( u1_type& U1 ) { *_u1 = U1; _u1_comp->cast_from( U1 ); } ;
	void set_u2( u2_type& U2 ) { *_u2 = U2; _u1_comp->cast_from( U2 ); } ;
	void set_u3( u3_type& U3 ) { *_u3 = U3; _u1_comp->cast_from( U3 ); } ;
	
	void get_core( t3_core_type& data_ ) const { data_ = _core; } ;
	void get_u1( u1_type& U1 ) const { U1 = *_u1; } ;
	void get_u2( u2_type& U2 ) const { U2 = *_u2; } ;
	void get_u3( u3_type& U3 ) const { U3 = *_u3; } ;
		
	void set_core_comp( t3_core_comp_type& core )  { _core_comp = t3_core_comp_type( core ); _core.cast_from( _core_comp ); } ;
	void set_u1_comp( u1_comp_type& U1 ) { *_u1_comp = U1; _u1->cast_from( U1 ); } ;
	void set_u2_comp( u2_comp_type& U2 ) { *_u2_comp = U2; _u1->cast_from( U2 ); } ;
	void set_u3_comp( u3_comp_type& U3 ) { *_u3_comp = U3; _u1->cast_from( U3 ); } ;

	void get_core_comp( t3_core_comp_type& data_ ) const { data_ = _core_comp; } ;
	void get_u1_comp( u1_comp_type& U1 ) const { U1 = *_u1_comp; } ;
	void get_u2_comp( u2_comp_type& U2 ) const { U2 = *_u2_comp; } ;
	void get_u3_comp( u3_comp_type& U3 ) const { U3 = *_u3_comp; } ;
		
	template< typename T >
	void export_to( std::vector< T >& data_ );
	template< typename T >
	void import_from( const std::vector< T >& data_ );
		
	//previous version, but works only with 16bit quantization
	void export_quantized_to(  std::vector<unsigned char>& data_out_  );
	void import_quantized_from( const std::vector<unsigned char>& data_in_ );
		
	//use this version, works with a better quantization for the core tensor:
	//logarithmic quantization and separate high energy core vale
    //suitable for voxelwise reconstruction
	void export_hot_quantized_to(  std::vector<unsigned char>& data_out_  );
	void import_hot_quantized_from( const std::vector<unsigned char>& data_in_ );
		
	//use this version for the ttm export/import (core: backward cyclic), without plain hot value 
	void export_ttm_quantized_to(  std::vector<unsigned char>& data_out_  );
	void import_ttm_quantized_from( const std::vector<unsigned char>& data_in_ );
		
	//get number of nonzeros for tensor decomposition
	size_t nnz() const;
	size_t nnz( const T_value& threshold ) const;	
	size_t nnz_core() const;
	size_t size_core() const;
	size_t size() const { return SIZE; } ;
	
	void threshold_core( const size_t& nnz_core_, size_t& nnz_core_is_ ); 
	void threshold_core( const T_coeff& threshold_value_, size_t& nnz_core_ ); 
	void reconstruct( t3_type& data_,
					 const T_internal& u_min_, const T_internal& u_max_,
					 const T_internal& core_min_, const T_internal& core_max_ ); 
	void reconstruct( t3_type& data_, 
					 const T_internal& u1_min_, const T_internal& u1_max_,
					 const T_internal& u2_min_, const T_internal& u2_max_,
					 const T_internal& u3_min_, const T_internal& u3_max_,
					 const T_internal& core_min_, const T_internal& core_max_ ); 
	void reconstruct( t3_type& data_ ); 
		
	void decompose( const t3_type& data_ );
	void decompose( const t3_type& data_, 
				   T_internal& u1_min_, T_internal& u1_max_,
				   T_internal& u2_min_, T_internal& u2_max_,
				   T_internal& u3_min_, T_internal& u3_max_,
				   T_internal& core_min_, T_internal& core_max_ ); 
	void decompose( const t3_type& data_, 
					   T_internal& u_min_, T_internal& u_max_,
					   T_internal& core_min_, T_internal& core_max_ ); 
		
	void tucker_als( const t3_type& data_ );	
		
		
	template< size_t K1, size_t K2, size_t K3>
	void reduce_ranks( const tucker3_tensor< K1, K2, K3, I1, I2, I3, T_value, T_coeff >& other ); //call TuckerJI.reduce_ranks(TuckerKI) K1 -> R1, K2 -> R2, K3 -> R3

	template< size_t K1, size_t K2, size_t K3>
	void subsampling( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor  );

	template< size_t K1, size_t K2, size_t K3>
	void subsampling_on_average( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor  );

	template< size_t K1, size_t K2, size_t K3>
	void region_of_interest( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, 
                                 const size_t& start_index1, const size_t& end_index1, 
                                 const size_t& start_index2, const size_t& end_index2, 
                                 const size_t& start_index3, const size_t& end_index3);

	friend std::ostream& operator << ( std::ostream& os, const tucker3_type& t3 )
	{
		t3_core_type core; t3.get_core( core );
		u1_type* u1 = new u1_type; t3.get_u1( *u1 );
		u2_type* u2 = new u2_type; t3.get_u2( *u2 );
		u3_type* u3 = new u3_type; t3.get_u3( *u3 );
		
		os << "U1: " << std::endl << *u1 << std::endl
		<< "U2: " << std::endl << *u2 << std::endl
		<< "U3: " << std::endl << *u3 << std::endl
		<< "core: " << std::endl << core << std::endl;
		
		delete u1;
		delete u2;
		delete u3;
		return os;
	}
		
		
		
protected:
		tucker3_type operator=( const tucker3_type& other ) { return (*this); };
        		
private:
        
        void cast_members();
        void cast_comp_members();
        void quantize_basis_matrices( T_internal& u_min_, T_internal& u_max_ );
        void quantize_basis_matrices( T_internal& u1_min_, T_internal& u1_max_, T_internal& u2_min_, T_internal& u2_max_, T_internal& u3_min_, T_internal& u3_max_ );
        void quantize_core( T_internal& core_min_, T_internal& core_max_ );
        void dequantize_basis_matrices( const T_internal& u1_min_, const T_internal& u1_max_, const T_internal& u2_min_, const T_internal& u2_max_, const T_internal& u3_min_, const T_internal& u3_max_ );
        void dequantize_core( const T_internal& core_min_, const T_internal& core_max_ );
		
        //t3_core_type* _core ;
        u1_type* _u1 ;
        u2_type* _u2 ;
        u3_type* _u3 ;
		t3_core_type _core ;
		
		//used only internally for computations to have a higher precision
        t3_core_comp_type _core_comp ;
        u1_comp_type* _u1_comp ;
        u2_comp_type* _u2_comp ;
        u3_comp_type* _u3_comp ;
		
		T_internal _hottest_core_value;
		tensor3< R1, R2, R3, char> _signs;
		
		bool _is_quantify_coeff; 
		bool _is_quantify_hot; 
		bool _is_quantify_log; 
		bool _is_quantify_linear; 
		
}; // class tucker3_tensor


#define VMML_TEMPLATE_STRING        template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T_value, typename T_coeff >
#define VMML_TEMPLATE_CLASSNAME     tucker3_tensor< R1, R2, R3, I1, I2, I3, T_value, T_coeff >


VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tucker3_tensor( )
	: _is_quantify_coeff( false ), _is_quantify_hot( false ), _hottest_core_value( 0 )
	, _is_quantify_linear( false ), _is_quantify_log( false )
{
	_core.zero();
	_u1 = new u1_type(); _u1->zero();
	_u2 = new u2_type(); _u2->zero();
	_u3 = new u3_type(); _u3->zero();	 
	_core_comp.zero();
	_u1_comp = new u1_comp_type(); _u1_comp->zero();
	_u2_comp = new u2_comp_type(); _u2_comp->zero();
	_u3_comp = new u3_comp_type(); _u3_comp->zero();	
	
	_signs.zero();
}
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tucker3_tensor( t3_core_type& core )
	: _is_quantify_coeff( false ), _is_quantify_hot( false ), _hottest_core_value( 0 )
	, _is_quantify_linear( false ), _is_quantify_log( false )
{
	_core = core;
	_u1 = new u1_type(); _u1->zero();
	_u2 = new u2_type(); _u2->zero();
	_u3 = new u3_type(); _u3->zero();	
	_u1_comp = new u1_comp_type(); _u1_comp->zero();
	_u2_comp = new u2_comp_type(); _u2_comp->zero();
	_u3_comp = new u3_comp_type(); _u3_comp->zero();	
	_core_comp.cast_from( core );
	
	_signs.zero();
}

VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tucker3_tensor( t3_core_type& core, u1_type& U1, u2_type& U2, u3_type& U3 )
	: _is_quantify_coeff( false ), _is_quantify_hot( false ), _hottest_core_value( 0 )
	, _is_quantify_linear( false ), _is_quantify_log( false )
{
	_core = core;
	_u1 = new u1_type( U1 );
	_u2 = new u2_type( U2 );
	_u3 = new u3_type( U3 );
	_u1_comp = new u1_comp_type(); 
	_u2_comp = new u2_comp_type(); 
	_u3_comp = new u3_comp_type(); 	
	cast_comp_members();
	
	_signs.zero();
}

VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::tucker3_tensor( const tucker3_type& other )
	: _is_quantify_coeff( false ), _is_quantify_hot( false ), _hottest_core_value( 0 )
	, _is_quantify_linear( false ), _is_quantify_log( false )
{
	_u1 = new u1_type();
	_u2 = new u2_type();
	_u3 = new u3_type();
	_u1_comp = new u1_comp_type(); 
	_u2_comp = new u2_comp_type(); 
	_u3_comp = new u3_comp_type(); 	
	
	other.get_core( _core );
	other.get_u1( *_u1 );
	other.get_u2( *_u2 );
	other.get_u3( *_u3 );

	cast_comp_members();
	
	_signs.zero();
}
		
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::cast_members()
{
	_u1->cast_from( *_u1_comp );
	_u2->cast_from( *_u2_comp );
	_u3->cast_from( *_u3_comp );	
	_core.cast_from( _core_comp);
}
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::cast_comp_members()
{
	_u1_comp->cast_from( *_u1 );
	_u2_comp->cast_from( *_u2 );
	_u3_comp->cast_from( *_u3 );	
	_core_comp.cast_from( _core);
}

	
VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::nnz_core() const
{	
	return _core_comp.nnz();
}

VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::size_core() const
{	
	return _core_comp.size();
}
	
	
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::quantize_basis_matrices(T_internal& u1_min_, T_internal& u1_max_,
												 T_internal& u2_min_, T_internal& u2_max_,
												 T_internal& u3_min_, T_internal& u3_max_ )
{
	_u1_comp->quantize( *_u1, u1_min_, u1_max_ );
	_u2_comp->quantize( *_u2, u2_min_, u2_max_ );
	_u3_comp->quantize( *_u3, u3_min_, u3_max_ );	
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::quantize_basis_matrices(T_internal& u_min_, T_internal& u_max_)
{
	u_min_ = _u1_comp->get_min();
	T_internal u2_min = _u2_comp->get_min();
	T_internal u3_min = _u3_comp->get_min();
	
	if ( u2_min < u_min_) {
		u_min_  = u2_min;
	}
	if ( u3_min < u_min_) {
		u_min_  = u3_min;
	}
	
	u_max_ = _u1_comp->get_max();
	T_internal u2_max = _u2_comp->get_max();
	T_internal u3_max = _u3_comp->get_max();
	
	if ( u2_max > u_max_ ) {
		u_max_  = u2_max;
	}
	if ( u3_max > u_max_ ) {
		u_max_  = u3_max;
	}
		
	_u1_comp->quantize_to( *_u1, u_min_, u_max_ );
	_u2_comp->quantize_to( *_u2, u_min_, u_max_ );
	_u3_comp->quantize_to( *_u3, u_min_, u_max_ );	
	
#if 0
	std::cout << "quantized (1u): " << std::endl << "u1-u3: " << std::endl
	<< *_u1 << std::endl << *_u1_comp << std::endl
	<< *_u2 << std::endl << *_u2_comp << std::endl
	<< *_u3 << std::endl << *_u3_comp << std::endl
	<< " core " << std::endl
	<< _core << std::endl
	<< " core_comp " << std::endl
	<< _core_comp << std::endl;
#endif
}	

	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::quantize_core( T_internal& core_min_, T_internal& core_max_ )
{
	if ( _is_quantify_hot ) {
		_hottest_core_value = _core_comp.at(0,0,0);
		_core_comp.at( 0, 0, 0 ) = 0;		
		_core_comp.quantize_log( _core, _signs, core_min_, core_max_, T_coeff(CORE_RANGE) );
	} else if ( _is_quantify_linear ) {
		_core_comp.quantize( _core, core_min_, core_max_ );
	} else if ( _is_quantify_log ) {
		_core_comp.quantize_log( _core, _signs, core_min_, core_max_, T_coeff(CORE_RANGE) );
	} else {
		_core_comp.quantize( _core, core_min_, core_max_ );
		std::cout << "quant.method not specified" << std::endl;
	}
}	


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::dequantize_basis_matrices( const T_internal& u1_min_, const T_internal& u1_max_, 
													 const T_internal& u2_min_, const T_internal& u2_max_, 
													 const T_internal& u3_min_, const T_internal& u3_max_ )
{
	_u1->dequantize( *_u1_comp, u1_min_, u1_max_ );
	_u2->dequantize( *_u2_comp, u2_min_, u2_max_ );
	_u3->dequantize( *_u3_comp, u3_min_, u3_max_ );	
}	

VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::dequantize_core( const T_internal& core_min_, const T_internal& core_max_ )
{
	if ( _is_quantify_hot ) {
		_core.dequantize_log( _core_comp, _signs, core_min_, core_max_ );
		_core.at(0,0,0) = _hottest_core_value;
		_core_comp.at(0,0,0) = _hottest_core_value;
	} else if ( _is_quantify_linear ) {
		_core.dequantize( _core_comp, core_min_, core_max_ );
	} else if ( _is_quantify_log ) {
		_core.dequantize_log( _core_comp, _signs, core_min_, core_max_ );
	} else {
		_core.dequantize( _core_comp, core_min_, core_max_ );
	}	
}		
	
VMML_TEMPLATE_STRING
VMML_TEMPLATE_CLASSNAME::~tucker3_tensor( )
{
	delete _u1;
	delete _u2;
	delete _u3;
	delete _u1_comp;
	delete _u2_comp;
	delete _u3_comp;
}
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruct( t3_type& data_, 
									 const T_internal& u_min_, const T_internal& u_max_, 
									 const T_internal& core_min_, const T_internal& core_max_ )
{
	dequantize_basis_matrices( u_min_, u_max_, u_min_, u_max_, u_min_, u_max_ );
	dequantize_core( core_min_, core_max_ );
	
#if 0
	std::cout << "dequantized (1u): " << std::endl << "u1-u3: " << std::endl
	<< *_u1 << std::endl << *_u1_comp << std::endl
	<< *_u2 << std::endl << *_u2_comp << std::endl
	<< *_u3 << std::endl << *_u3_comp << std::endl
	<< " core " << std::endl
	<< _core << std::endl
	<< " core_comp " << std::endl
	<< _core_comp << std::endl;
#endif
	
	reconstruct( data_ );
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruct( t3_type& data_, 
									 const T_internal& u1_min_, const T_internal& u1_max_,
									 const T_internal& u2_min_, const T_internal& u2_max_,
									 const T_internal& u3_min_, const T_internal& u3_max_,
									 const T_internal& core_min_, const T_internal& core_max_ )
{
	dequantize_basis_matrices( u1_min_, u1_max_, u2_min_, u2_max_, u3_min_, u3_max_ );
	dequantize_core( core_min_, core_max_ );
	
    reconstruct( data_ );
}
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::reconstruct( t3_type& data_ )
{
	t3_comp_type data;
	data.cast_from( data_ );
	data.full_tensor3_matrix_multiplication( _core_comp, *_u1_comp, *_u2_comp, *_u3_comp );
	
	//convert reconstructed data, which is in type T_internal (double, float) to T_value (uint8 or uint16)
	if( (sizeof(T_value) == 1) || (sizeof(T_value) == 2) ){
		data_.float_t_to_uint_t( data );
	} else {
		data_.cast_from( data );
	}
}
	
	
	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::threshold_core( const size_t& nnz_core_, size_t& nnz_core_is_ )
{
	nnz_core_is_ = _core_comp.nnz();
	T_coeff threshold_value = 0.00001;
	while( nnz_core_is_ > nnz_core_ ) {
		_core_comp.threshold( threshold_value );
		nnz_core_is_ = _core_comp.nnz();
		
		//threshold value scheme
		if( threshold_value < 0.01) {
			threshold_value *= 10;
		} else if ( threshold_value < 0.2) {
			threshold_value += 0.05;
		} else if ( threshold_value < 1) {
			threshold_value += 0.25;
		} else if (threshold_value < 10 ) {
			threshold_value += 1;
		} else if (threshold_value < 50 ) {
			threshold_value += 10;
		} else if (threshold_value < 200 ) {
			threshold_value += 50;
		} else if (threshold_value < 500 ) {
			threshold_value += 100;
		} else if (threshold_value < 2000 ) {
			threshold_value += 500;
		} else if (threshold_value < 5000 ) {
			threshold_value += 3000;
		} else if (threshold_value >= 5000 ){
			threshold_value += 5000;
		}
	}
	_core.cast_from( _core_comp);
}



VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::threshold_core( const T_coeff& threshold_value_, size_t& nnz_core_ )
{
	_core_comp.threshold( threshold_value_ );
	nnz_core_ = _core_comp.nnz();
	_core.cast_from( _core_comp);
}
	
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::decompose( const t3_type& data_ ) 

{
	tucker_als( data_ );
}
	
VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::decompose( const t3_type& data_, 
								   T_internal& u1_min_, T_internal& u1_max_,
								   T_internal& u2_min_, T_internal& u2_max_,
								   T_internal& u3_min_, T_internal& u3_max_,
								   T_internal& core_min_, T_internal& core_max_ ) 
	
{
    decompose( data_ );
	
	quantize_basis_matrices( u1_min_, u1_max_, u2_min_, u2_max_, u3_min_, u3_max_ );
	quantize_core(core_min_, core_max_ );			
}

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::decompose( const t3_type& data_, 
								   T_internal& u_min_, T_internal& u_max_,
								   T_internal& core_min_, T_internal& core_max_ ) 

{
	decompose( data_ );
	
	quantize_basis_matrices( u_min_, u_max_ );
	quantize_core(core_min_, core_max_ );		
}
	

VMML_TEMPLATE_STRING
void 
VMML_TEMPLATE_CLASSNAME::tucker_als( const t3_type& data_ )
{
	t3_comp_type data;
	data.cast_from( data_ );

	t3_hooi< R1, R2, R3, I1, I2, I3, T_internal >::als( data, *_u1_comp, *_u2_comp, *_u3_comp, _core_comp ); 

	cast_members();
}



VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::reduce_ranks( const tucker3_tensor< K1, K2, K3, I1, I2, I3, T_value, T_coeff >& other )
//TuckerJI.rank_recuction(TuckerKI) K1 -> R1, K2 -> R2, K3 -> R3; I1, I2, I3 stay the same
{
     assert(R1 <= K1);
     assert(R2 <= K2);
     assert(R3 <= K3);	
             
     //reduce basis matrices
     matrix< I1, K1, T_coeff >* u1 = new matrix< I1, K1, T_coeff >();
     other.get_u1( *u1);
     for( size_t r1 = 0; r1 < R1; ++r1 ) 
     {
             _u1->set_column( r1, u1->get_column( r1 ));
     }
     
     matrix< I2, K2, T_coeff >* u2 = new matrix< I2, K2, T_coeff >();
     other.get_u2( *u2 );
     for( size_t r2 = 0; r2 < R2; ++r2) 
     {
             _u2->set_column( r2, u2->get_column( r2 ));
     }
     
     matrix< I3, K3, T_coeff >* u3 = new matrix< I3, K3, T_coeff >();
     other.get_u3( *u3 );
     for( size_t r3 = 0; r3 < R3; ++r3) 
     {
             _u3->set_column( r3, u3->get_column( r3 ));
     }
     
     //reduce core
	tensor3<K1, K2, K3, T_coeff > other_core;
     other.get_core( other_core );

     for( size_t r3 = 0; r3 < R3; ++r3 ) 
     {
          for( size_t r1 = 0; r1 < R1; ++r1 ) 
          {
              for( size_t r2 = 0; r2 < R2; ++r2 ) 
              {
                      _core.at( r1, r2, r3 ) = other_core.at( r1, r2, r3 );
              }
          }
     }

	
	cast_comp_members();

	delete u1;
	delete u2;
	delete u3;
}



VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::subsampling( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor )
{
     assert(I1 <= K1);
     assert(I1 <= K2);
     assert(I1 <= K3);	
     
     //subsample basis matrices
     matrix< K1, R1, T_coeff >* u1 = new matrix< K1, R1, T_coeff >();
     other.get_u1( *u1 );
     for( size_t i1 = 0, i = 0; i1 < K1; i1 += factor, ++i ) 
     {
             _u1->set_row( i, u1->get_row( i1 ));
     }
     
     matrix< K2, R2, T_coeff >* u2 = new matrix< K2, R2, T_coeff >();
     other.get_u2( *u2 );
     for( size_t i2 = 0,  i = 0; i2 < K2; i2 += factor, ++i) 
     {
             _u2->set_row( i, u2->get_row( i2 ));
     }
     
     matrix< K3, R3, T_coeff >* u3 = new matrix< K3, R3, T_coeff >() ;
     other.get_u3( *u3 );
     for( size_t i3 = 0,  i = 0; i3 < K3; i3 += factor, ++i) 
     {
             _u3->set_row( i, u3->get_row( i3 ));
     }
     
     other.get_core( _core );
	
	cast_comp_members();
	delete u1;
	delete u2;
	delete u3;
}


VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::subsampling_on_average( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, const size_t& factor )
{
    assert(I1 <= K1);
    assert(I1 <= K2);
    assert(I1 <= K3);	
    
    
    //subsample basis matrices
    matrix< K1, R1, T_coeff >* u1 = new matrix< K1, R1, T_coeff >();
    other.get_u1( *u1 );
    for( size_t i1 = 0, i = 0; i1 < K1; i1 += factor, ++i )
    {
            vector< R1, T_internal > tmp_row = u1->get_row( i1 );
            T_internal num_items_averaged = 1;
            for( size_t j = i1+1; (j < (factor+i1)) & (j < K1); ++j, ++num_items_averaged )
                    tmp_row += u1->get_row( j );

            tmp_row /= num_items_averaged;
            _u1->set_row( i, tmp_row);
    }
    
    matrix< K2, R2, T_coeff >* u2 = new matrix< K2, R2, T_coeff >();
    other.get_u2( *u2 );
    for( size_t i2 = 0,  i = 0; i2 < K2; i2 += factor, ++i) 
    {
            vector< R2, T_internal > tmp_row = u2->get_row( i2 );
            T_internal num_items_averaged = 1;
            for( size_t j = i2+1; (j < (factor+i2)) & (j < K2); ++j, ++num_items_averaged )
                    tmp_row += u2->get_row( j );

            tmp_row /= num_items_averaged;
            _u2->set_row( i, u2->get_row( i2 ));
    }
    
    matrix< K3, R3, T_coeff >* u3  = new matrix< K3, R3, T_coeff >();
    other.get_u3( *u3 );
    for( size_t i3 = 0,  i = 0; i3 < K3; i3 += factor, ++i) 
    {
            vector< R3, T_internal > tmp_row = u3->get_row( i3 );
            T_internal num_items_averaged = 1;
            for( size_t j = i3+1; (j < (factor+i3)) & (j < K3); ++j, ++num_items_averaged )
                    tmp_row += u3->get_row( j );
            
            tmp_row /= num_items_averaged;
            _u3->set_row( i, u3->get_row( i3 ));
    }
    
	other.get_core( _core );
	cast_comp_members();
	delete u1;
	delete u2;
	delete u3;
}




VMML_TEMPLATE_STRING
template< size_t K1, size_t K2, size_t K3>
void 
VMML_TEMPLATE_CLASSNAME::region_of_interest( const tucker3_tensor< R1, R2, R3, K1, K2, K3, T_value, T_coeff >& other, 
                                             const size_t& start_index1, const size_t& end_index1, 
                                             const size_t& start_index2, const size_t& end_index2, 
                                             const size_t& start_index3, const size_t& end_index3)
{
    assert(I1 <= K1);
    assert(I1 <= K2);
    assert(I1 <= K3);
    assert(start_index1 < end_index1);
    assert(start_index2 < end_index2);
    assert(start_index3 < end_index3);
    assert(end_index1 < K1);
    assert(end_index2 < K2);
    assert(end_index3 < K3);
    
    //region_of_interes of basis matrices
    matrix< K1, R1, T_internal >* u1 = new matrix< K1, R1, T_internal >();
    other.get_u1_comp( *u1 );
    for( size_t i1 = start_index1,  i = 0; i1 < end_index1; ++i1, ++i ) 
    {
            _u1_comp->set_row( i, u1->get_row( i1 ));
    }
    
    matrix< K2, R2, T_internal>* u2 = new matrix< K2, R2, T_internal>();
    other.get_u2_comp( *u2 );
    for( size_t i2 = start_index2,  i = 0; i2 < end_index2; ++i2, ++i) 
    {
            _u2_comp->set_row( i, u2->get_row( i2 ));
    }
	
    matrix< K3, R3, T_internal >* u3  = new matrix< K3, R3, T_internal>();
    other.get_u3_comp( *u3 );
    for( size_t i3 = start_index3,  i = 0; i3 < end_index3; ++i3, ++i) 
    {
            _u3_comp->set_row( i, u3->get_row( i3 ));
    }
    
    other.get_core_comp( _core_comp );
	
	//cast_comp_members();
	delete u1;
	delete u2;
	delete u3;
}
	
	
VMML_TEMPLATE_STRING
template< typename T >
void
VMML_TEMPLATE_CLASSNAME::export_to( std::vector< T >& data_ )
{
	
	data_.clear();
    
	cast_members();
	u1_const_iterator  it = _u1->begin(),
    it_end = _u1->end();
    for( ; it != it_end; ++it )
    {
        data_.push_back( static_cast< T >( *it) );
    }
    
    u2_const_iterator  u2_it = _u2->begin(),
    u2_it_end = _u2->end();
    for( ; u2_it != u2_it_end; ++u2_it )
    {
        data_.push_back(static_cast< T >(*u2_it) );
    }

    u3_const_iterator  u3_it = _u3->begin(),
    u3_it_end = _u3->end();
    for( ; u3_it != u3_it_end; ++u3_it )
    {
        data_.push_back(static_cast< T >( *u3_it) );
    }
    
    t3_core_iterator  it_core = _core.begin(),
    it_core_end = _core.end();
    for( ; it_core != it_core_end; ++it_core )
    {
        data_.push_back(static_cast< T >( *it_core) );
    }
}
	
	
VMML_TEMPLATE_STRING
template< typename T >
void
VMML_TEMPLATE_CLASSNAME::import_from( const std::vector< T >& data_ )
{
    size_t i = 0; //iterator over data_
    size_t data_size = (size_t) data_.size();

    if ( data_size != SIZE  )
        VMMLIB_ERROR( "import_from: the input data must have the size R1xR2xR3 + R1xI1 + R2xI2 + R3xI3 ", VMMLIB_HERE );
	
    u1_iterator  it = _u1->begin(),
    it_end = _u1->end();
    for( ; it != it_end; ++it, ++i )
    {
            *it = static_cast< T >( data_.at(i));
    }
    
    u2_iterator  u2_it = _u2->begin(),
    u2_it_end = _u2->end();
    for( ; u2_it != u2_it_end; ++u2_it, ++i )
    {
            *u2_it = static_cast< T >( data_.at(i));
    }
    
    u3_iterator  u3_it = _u3->begin(),
    u3_it_end = _u3->end();
    for( ; u3_it != u3_it_end; ++u3_it, ++i )
    {
            *u3_it = static_cast< T >( data_.at(i));
    }
    
    t3_core_iterator  it_core = _core.begin(),
    it_core_end = _core.end();
    for( ; it_core != it_core_end; ++it_core, ++i )
    {
            *it_core = static_cast< T >( data_.at(i));
    }
	
	cast_comp_members();
}
	
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::export_quantized_to( std::vector<unsigned char>& data_out_ )
{
	enable_quantify_coeff();
	//quantize tucker3 components (u1-u3 and core)
	size_t len_export_data = SIZE * sizeof(T_coeff) + 8*sizeof(T_internal);
	char * data = new char[ len_export_data ];
	size_t end_data = 0;
	size_t len_t_comp = sizeof( T_internal );
	
	//quantize basis matrices and copy min-max values
#if CODE_ALL_U_MIN_MAX	
	T_internal u1_min, u1_max, u2_min, u2_max, u3_min, u3_max;
	quantize_basis_matrices( u1_min, u1_max, u2_min, u2_max, u3_min, u3_max );
	memcpy( data, &u1_min, len_t_comp ); end_data = len_t_comp;
	memcpy( data + end_data, &u1_max, len_t_comp ); end_data += len_t_comp;
	memcpy( data + end_data, &u2_min, len_t_comp ); end_data += len_t_comp;
	memcpy( data + end_data, &u2_max, len_t_comp ); end_data += len_t_comp;
	memcpy( data + end_data, &u3_min, len_t_comp ); end_data += len_t_comp;
	memcpy( data + end_data, &u3_max, len_t_comp ); end_data += len_t_comp;
#else
	T_internal u_min, u_max;
	quantize_basis_matrices( u_min, u_max);
	memcpy( data, &u_min, len_t_comp ); end_data = len_t_comp;
	memcpy( data + end_data, &u_max, len_t_comp ); end_data += len_t_comp;
#endif
	
	//quantize core and copy min-max values
	T_internal core_min, core_max;
	quantize_core( core_min, core_max );		
	memcpy( data + end_data, &core_min, len_t_comp ); end_data += len_t_comp;
	memcpy( data + end_data, &core_max, len_t_comp ); end_data += len_t_comp;
	
	//copy data for u1
	size_t len_u1 = I1 * R1 * sizeof( T_coeff );
	memcpy( data + end_data, _u1, len_u1 ); end_data += len_u1;
	
	//copy data for u2
	size_t len_u2 = I2 * R2 * sizeof( T_coeff );
	memcpy( data + end_data, _u2, len_u2 ); end_data += len_u2;
	
	//copy data for u3
	size_t len_u3 = I3 * R3 * sizeof( T_coeff );
	memcpy( data + end_data, _u3, len_u3 ); end_data += len_u3;

	//copy data for core
	size_t len_core_slice = R1 * R2 * sizeof( T_coeff );
	for (size_t r3 = 0; r3 < R3; ++r3 ) {
		memcpy( data + end_data, _core.get_frontal_slice_fwd( r3 ), len_core_slice );
		end_data += len_core_slice;
	}
	
	data_out_.clear();
	for( size_t byte = 0; byte < len_export_data; ++byte )
	{
		data_out_.push_back( data[byte] );
	}
	delete[] data;
	
}


VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_quantized_from( const std::vector<unsigned char>& data_in_  )
{
	enable_quantify_coeff();
	size_t end_data = 0;
	size_t len_t_comp = sizeof( T_internal );
	size_t len_export_data = SIZE * sizeof(T_coeff) + 8*sizeof(T_internal);
	unsigned char * data = new unsigned char[ len_export_data ];
	for( size_t byte = 0; byte < len_export_data; ++byte )
	{
		data[byte] = data_in_.at(byte);
	}
	
	//copy min and max values: u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max
#if CODE_ALL_U_MIN_MAX	
	T_internal u1_min = 0; T_internal u1_max = 0;
	T_internal u2_min = 0; T_internal u2_max = 0;
	T_internal u3_min = 0; T_internal u3_max = 0;
	memcpy( &u1_min, data, len_t_comp ); end_data = len_t_comp;
	memcpy( &u1_max, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &u2_min, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &u2_max, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &u3_min, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &u3_max, data + end_data, len_t_comp ); end_data += len_t_comp;
#else
	T_internal u_min = 0; T_internal u_max = 0;
	memcpy( &u_min, data, len_t_comp ); end_data = len_t_comp;
	memcpy( &u_max, data + end_data, len_t_comp ); end_data += len_t_comp;
#endif
	
	T_internal core_min = 0; T_internal core_max = 0;
	memcpy( &core_min, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &core_max, data + end_data, len_t_comp ); end_data += len_t_comp;
		
	//copy data to u1
	size_t len_u1 = I1 * R1 * sizeof( T_coeff );
	memcpy( _u1, data + end_data, len_u1 ); end_data += len_u1;
	
	//copy data to u2
	size_t len_u2 = I2 * R2 * sizeof( T_coeff );
	memcpy( _u2, data + end_data, len_u2 ); end_data += len_u2;
	
	//copy data to u3
	size_t len_u3 = I3 * R3 * sizeof( T_coeff );
	memcpy( _u3, data + end_data, len_u3 ); end_data += len_u3;
	
	//copy data to core
	size_t len_core_slice = R1 * R2 * sizeof( T_coeff );
	front_core_slice_type* slice = new front_core_slice_type();
	for (size_t r3 = 0; r3 < R3; ++r3 ) {
		memcpy( slice, data + end_data, len_core_slice );
		_core.set_frontal_slice_fwd( r3, *slice );
		end_data += len_core_slice;
	}
	delete slice;
	delete[] data;
	
	//dequantize tucker3 components (u1-u3 and core)
#if CODE_ALL_U_MIN_MAX	
	dequantize_basis_matrices( u1_min, u1_max, u2_min, u2_max, u3_min, u3_max );
#else
	dequantize_basis_matrices( u_min, u_max, u_min, u_max, u_min, u_max  );
#endif
	
	dequantize_core( core_min, core_max );	
	
#if 0
        std::cout << "dequantized: " << std::endl << "u1-u3: " << std::endl
        << *_u1 << std::endl << *_u1_comp << std::endl
        << *_u2 << std::endl << *_u2_comp << std::endl
        << *_u3 << std::endl << *_u3_comp << std::endl
        << " core " << std::endl
        << _core << std::endl
        << " core_comp " << std::endl
        << _core_comp << std::endl;
#endif
}

VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::export_hot_quantized_to( std::vector<unsigned char>& data_out_ )
{
	enable_quantify_hot();
	//quantize tucker3 components (u1-u3 and core)
	size_t len_export_data = R1*R2*R3 + (R1*I1 + R2*I2 + R3*I3) * sizeof(T_coeff) + 4*sizeof(T_internal);
	char * data = new char[ len_export_data ];
	size_t end_data = 0;
	size_t len_t_comp = sizeof( T_internal );
	
	//quantize basis matrices and copy min-max values
	T_internal u_min, u_max;
	quantize_basis_matrices( u_min, u_max);
	memcpy( data, &u_min, len_t_comp ); end_data = len_t_comp;
	memcpy( data + end_data, &u_max, len_t_comp ); end_data += len_t_comp;
	
	//quantize core and copy min-max values
	T_internal core_min, core_max;
	quantize_core( core_min, core_max );		
	//memcpy( data + end_data, &core_min, len_t_comp ); end_data += len_t_comp; min_value is always zero in log quant
	memcpy( data + end_data, &core_max, len_t_comp ); end_data += len_t_comp;
	
	//copy first value of core tensor separately as a float
	memcpy( data + end_data, &_hottest_core_value, len_t_comp ); end_data += len_t_comp;
	
	//copy data for u1
	size_t len_u1 = I1 * R1 * sizeof( T_coeff );
	memcpy( data + end_data, _u1, len_u1 ); end_data += len_u1;
	
	//copy data for u2
	size_t len_u2 = I2 * R2 * sizeof( T_coeff );
	memcpy( data + end_data, _u2, len_u2 ); end_data += len_u2;
	
	//copy data for u3
	size_t len_u3 = I3 * R3 * sizeof( T_coeff );
	memcpy( data + end_data, _u3, len_u3 ); end_data += len_u3;
	
	//copy data for core
	size_t len_core_el = 1; //currently 1 bit for sign and 7 bit for values
	
	//colume-first iteration
	unsigned char core_el;
	for (size_t r3 = 0; r3 < R3; ++r3 ) {
		for (size_t r2 = 0; r2 < R2; ++r2 ) {
			for (size_t r1 = 0; r1 < R1; ++r1 ) {
				core_el = (_core.at( r1, r2, r3 ) | (_signs.at( r1, r2, r3) * 0x80 ));
				/*std::cout << "value: " << int(_core.at( r1, r2, r3 )) << " bit " << int( core_el ) 
				<< " sign: " << int(_signs.at( r1, r2, r3)) << std::endl;*/
				memcpy( data + end_data, &core_el, len_core_el );
				++end_data;
			}
		}
	} 
		
	data_out_.clear();
	for( size_t byte = 0; byte < len_export_data; ++byte )
	{
		data_out_.push_back( data[byte] );
	}
	delete[] data;
}
	
	
	
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_hot_quantized_from( const std::vector<unsigned char>& data_in_  )
{
	enable_quantify_hot();
	size_t end_data = 0;
	size_t len_t_comp = sizeof( T_internal );
	size_t len_export_data = R1*R2*R3 + (R1*I1 + R2*I2 + R3*I3) * sizeof(T_coeff) + 4*sizeof(T_internal);
	unsigned char * data = new unsigned char[ len_export_data ];
	for( size_t byte = 0; byte < len_export_data; ++byte )
	{
		data[byte] = data_in_.at(byte);
	}
	
	//copy min and max values: u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max
	T_internal u_min = 0; T_internal u_max = 0;
	memcpy( &u_min, data, len_t_comp ); end_data = len_t_comp;
	memcpy( &u_max, data + end_data, len_t_comp ); end_data += len_t_comp;
	
	T_internal core_min = 0; T_internal core_max = 0; //core_min is 0
	//memcpy( &core_min, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &core_max, data + end_data, len_t_comp ); end_data += len_t_comp;
	//copy first value of core tensor separately as a float
	memcpy( &_hottest_core_value, data + end_data, len_t_comp ); end_data += len_t_comp;
	
	//copy data to u1
	size_t len_u1 = I1 * R1 * sizeof( T_coeff );
	memcpy( _u1, data + end_data, len_u1 ); end_data += len_u1;
	
	//copy data to u2
	size_t len_u2 = I2 * R2 * sizeof( T_coeff );
	memcpy( _u2, data + end_data, len_u2 ); end_data += len_u2;
	
	//copy data to u3
	size_t len_u3 = I3 * R3 * sizeof( T_coeff );
	memcpy( _u3, data + end_data, len_u3 ); end_data += len_u3;
	
	//copy data to core
	size_t len_core_el = 1; //currently 1 bit for sign and 7 bit for values

	unsigned char core_el;
	for (size_t r3 = 0; r3 < R3; ++r3 ) {
		for (size_t r2 = 0; r2 < R2; ++r2 ) {
			for (size_t r1 = 0; r1 < R1; ++r1 ) {
				memcpy( &core_el, data + end_data, len_core_el );
				_signs.at( r1, r2, r3 ) = (core_el & 0x80)/128;
				_core.at( r1, r2, r3 ) = core_el & 0x7f ;
				++end_data;
			}
		}
	} 
	//std::cout << "signs: " << _signs << std::endl;
	//std::cout << "_core: " << _core << std::endl;
	
	delete[] data;
	
	//dequantize tucker3 components (u1-u3 and core)
	dequantize_basis_matrices( u_min, u_max, u_min, u_max, u_min, u_max  );
	
	dequantize_core( core_min, core_max );	
	
#if 0
	std::cout << "dequantized: " << std::endl << "u1-u3: " << std::endl
	<< *_u1 << std::endl << *_u1_comp << std::endl
	<< *_u2 << std::endl << *_u2_comp << std::endl
	<< *_u3 << std::endl << *_u3_comp << std::endl
	<< " core " << std::endl
	<< _core << std::endl
	<< " core_comp " << std::endl
	<< _core_comp << std::endl;
#endif
}
	
VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::export_ttm_quantized_to( std::vector<unsigned char>& data_out_ )
{
	enable_quantify_log();
	//quantize tucker3 components (u1-u3 and core)
	size_t len_export_data = R1*R2*R3 + (R1*I1 + R2*I2 + R3*I3) * sizeof(T_coeff) + 3*sizeof(T_internal);
	char * data = new char[ len_export_data ];
	size_t end_data = 0;
	size_t len_t_comp = sizeof( T_internal );
	
	//quantize basis matrices and copy min-max values
	T_internal u_min, u_max;
	quantize_basis_matrices( u_min, u_max);
	memcpy( data, &u_min, len_t_comp ); end_data = len_t_comp;
	memcpy( data + end_data, &u_max, len_t_comp ); end_data += len_t_comp;
	
	//quantize core and copy min-max values
	T_internal core_min, core_max;
	quantize_core( core_min, core_max );		
	//memcpy( data + end_data, &core_min, len_t_comp ); end_data += len_t_comp; min_value is always zero in log quant
	memcpy( data + end_data, &core_max, len_t_comp ); end_data += len_t_comp;
	
	//copy data for u1
	size_t len_u1 = I1 * R1 * sizeof( T_coeff );
	memcpy( data + end_data, _u1, len_u1 ); end_data += len_u1;
	
	//copy data for u2
	size_t len_u2 = I2 * R2 * sizeof( T_coeff );
	memcpy( data + end_data, _u2, len_u2 ); end_data += len_u2;
	
	//copy data for u3
	size_t len_u3 = I3 * R3 * sizeof( T_coeff );
	memcpy( data + end_data, _u3, len_u3 ); end_data += len_u3;
	
	//copy data for core
	size_t len_core_el = 1; //currently 1 bit for sign and 7 bit for values
	
	//colume-first iteration
	//backward cylcling after lathauwer et al. 
	unsigned char core_el;
	for (size_t r2 = 0; r2 < R2; ++r2 ) {
		for (size_t r3 = 0; r3 < R3; ++r3 ) {
			for (size_t r1 = 0; r1 < R1; ++r1 ) {
				core_el = (_core.at( r1, r2, r3 ) | (_signs.at( r1, r2, r3) * 0x80 ));
				/*std::cout << "value: " << int(_core.at( r1, r2, r3 )) << " bit " << int( core_el ) 
				 << " sign: " << int(_signs.at( r1, r2, r3)) << std::endl;*/
				memcpy( data + end_data, &core_el, len_core_el );
				++end_data;
			}
		}
	} 
	
	data_out_.clear();
	for( size_t byte = 0; byte < len_export_data; ++byte )
	{
		data_out_.push_back( data[byte] );
	}
	delete[] data;
}

VMML_TEMPLATE_STRING
void
VMML_TEMPLATE_CLASSNAME::import_ttm_quantized_from( const std::vector<unsigned char>& data_in_  )
{
	enable_quantify_log();
	size_t end_data = 0;
	size_t len_t_comp = sizeof( T_internal );
	size_t len_export_data = R1*R2*R3 + (R1*I1 + R2*I2 + R3*I3) * sizeof(T_coeff) + 3*sizeof(T_internal);
	unsigned char * data = new unsigned char[ len_export_data ];
	for( size_t byte = 0; byte < len_export_data; ++byte )
	{
		data[byte] = data_in_.at(byte);
	}
	
	//copy min and max values: u1_min, u1_max, u2_min, u2_max, u3_min, u3_max, core_min, core_max
	T_internal u_min = 0; T_internal u_max = 0;
	memcpy( &u_min, data, len_t_comp ); end_data = len_t_comp;
	memcpy( &u_max, data + end_data, len_t_comp ); end_data += len_t_comp;
	
	T_internal core_min = 0; T_internal core_max = 0; //core_min is 0
	//memcpy( &core_min, data + end_data, len_t_comp ); end_data += len_t_comp;
	memcpy( &core_max, data + end_data, len_t_comp ); end_data += len_t_comp;
	
	//copy data to u1
	size_t len_u1 = I1 * R1 * sizeof( T_coeff );
	memcpy( _u1, data + end_data, len_u1 ); end_data += len_u1;
	
	//copy data to u2
	size_t len_u2 = I2 * R2 * sizeof( T_coeff );
	memcpy( _u2, data + end_data, len_u2 ); end_data += len_u2;
	
	//copy data to u3
	size_t len_u3 = I3 * R3 * sizeof( T_coeff );
	memcpy( _u3, data + end_data, len_u3 ); end_data += len_u3;
	
	//copy data to core
	size_t len_core_el = 1; //currently 1 bit for sign and 7 bit for values
	
	//backward cyclic after lathauwer et al. 
	unsigned char core_el;
	for (size_t r2 = 0; r2 < R2; ++r2 ) {
		for (size_t r3 = 0; r3 < R3; ++r3 ) {
			for (size_t r1 = 0; r1 < R1; ++r1 ) {
				memcpy( &core_el, data + end_data, len_core_el );
				_signs.at( r1, r2, r3 ) = (core_el & 0x80)/128;
				_core.at( r1, r2, r3 ) = core_el & 0x7f ;
				++end_data;
			}
		}
	} 
	//std::cout << "signs: " << _signs << std::endl;
	//std::cout << "_core: " << _core << std::endl;
	
	delete[] data;
	
	//dequantize tucker3 components (u1-u3 and core)
	dequantize_basis_matrices( u_min, u_max, u_min, u_max, u_min, u_max  );
	dequantize_core( core_min, core_max );	
}
	
	
	
VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::nnz() const
{
	size_t counter = 0;
	
	counter += _u1_comp->nnz();
	counter += _u2_comp->nnz();
	counter += _u3_comp->nnz();
	counter += _core_comp.nnz();
	
	return counter;
}
	
VMML_TEMPLATE_STRING
size_t
VMML_TEMPLATE_CLASSNAME::nnz( const T_value& threshold ) const
{
	size_t counter = 0;
	
	counter += _u1_comp->nnz( threshold );
	counter += _u2_comp->nnz( threshold );
	counter += _u3_comp->nnz( threshold );
	counter += _core_comp.nnz( threshold );

	return counter;
}

	
#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME
	
} // namespace vmml

#endif
