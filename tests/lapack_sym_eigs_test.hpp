#ifndef __VMML__LAPACK_SYM_EIGS_TEST__HPP__
#define __VMML__LAPACK_SYM_EIGS_TEST__HPP__

#include "unit_test.hpp"

namespace vmml
{
	
	class lapack_sym_eigs_test : public unit_test
	{
	public:
		lapack_sym_eigs_test() : unit_test( "symmetric eigenvalue decomposition using lapack" ) {}
		virtual bool run();
		
	protected:
		
	}; // class lapack_sym_eigs_test
	
} // namespace vmml

#endif