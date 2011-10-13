#ifndef __VMML__TUCK3_EXPORTER_IMPORTER_TEST__HPP__
#define __VMML__TUCK3_EXPORTER_IMPORTER_TEST__HPP__

#include "unit_test.hpp"
#include <string>

namespace vmml
{
	
	class tucker3_exporter_importer_test : public unit_test
	{
	public:
		tucker3_exporter_importer_test() : unit_test( "tucker3 tensor exporter/importer test" ) {}
		bool run();
		
	protected:
		
	}; // class tucker3_exporter_importer_test
	
} // namespace vmml

#endif


