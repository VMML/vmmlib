/***************************************************************************
*   Copyright (C) 2004 by jonas boesch									  *
*   boesch@ifi.unizh.ch													*
*                                                                         *
***************************************************************************/

#include <vmmlib/vmmlib.h>

#include <iostream>
#include <vector>

#include "Matrix3Test.h"
#include "Matrix4Test.h"
#include "Vector2Test.h"
#include "Vector3Test.h"
#include "Vector4Test.h"
#include "QuaternionTest.h"
#include "SVDTest.h"
#include "JacobiSolverTest.h"

using std::cout;
using std::endl;

int main()
{	
    vmml::Matrix3Test m3t;
    bool ok = m3t.test();

    vmml::Matrix4Test m4t;
    if ( ok ) 
        ok = m4t.test();
		
	vmml::Vector2Test v2t;
	if ( ok ) 
        ok = v2t.test();
		
	vmml::Vector3Test v3t;
	if ( ok ) 
        ok = v3t.test();
	
	vmml::Vector4Test v4t;
	if ( ok ) 
        ok = v4t.test();
      
	vmml::QuaternionTest qtest;
	if ( ok )
		ok = qtest.test();
		
	vmml::SVDTest stest;
	if ( ok )
		ok = stest.test();
		
	vmml::JacobiSolverTest jtest;
	if ( ok )
		ok = jtest.test();

    if ( ok ) 
        cout << "All tests passed!" << endl;
    else
        cout << "Some tests have failed!" << endl;    
    return 0;
}
