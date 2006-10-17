/***************************************************************************
*   Copyright (C) 2004 by jonas boesch									  *
*   boesch@ifi.unizh.ch													*
*                                                                         *
***************************************************************************/

#include "VMMLib.h"

#include "JacobiSolver.h"
#include "SingularValueDecomposition.h"
#include "AxisAlignedBoundingBox.h"
#include "OctreeNode.h"

#include <iostream>
#include <vector>

#include "Matrix3Test.h"
#include "Matrix4Test.h"

using namespace std;
int main()
{
    vmml::Matrix3Test m3t;
    bool ok = m3t.test();

    vmml::Matrix4Test m4t;
    if ( ok ) 
        ok = m4t.test();
        
    if ( ok ) 
        cout << "All tests passed!" << endl;
    else
        cout << "Some tests have failed!" << endl;    
    return 0;
}