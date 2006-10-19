#include "VMMLib.h"

#include "AxisAlignedBoundingBox.h"
#include "QuadtreeNode.h"
#include "OctreeNode.h"

#include <vector>

typedef std::vector< vmml::Vector3f > PointVector;
typedef vmml::OctreeNode< PointVector > ONode;
typedef vmml::QuadtreeNode< PointVector > QuadNode;

using vmml::OctreeNode;
using vmml::QuadtreeNode;
using vmml::Vector3f;

int main( int argc, char** argv )
{
    QuadNode qroot( .5f, .5f, 1.f );
    qroot.createChildren( 5 );
    Vector3f p0( .9f, .9f, .3f );
    Vector3f p1( .1f, .3f, .7f );

    QuadNode* qn0 = qroot.getLeaf( p0 );
    QuadNode* qn1 = qroot.getLeaf( p1 );
    std::cout << "Point 0 : " << p0 << std::endl;
    qn0->spam();
    std::cout << "Point 1 : " << p0 << std::endl;
    qn1->spam();
    std::cout << std::endl;

    ONode oroot( .5f, .5f, .5f, 1.f );
    oroot.createChildren( 5 );

    ONode* on0 = oroot.getLeaf( p0 );
    ONode* on1 = oroot.getLeaf( p1 );
    std::cout << "Point 0 : " << p0 << std::endl;
    on0->spam();
    std::cout << "Point 1 : " << p0 << std::endl;
    on1->spam();
    std::cout << std::endl;
    

};