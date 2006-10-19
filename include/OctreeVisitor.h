#ifndef _vmml_OctreeVisitor_H_
#define _vmml_OctreeVisitor_H_

namespace vmml
{

template< typename T >
class OctreeNode;

template< typename T >
class OctreeVisitor
{
public:
    virtual void visit( OctreeNode< T >* node ) = 0;
};

}; // namespace vmml

#endif
