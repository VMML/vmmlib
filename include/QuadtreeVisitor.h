#ifndef _vmml_QuadtreeVisitor_H_
#define _vmml_QuadtreeVisitor_H_

namespace vmml
{
template< typename T >
class QuadtreeNode;

template< typename T >
class QuadtreeVisitor
{
public:
    virtual void visit( QuadtreeNode< T >* node ) = 0;
};

}; // namespace vmml

#endif
