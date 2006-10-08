#ifndef _vmml_OctreeNode_H_
#define _vmml_OctreeNode_H_

#include <iostream>
#include <sstream>
#include <vector>
#include "VMMLib.h"
//#include "Pool.h"
#include "AxisAlignedBoundingBox.h"

/** 
* an octree class 
* 
* @author jonas boesch
**/

namespace vmml
{

//
// declaration
//

template< typename T >
class OctreeNode
{
    static const unsigned char ROOT = 255;
    static const unsigned char XBIT = 4;
    static const unsigned char YBIT = 2;
    static const unsigned char ZBIT = 1;

//    typedef Pool< OctreeNode< T > > ONodePool;

public:
    // root node ctor
    // creates only the root node
    OctreeNode( const Aabbf& aabb );

    // child node ctor
    OctreeNode( OctreeNode< T >* parent, unsigned char pIndex );
    
    // dtor
	~OctreeNode();
    
    // faster version
    // PRECONDITION: pos must be in the node's aabb
    // (this is not checked).
    OctreeNode< T >* getLeaf( const Vector3f& pos );

    // safe version of getLeaf
    // returns 0 if pos is outside of the root aabb
    OctreeNode< T >* findLeaf( const Vector3f& pos );
   
    //OctreeNode< T >* findNode( const Sphere< float >& sphere );
     
    // getters / setters
    inline void setLoadPtr( T* load ) { _load = load; };
    inline T* getLoadPtr() { return _load; };

    inline const Aabbf& getAabb() { return _aabb; };
    inline const Vector3f& getCenter() { return _center; };

    inline unsigned char getIndex() { return _index; };
    inline unsigned char getDepth() { return _depth; };
    inline bool isRoot() { return ! _parent; };
    inline bool isLeaf() { return ! _children.size(); };
    inline bool hasChildren() { return _children.size(); };

    inline OctreeNode< T >* getChild( size_t index );
    inline OctreeNode< T >* operator[]( size_t index );
    
    void createChildren( size_t maxDepth );
    
    // warning: slow!    
    void spam();

    // warning: slow! 
    // returns the path-index as string...
    void getId( std::string& nodeId );

protected:
    // dangerous ctor/setPackage package for pool
    OctreeNode(); 
    void _setParent( OctreeNode< T >* parent, unsigned char index );

    OctreeNode< T >* _parent;
    std::vector< OctreeNode< T >* > _children;
    Vector3f _center;
    Aabbf _aabb;
    T* _load;
    unsigned char _index;
    unsigned char _depth;

    //ONodePool _pool;
};

//
// implementation
//

template< typename T >
OctreeNode< T >::OctreeNode( const Aabbf& aabb )
    : _parent( 0 )
    , _children( 0 )
    , _aabb( aabb )
    , _load( 0 )
    , _index( ROOT )
    , _depth( 0 )
{
    _center = _aabb.getMin();
    float edgelen = _aabb.getMax().x - _center.x;
    edgelen *= 0.5f;
    _center += edgelen;
}

template< typename T >
OctreeNode< T >::OctreeNode()
    : _parent( 0 )
    , _children( 0 )
    , _aabb()
    , _load( 0 )
    , _index( ROOT )
    , _depth( 0)
{}

template< typename T >
OctreeNode< T >::OctreeNode( OctreeNode< T >* parent, unsigned char pIndex )
    : _parent( parent )
    , _children( 0 )
    , _aabb()
    , _load( 0 )
    , _index( pIndex )
    , _depth( parent->getDepth() + 1 )

{
    assert( _parent );
    const Aabbf& pbox = parent->getAabb();
    Vector3f min = pbox.getMin();
    Vector3f max = pbox.getMax();
    // in this context, the aabb is always a cube.
    float edgelen = max.x - min.x;
    edgelen *= 0.5f;
    // compute the child aabb
    if ( _index & XBIT )
        min.x += edgelen;
    else 
        max.x -= edgelen;
        
    if ( _index & YBIT )
        min.y += edgelen;
    else
        max.y -= edgelen;
    
    if ( _index & ZBIT )
        min.z += edgelen;
    else
        max.z -= edgelen;
    _aabb.set( min, max );

    //std::cout << "OctreeNode< T > _depth " << (int) _depth << " _index " 
    //    << (int) _index << " min " << min << " max " << max << std::endl;
    // computing the center (split) point
    edgelen *= 0.5f;
    _center = min + edgelen;
    //std::cout << "OctreeNode< T > " << (int)_depth << "." << (int) _index << " center " << min << std::endl;
}

template< typename T >
void OctreeNode< T >::_setParent( OctreeNode< T >* parent, unsigned char index )
{
    assert( _parent );
    _parent = parent;

    _index =  index;
    _depth = parent->getDepth() + 1;
    
    const Aabbf& pbox = parent->getAabb();
    Vector3f min = pbox.getMin();
    Vector3f max = pbox.getMax();
    // in this context, the aabb is always a cube.
    float edgelen = max.x - min.x;
    edgelen *= 0.5f;
    // compute the child aabb
    if ( _index & XBIT )
        min.x += edgelen;
    else 
        max.x -= edgelen;
        
    if ( _index & YBIT )
        min.y += edgelen;
    else
        max.y -= edgelen;
    
    if ( _index & ZBIT )
        min.z += edgelen;
    else
        max.z -= edgelen;
    _aabb.set( min, max );

    //std::cout << "OctreeNode< T > _depth " << (int) _depth << " _index " 
    //    << (int) _index << " min " << min << " max " << max << std::endl;
    // computing the center (split) point
    edgelen *= 0.5f;
    _center = min + edgelen;
    //std::cout << "OctreeNode< T > " << (int)_depth << "." << (int) _index << " center " << min << std::endl;
}

template< typename T >
OctreeNode< T >::~OctreeNode()
{
    if ( _children.size() > 0 )
    {
        typename std::vector< OctreeNode< T >* >::iterator it = _children.begin();
        for ( ; it != _children.end(); ++it )
        {
            delete *it;
            //_pool.destroy( *it );
        }
    }
}

template< typename T >
OctreeNode< T >* OctreeNode< T >::getLeaf( const Vector3f& pos )
{
    // if this is a leaf node, return this;
    if ( _children.size() == 0 )
        return this;
    // else determine which child is requested
    unsigned char child = 0;
    if ( pos.x > _center.x )
        child += XBIT;
    if ( pos.y > _center.y ) 
        child += YBIT;
    if ( pos.z > _center.z ) 
        child += ZBIT;
    return _children[child]->getLeaf( pos );
}

template< typename T >
OctreeNode< T >* OctreeNode< T >::findLeaf( const Vector3f& pos )
{
    if ( ! _aabb.isIn( pos ) )
    { // pos is not in this node's aabb
        if ( _parent )
            _parent->findLeaf( pos );
        else // this is the root node and pos is outside
            return 0;
    }
    // if this is a leaf node, return this;
    if ( _children.size() == 0 )
        return this;
    // else determine which child is requested
    unsigned char child = 0;
    if ( pos.x > _center.x )
        child += XBIT;
    if ( pos.y > _center.y ) 
        child += YBIT;
    if ( pos.z > _center.z ) 
        child += ZBIT;
    return _children[child]->getLeaf( pos );
}

template< typename T>
void OctreeNode< T >::createChildren( size_t maxDepth )
{
    if ( _depth >= maxDepth )
        return;
    assert( _children.size() == 0 );
    for ( unsigned char i = 0; i < 8; ++i )
    {
        // FIXME pool
        OctreeNode< T >* child = /*_pool.create(); */ new OctreeNode< T >( this, i ); 
        child->_setParent( this, i );
        _children.push_back( child );
        child->createChildren( maxDepth );
    }
}

template< typename T>
void OctreeNode< T >::getId( std::string& nodeId )
{
    if ( _depth >= nodeId.size() )
    {
        nodeId = "";
        for ( ssize_t i = -1; i < _depth; ++i )
        {
            nodeId.push_back( 'x' );
        }
    }
    if ( _parent )
    {
        std::stringstream ss;
        ss << (unsigned int)_index;
        nodeId.replace( _depth, 1, ss.str() );
        _parent->getId( nodeId );    
    }
    else
       nodeId.replace( _depth, 1, "R" ); 
}


template< typename T>
void OctreeNode< T >::spam()
{
    std::string nodeId;
    getId( nodeId );
    std::cout << "OctreeNode " << nodeId << std::endl 
        << "  maximum: " << _aabb.getMax() << std::endl 
        << "  minimum: " << _aabb.getMin() << std::endl
        << "  center:  " << _center << std::endl;
}

template< typename T >
inline OctreeNode< T >* OctreeNode< T >::getChild( size_t index )
{ 
    assert( index < 8 ); 
    return _children[index]; 
};

template< typename T >    
inline OctreeNode< T >* OctreeNode< T >::operator[]( size_t index )
{ 
    assert( index < 8 ); 
    return _children[index]; 
};
    



};

#endif
