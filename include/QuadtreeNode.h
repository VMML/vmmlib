#ifndef _vmml_QuadtreeNode_H_
#define _vmml_QuadtreeNode_H_

#include <iostream>
#include <sstream>
#include <vector>
#include "VMMLib.h"
//#include "Pool.h"
#include "AxisAlignedBoundingBox.h"
#include "QuadtreeVisitor.h"


/** 
* an balanced quadtree class
*
* note: uses Vector3<> instead of Vector2<> ( and just ignores z )
*
* @author jonas boesch
**/

namespace vmml
{

//
// declaration
//

template< typename T >
class QuadtreeNode
{
    static const unsigned char ROOT = 255;
    static const unsigned char XBIT = 1;
    static const unsigned char YBIT = 2;

//    typedef Pool< QuadtreeNode< T > > ONodePool;

public:

    // root node ctor
    // creates only the root node
    QuadtreeNode( const Aabbf& aabb );
    QuadtreeNode( float cx, float cy, float size );

    // child node ctor
    QuadtreeNode( QuadtreeNode< T >* parent, unsigned char pIndex );
    
    // dtor
	~QuadtreeNode();
    
    // faster version
    // PRECONDITION: pos must be in the node's aabb
    // (this is not checked).
    QuadtreeNode< T >* getLeaf( const Vector3f& pos );

    // safe version of getLeaf
    // returns 0 if pos is outside of the root aabb
    QuadtreeNode< T >* findLeaf( const Vector3f& pos );
   
    //QuadtreeNode< T >* findNode( const Sphere< float >& sphere );
     
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

    inline QuadtreeNode< T >* getChild( size_t index );
    inline QuadtreeNode< T >* operator[]( size_t index );
    
    void createChildren( size_t maxDepth );
    void setupLeaves();
    
    void accept( QuadtreeVisitor< T >* visitor );
    
    // warning: slow!    
    void spam();

    // warning: slow! 
    // returns the path-index as string...
    void getId( std::string& nodeId );

protected:
    QuadtreeNode(); 
    void _setupRoot( const Aabbf& aabb ); 
    void _setParent( QuadtreeNode< T >* parent, unsigned char index );

    QuadtreeNode< T >* _parent;
    std::vector< QuadtreeNode< T >* > _children;
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
QuadtreeNode< T >::QuadtreeNode( const Aabbf& aabb )
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
QuadtreeNode< T >::QuadtreeNode( float cx, float cy, float size )
    : _parent( 0 )
    , _children( 0 )
    , _aabb( cx, cy, 0, size )
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
QuadtreeNode< T >::QuadtreeNode( QuadtreeNode< T >* parent, unsigned char pIndex )
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

    _aabb.set( min, max );

    //std::cout << "QuadtreeNode< T > _depth " << (int) _depth << " _index " 
    //    << (int) _index << " min " << min << " max " << max << std::endl;
    // computing the center (split) point
    edgelen *= 0.5f;
    _center = min + edgelen;
    //std::cout << "QuadtreeNode< T > " << (int)_depth << "." << (int) _index << " center " << min << std::endl;
}

template< typename T >
QuadtreeNode< T >::QuadtreeNode()
    : _parent( 0 )
    , _children( 0 )
    , _aabb()
    , _load( 0 )
    , _index( ROOT )
    , _depth( 0)
{}

template< typename T >
void QuadtreeNode< T >::_setupRoot( const Aabbf& aabb )
{
    _aabb = aabb;
    _center = _aabb.getMin();
    float edgelen = _aabb.getMax().x - _center.x;
    edgelen *= 0.5f;
    _center += edgelen;
}

template< typename T >
void QuadtreeNode< T >::_setParent( QuadtreeNode< T >* parent, unsigned char index )
{
    assert( _parent );
    _parent = parent;

    _index =  index;
    _depth = parent->getDepth() + 1;
    
    _aabb = parent->getAabb();
    Vector3f min = _aabb.getMin();
    Vector3f max = _aabb.getMax();
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
    _aabb.set( min, max );

    //std::cout << "QuadtreeNode< T > _depth " << (int) _depth << " _index " 
    //    << (int) _index << " min " << min << " max " << max << std::endl;
    // computing the center (split) point
    edgelen *= 0.5f;
    _center = min + edgelen;
    //std::cout << "QuadtreeNode< T > " << (int)_depth << "." << (int) _index << " center " << min << std::endl;
}

template< typename T >
QuadtreeNode< T >::~QuadtreeNode()
{
    if ( _children.size() > 0 )
    {
        typename std::vector< QuadtreeNode< T >* >::iterator it = _children.begin();
        for ( ; it != _children.end(); ++it )
        {
            delete *it;
            //_pool.destroy( *it );
        }
    }
}

template< typename T >
QuadtreeNode< T >* QuadtreeNode< T >::getLeaf( const Vector3f& pos )
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
    return _children[child]->getLeaf( pos );
}

template< typename T >
QuadtreeNode< T >* QuadtreeNode< T >::findLeaf( const Vector3f& pos )
{
    if ( ! _aabb.isIn2d( pos ) )
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
    return _children[child]->getLeaf( pos );
}

template< typename T>
void QuadtreeNode< T >::createChildren( size_t maxDepth )
{
    if ( _depth >= maxDepth )
        return;
    assert( _children.size() == 0 );
    for ( unsigned char i = 0; i < 4; ++i )
    {
        // FIXME pool
        QuadtreeNode< T >* child = /*_pool.create(); */ new QuadtreeNode< T >( this, i ); 
        child->_setParent( this, i );
        _children.push_back( child );
        child->createChildren( maxDepth );
    }
}

template< typename T>
void QuadtreeNode< T >::getId( std::string& nodeId )
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
void QuadtreeNode< T >::spam()
{
    std::string nodeId;
    getId( nodeId );
    std::cout << "QuadtreeNode " << nodeId << std::endl 
        << "  maximum: " << _aabb.getMax() << std::endl 
        << "  minimum: " << _aabb.getMin() << std::endl
        << "  center:  " << _center << std::endl;
}

template< typename T >
inline QuadtreeNode< T >* QuadtreeNode< T >::getChild( size_t index )
{ 
    assert( index < 4 ); 
    return _children[index]; 
};

template< typename T >    
inline QuadtreeNode< T >* QuadtreeNode< T >::operator[]( size_t index )
{ 
    assert( index < 4 ); 
    return _children[index]; 
};
    
template< typename T >    
void QuadtreeNode< T >::setupLeaves()
{ 
    if ( _children.size() > 0 )
    {
        for ( size_t i = 0; i < 4; ++i )
        {
            _children[i]->setupLeaves();
        }
    }
    //else if ( _load == 0 )
    //{
    //    _load = new T();
    //}
};

template< typename T > 
void QuadtreeNode< T >::accept( QuadtreeVisitor< T >* visitor )
{
    visitor->visit( this );
}


};

#endif
