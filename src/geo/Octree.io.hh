#ifndef D6_OCTREE_IO_HH
#define D6_OCTREE_IO_HH

#include "Octree.hh"

#include <boost/serialization/split_member.hpp>

namespace d6 {

template < typename Archive >
void Octree::save( Archive &ar, unsigned int ) const
{
    ar << m_dim ;
    ar << m_dx ;
    ar << m_maxDepth ;
    ar << m_trees ;
}

template < typename Archive >
void Octree::load( Archive &ar, unsigned int )
{
    ar >> m_dim ;
    ar >> m_dx ;
    ar >> m_maxDepth ;
    ar >> m_trees ;

    rebuild() ;
}


template < typename Archive >
void Octree::serialize( Archive &ar, unsigned int file_version )
{
    boost::serialization::split_member(ar, *this, file_version);
}



}

#endif
