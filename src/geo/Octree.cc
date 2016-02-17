#include "Octree.hh"

#include <iostream>

namespace d6 {

struct Octree::SubTree::SubCell {
    Index res ;
    ArrWi coords ;
};

Octree::SubTree::SubTree()
    : m_offset(0), m_nodes{{1,1}}
{}

Index Octree::SubTree::nLeafs() const
{
    return m_nodes[0].nLeafs ;
}

void Octree::SubTree::find( const Vec& pos, Index &leafIndex, Coords& coords ) const  {
    Arr rel = pos ;
    leafIndex = 0 ;
    find( 0, rel, leafIndex, coords ) ;
    leafIndex += m_offset ;
}

void Octree::SubTree::compute_geo( Index leafIndex, Voxel &geo ) const  {
    SubCell cell ;
    subcell( leafIndex, cell ) ;
    geo.origin = cell.coords.cast< Scalar >() / cell.res ;
    geo.box    = Vec::Constant(1./cell.res) ;
}

void Octree::SubTree::upres_cell( Index leafIndex, Index target_res, ArrWi& hires_cell ) const
{
    SubCell cell ;
    subcell( leafIndex, cell ) ;
    hires_cell = cell.coords * target_res / cell.res ;
}


void Octree::SubTree::find( Index offset, Arr& rel, Index &leafIndex, Coords& coords ) const
{

    if( m_nodes[ offset ].subSize == 1 ) {
        coords = rel ;
    } else {
        ArrWi sub = (2*rel).cast< Index >().min(1) ;
        Index subi = Voxel::cornerIndex( sub ) ;

        ++offset ;
        for (Index i = 0 ; i < subi ; ++i) {
            offset    += m_nodes[ offset ].subSize ;
            leafIndex += m_nodes[ offset ].nLeafs ;
        }

        rel = 2*rel - sub.cast<Scalar>() ;

        find( offset, rel, leafIndex, coords ) ;
    }
}

void Octree::SubTree::subcell ( Index leafIndex, SubCell& cell,
               Index offset = 0, Index leafOffset = 0 ) const
{

    const Node& cur = m_nodes[offset] ;
    assert( leafIndex < leafOffset + cur.nLeafs ) ;
    if( cur.subSize == 1 ) {
        cell.res = 1 ;
        cell.coords.setZero() ;
    } else {
        ++offset ;

        Index subi = 0 ;
        while( leafOffset + m_nodes[ offset ].nLeafs <= leafIndex ) {
            ++ subi ;
            offset     += m_nodes[ offset ].subSize ;
            leafOffset += m_nodes[ offset ].nLeafs ;
        }

        const ArrWi sub = Voxel::corner( subi ) ;
        subcell( leafIndex, cell, offset, leafOffset ) ;

        cell.coords += cell.res * sub ;
        cell.res *= 2 ;
    }
}

void Octree::SubTree::split( Index leafIndex )
{
    leafIndex -= m_offset ;
    Index offset = 0 ;
    split( leafIndex, offset, 0 ) ;

    const std::size_t orig_size = m_nodes.size() ;
    m_nodes.resize( orig_size + Voxel::NV );
    std::copy_backward( m_nodes.begin()+m_offset+1, m_nodes.begin() + orig_size, m_nodes.end());

    for( Index k = 0 ; k < Voxel::NV ; ++k ) {
        m_nodes[ m_offset+1+k ].subSize = 1 ;
        m_nodes[ m_offset+1+k ].nLeafs  = 1 ;
    }
}

void Octree::SubTree::split ( Index leafIndex,  Index &offset, Index leafOffset )
{
    Node& cur = m_nodes[offset] ;
    cur.subSize += Voxel::NV ;
    cur.nLeafs  += Voxel::NV - 1 ;

    if( cur.subSize > Voxel::NV  + 1 ) {
        ++offset ;

        while( leafOffset + m_nodes[ offset ].nLeafs <= leafIndex ) {
            offset     += m_nodes[ offset ].subSize ;
            leafOffset += m_nodes[ offset ].nLeafs ;
        }

        split( leafIndex, offset, leafOffset ) ;
    }
}


}
