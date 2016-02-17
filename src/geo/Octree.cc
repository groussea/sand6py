#include "Octree.hh"

#include "BoundaryInfo.hh"

#include <iostream>

namespace d6 {

OctreeIterator& OctreeIterator::operator ++()
{
	if( ++cell[2] == grid.subtree(cell).nLeafs() ) {
		cell[2] = 0 ;
		++cell[1] ;
		if(cell[1] == grid.dim()[1]) {
			cell[1] = 0 ;
			++cell[0] ;
		}
	}

	return *this ;
}

/// 2D ONLY

void Octree::boundaryInfo( const Location &loc, const BoundaryMapper& mapper, BoundaryInfo &info ) const
{
	constexpr Scalar eps = 1.e-6 ;

	const Vec &p = pos( loc ) ;
	const Vec &b = box() ;

	info.bc = BoundaryInfo::Interior ;

	if( p[0] < eps )
		info.combine( mapper( "left"), Vec(-1,0) ) ;
	if( p[0] > b[0] - eps )
		info.combine( mapper("right"), Vec( 1,0) ) ;

	if( p[1] < eps )
		info.combine( mapper("bottom"), Vec(0,-1) ) ;
	if( p[1] > b[1] - eps )
		info.combine( mapper(   "top"), Vec(0, 1) ) ;
}


/// END 2d

Index OctreeIterator::index() const
{
	return grid.cellIndex( cell ) ;
}

Octree::Octree(const Vec &box, const VecWi &res)
	: Base(),
	  m_maxDepth(4)
{
	m_dim = res ;
	set_box( box ) ;
	rebuild() ;
}

void Octree::set_box( const Vec& box )
{
	m_dx.array() = box.array()/m_dim.array().cast< Scalar >() ;
}

void Octree::rebuild()
{

	Index nTrees = m_dim.prod() ;
	m_trees.resize( nTrees );

	Index offset = 0 ;
	for( SubTree& tree : m_trees) {
		tree.offset( offset ) ;
		offset += tree.nLeafs() ;
	}

	Index maxRes = 1 << m_maxDepth ;

	m_nodeIds.clear();
	for( auto it = cellBegin() ; it != cellEnd() ; ++it ) {
		const Cell& cell = *it ;
		ArrWi hires ; Index size ;
		subtree( cell ).upres_cell( cell[WD], maxRes, hires, size );
		hires += maxRes * cell.head<WD>() ;

		for( Index k = 0 ; k < Voxel::NV ; ++k ) {
			const ArrWi corner = Voxel::corner( k ) ;
			const Index nodeIdx = hrNodeIndex( hires + size*corner, maxRes ) ;

			if( m_nodeIds.find(nodeIdx) == m_nodeIds.end() ) {
				m_nodeIds[nodeIdx] = m_nodeIds.size() ;
			}
		}
	}
}

void Octree::locate(const Vec &x, Location &loc) const
{
	Vec rel = x.array()/m_dx.array() ;
	loc.cell.head<WD>() = rel.cast< Index >();
	clamp_cell( loc.cell) ;

	rel -= loc.cell.head< WD >().cast< Scalar >().matrix() ;

	subtree(loc.cell).find( rel, loc.cell[WD], loc.coords );
}

void Octree::get_geo( const Cell &cell, CellGeo& geo ) const
{
	subtree(cell).compute_geo(cell[WD], geo);
	geo.box.array() *= m_dx ;
	geo.origin += firstCorner( cell.head<WD>() ) ;
}

bool Octree::split(const Cell &cell)
{
	return subtree(cell).split(cell[WD], m_maxDepth) ;
}

/////////////////////////
////////////////////////

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
}

void Octree::SubTree::compute_geo( Index leafIndex, Voxel &geo ) const  {
	SubCell cell ;
	subcell( leafIndex, cell, 0, 0 ) ;
	geo.origin = cell.coords.cast< Scalar >() / cell.res ;
	geo.box    = Vec::Constant(1./cell.res) ;
}

void Octree::SubTree::upres_cell( Index leafIndex, Index target_res, ArrWi& hires_cell, Index &size ) const
{
	SubCell cell ;
	subcell( leafIndex, cell, 0, 0 ) ;
	size = target_res / cell.res ;
	hires_cell = cell.coords * size ;
}


void Octree::SubTree::find( Index offset, Arr& rel, Index &leafIndex, Coords& coords ) const
{
//	std::cout << offset << " "
//			  << m_nodes[offset].subSize << " "
//			  << m_nodes[offset].nLeafs << " -> "
//			  << leafIndex << " "
//			  << std::endl ;

	if( m_nodes[ offset ].subSize == 1 ) {
		coords = rel ;
	} else {
		ArrWi sub = (2*rel).cast< Index >().min(1) ;
		Index subi = Voxel::cornerIndex( sub ) ;

		++offset ;
		for (Index i = 0 ; i < subi ; ++i) {
			leafIndex += m_nodes[ offset ].nLeafs ;
			offset    += m_nodes[ offset ].subSize ;
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
			leafOffset += m_nodes[ offset ].nLeafs ;
			offset     += m_nodes[ offset ].subSize ;
		}

		const ArrWi sub = Voxel::corner( subi ) ;
		subcell( leafIndex, cell, offset, leafOffset ) ;

		cell.coords += cell.res * sub ;
		cell.res *= 2 ;
	}
}

bool Octree::SubTree::split(Index leafIndex, Index maxDepth )
{
	Index offset = 0 ;
	bool ok = split( leafIndex, offset, 0, maxDepth ) ;

	if( !ok ) return false ;

	const std::size_t orig_size = m_nodes.size() ;
	m_nodes.resize( orig_size + Voxel::NV );
	std::copy_backward( m_nodes.begin()+offset+1, m_nodes.begin() + orig_size, m_nodes.end());

	for( Index k = 0 ; k < Voxel::NV ; ++k ) {
		m_nodes[ offset+1+k ].subSize = 1 ;
		m_nodes[ offset+1+k ].nLeafs  = 1 ;
	}

	return true ;
}

bool Octree::SubTree::split ( Index leafIndex,  Index &offset,
							  Index leafOffset, Index maxDepth )
{
	if( 0 == maxDepth ) return false ;

	Node& cur = m_nodes[offset] ;

//	std::cout << offset << " [spli] "
//			  << cur.subSize<< " "
//			  << cur.nLeafs<< " "
//			  << std::endl ;

	if( cur.subSize > 1 ) {
		++offset ;

		while( leafOffset + m_nodes[ offset ].nLeafs <= leafIndex ) {
			leafOffset += m_nodes[ offset ].nLeafs ;
			offset     += m_nodes[ offset ].subSize ;
		}

		if( ! split( leafIndex, offset, leafOffset, maxDepth-1 ) )
			return false ;

	}

	cur.subSize += Voxel::NV ;
	cur.nLeafs  += Voxel::NV - 1 ;

	return true ;

}


}
