#ifndef D6_OCTREE_HH
#define D6_OCTREE_HH

#include "MeshBase.hh"
#include "Voxel.hh"

#include <unordered_map>

namespace d6 {

class Octree ;

struct OctreeIterator
{
	typedef MeshTraits< Octree > Traits ;
	typedef typename Traits::Cell    Cell ;
	typedef typename Traits::CellGeo CellGeo ;

	const Octree& grid ;
	Cell cell ;

	OctreeIterator( const Octree& g, const Cell& c )
		: grid(g), cell(c)
	{}

	OctreeIterator& operator++() ;

	bool operator==( const OctreeIterator& o ) const
	{
		return (cell == o.cell).all() ;
	}
	bool operator!=( const OctreeIterator& o ) const
	{
		return (cell != o.cell).any() ;
	}

	Cell operator*() const {
		return cell ;
	}

	Index index() const ;

};


class Octree : public MeshBase< Octree >
{
public:

	typedef MeshBase< Octree > Base ;

	typedef typename Base::Cell Cell ;
	typedef typename Base::CellGeo CellGeo ;
	typedef VecWi 				Vertex ;

	Octree( const Vec& box, const VecWi &res ) ;

	void rebuild() ;

	void set_box( const Vec& box ) ;

	Index nNodes() const {
		return m_nodeIds.size() ;
	}

	Index nCells() const {
		return m_trees.back().offset() + m_trees.back().nLeafs() ;
	}
	Index nSubCells( const Cell& cell ) const {
		return subtree(cell).nLeafs() ;
	}

	Index cellIndex( const Cell& cell ) const {
		return subtree(cell).offset() + cell[WD] ;
	}

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	using Base::locate ;
	void locate( const Vec &x, Location& loc ) const ;

	CellIterator cellBegin() const {
		return OctreeIterator( *this, Cell::Zero() ) ;
	}
	CellIterator cellEnd() const {
		Cell cell = Cell::Zero() ;
		cell[0] = m_dim[0] ;
		return OctreeIterator( *this, cell ) ;
	}

	void get_geo( const Cell &cell, CellGeo& geo ) const ;

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & m_dim ;
		ar &  m_dx ;
	}

	bool onBoundary( const Cell &cell ) const {
		//TODO remove false-positives
		return cell.head<WD>().minCoeff() == 0 || (m_dim - cell.head<WD>()).minCoeff() == 1 ;
	}

	void boundaryInfo( const Location &loc, const BoundaryMapper& mapper, BoundaryInfo &info ) const ;

	void clamp_cell( Cell& cell ) const {
		cell.head<WD>() = ArrWi::Zero().max(cell.head<WD>()).min(m_dim.array()-ArrWi::Ones()) ;
	}

	bool split( const Cell& cell ) ;

	const ArrWi& dim() const { return m_dim ; }
	const Arr&    dx() const { return  m_dx ; }

	struct SubTree {

		SubTree() ;

		void find( const Vec& pos, Index &leafIndex, Coords& coords ) const ;
		void compute_geo( Index leafIndex, Voxel &geo ) const ;
		void upres_cell(Index leafIndex, Index target_res, ArrWi& hires_cell, Index &size) const ;

		Index nLeafs() const ;

		Index offset() const { return m_offset ; }
		void offset( Index off ) { m_offset = off ; }

		bool split( Index leafIndex, Index maxDepth ) ;

	private:
		struct Node {
			Index subSize ;
			Index nLeafs  ;
		} ;

		Index m_offset ;
		std::vector< Node > m_nodes ;

		struct SubCell ;

		void find( Index offset, Arr& rel, Index &leafIndex, Coords& coords ) const ;

		void subcell ( Index leafIndex, SubCell& cell,
					   Index offset, Index leafOffset ) const  ;

		bool split( Index leafIndex,  Index &offset, Index leafOffset, Index maxDepth ) ;

	} ;

private:

	void get_corner( const ArrWi &cell, Vec& corner ) const {
		corner = (cell.head<WD>().array().cast< Scalar >() * m_dx.array()).matrix() ;
	}

	Vec firstCorner( const ArrWi &cell ) const
	{ Vec corner ; get_corner( cell, corner ) ; return corner ; }


	Index subtreeIndex( const Cell& cell ) const
	{
		Index idx =  (m_dim[1]) * cell[0]	+ cell[1] ;
		if( WD == 2 )
			return idx ;
		return idx*m_dim[2] + cell[2] ;
	}
	SubTree& subtree( const Cell& cell )
	{ return m_trees[ subtreeIndex(cell)] ; }
	const SubTree& subtree( const Cell& cell ) const
	{ return m_trees[ subtreeIndex(cell)] ; }

	Index hrNodeIndex( const ArrWi& node, Index res ) const
	{
		Index idx = (res*m_dim[1]+1) * node[0] + node[1] ;
		if( Cell::RowsAtCompileTime == 2 )
			return idx ;
		return idx*(res*m_dim[2]+1) + node[2] ;
	}

	ArrWi m_dim ;
	Arr   m_dx  ;
	Index m_maxDepth ;

	std::vector< SubTree > m_trees ;

	std::unordered_map< Index, Index > m_nodeIds ;

	friend struct OctreeIterator ;
} ;

} //d6
#endif
