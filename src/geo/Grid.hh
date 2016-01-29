#ifndef D6_GRID_HH
#define D6_GRID_HH

#include "MeshBase.hh"
#include "Voxel.hh"

namespace d6 {

class Grid ;

struct GridIterator
{
	typedef MeshTraits< Grid > Traits ;
	typedef typename Traits::Cell    Cell ;
	typedef typename Traits::CellGeo CellGeo ;

	const Grid& grid ;
	Cell cell ;

	GridIterator( const Grid& g, const Cell& c )
		: grid(g), cell(c)
	{}

	GridIterator& operator++() ;

	bool operator==( const GridIterator& o ) const
	{
		return (cell == o.cell).all() ;
	}
	bool operator!=( const GridIterator& o ) const
	{
		return (cell != o.cell).any() ;
	}

	Cell operator*() const {
		return cell ;
	}

	Index index() const ;

};


class Grid : public MeshBase< Grid >
{
public:

	typedef MeshBase< Grid > Base ;

	typedef typename Base::Cell Cell ;
	typedef typename Base::CellGeo CellGeo ;
	typedef VecWi 				Vertex ;

	Grid( const Vec& box, const VecWi &res ) ;

	void set_box( const Vec& box ) ;

	Index nNodes() const
	{ return (m_dim+1).prod() ; }

	Index nCells() const
	{ return m_dim.prod(); }

	Index cellIndex( const Cell& cell ) const
	{
		Index idx =  (m_dim[1]) * cell[0]	+ cell[1] ;
		if( Cell::RowsAtCompileTime == 2 )
			return idx ;
		return idx*m_dim[2] + cell[2] ;
	}
	Index nodeIndex( const Vertex& node ) const
	{
		Index idx = (m_dim[1]+1) * node[0] + node[1] ;
		if( Cell::RowsAtCompileTime == 2 )
			return idx ;
		return idx*(m_dim[2]+1) + node[2] ;
	}

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	using Base::locate ;
	void locate( const Vec &x, Location& loc ) const
	{
		loc.coords = x.array()/m_dx.array() ;
		loc.cell = loc.coords.cast< Index >();
		clamp_cell( loc.cell) ;

		loc.coords -= loc.cell.cast< Scalar >().matrix() ;
	}

	CellIterator cellBegin() const {
		return GridIterator( *this, VecWi::Zero() ) ;
	}
	CellIterator cellEnd() const {
		VecWi cell = m_dim ;
		cell.tail<WD-1>().setZero() ;
		return GridIterator( *this, cell ) ;
	}

	void get_geo( const Cell &cell, CellGeo& geo ) const {
		get_corner( cell, geo.origin );
		geo.box = m_dx ;
	}

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & m_dim ;
		ar &  m_dx ;
	}

	void make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const ;

	Index nAdjacent( Index ) const {
		return NV ;
	}

	Vec nodePosition( const Vertex& node ) const
	{
		return firstCorner( node ) ;
	}

	const ArrWi& dim() const { return m_dim ; }
	const Arr&    dx() const { return  m_dx ; }

private:

	void get_corner( const Cell &cell, Vec& corner ) const {
		corner = (cell.array().cast< Scalar >() * m_dx.array()).matrix() ;
	}

	void clamp_cell( Cell& cell ) const {
		cell = Cell::Zero().max(cell).min(m_dim.array()-Cell::Ones()) ;
	}

	Vec firstCorner( const Cell &cell ) const
	{ Vec corner ; get_corner( cell, corner ) ; return corner ; }

	ArrWi m_dim ;
	Arr   m_dx  ;

	friend struct GridIterator ;
} ;

} //d6
#endif
