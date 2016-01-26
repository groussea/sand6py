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
	{ return (m_dim[0]+1) * (m_dim[1]+1) ; }

	Index nCells() const
	{ return (m_dim[0]) * (m_dim[1]); }

	Index cellIndex( const Cell& cell ) const
	{
		return (m_dim[1]) * cell[0]	+ cell[1] ;
	}

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	void locate( const Vec &x, Location& loc ) const ;

	using Base::interpolate ;
	void interpolate( const Location &loc , Interpolation& itp ) const ;

	void get_derivatives( const Location& loc, Derivatives& dc_dx ) const ;

	CellIterator cellBegin() const {
		return GridIterator( *this, VecWi::Zero() ) ;
	}
	CellIterator cellEnd() const {
		return GridIterator( *this, VecWi(m_dim[0],0) ) ;
	}

	void get_geo( const Cell &cell, CellGeo& geo ) const {
		get_corner( cell, geo.origin );
		geo.box = m_dx ;
	}

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & m_dim ;
		ar &  m_dx ;
		ar & m_idx ;
	}

	void list_nodes( const Cell& cell, NodeList& list ) const {
		Location loc { cell, Coords::Zero() } ;
		Interpolation itp ;
		interpolate( loc, itp );
		list = itp.nodes ;
	}

	void make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const ;

	Index nAdjacent( Index ) const {
		return NV ;
	}

	Vec nodePosition( const Vertex& node ) const
	{
		return firstCorner( node ) ;
	}

	Index nodeIndex( const Vertex& node ) const
	{
		return (m_dim[1]+1) * node[0] + node[1] ;
	}

	const ArrWi& dim() const { return m_dim ; }

private:

	const Arr&    dx() const { return  m_dx ; }
	const Arr&   idx() const { return m_idx ; }

	void get_corner( const Cell &cell, Vec& corner ) const {
		corner = (cell.array().cast< Scalar >() * m_dx.array()).matrix() ;
	}

	void clamp_cell( Cell& cell ) const ;

	Vec firstCorner( const Cell &cell ) const
	{ Vec corner ; get_corner( cell, corner ) ; return corner ; }

	ArrWi m_dim ;
	Arr   m_dx  ;
	Arr   m_idx  ;

	friend struct GridIterator ;
} ;

} //d6
#endif
