#ifndef D6_TETGRID_HH
#define D6_TETGRID_HH

#include "MeshBase.hh"
#include "Tet.hh"

namespace d6 {

class TetGrid ;

struct TetGridIterator
{
	typedef MeshTraits< TetGrid > Traits ;
	typedef typename Traits::Cell    Cell ;
	typedef typename Traits::CellGeo CellGeo ;

	const TetGrid& grid ;
	Cell cell ;

	TetGridIterator( const TetGrid& g, const Cell& c )
		: grid(g), cell(c)
	{}

	TetGridIterator& operator++() ;

	bool operator==( const TetGridIterator& o ) const
	{
		return cell.isApprox( o.cell );
	}
	bool operator!=( const TetGridIterator& o ) const
	{
		return ! cell.isApprox( o.cell );
	}

	Cell operator*() const {
		return cell ;
	}

	Index index() const ;

};


class TetGrid : public MeshBase< TetGrid >
{
public:
	typedef MeshBase< TetGrid > Base ;

	typedef typename Base::Cell Cell ;
	typedef typename Base::CellGeo CellGeo ;
	typedef Vec3i 				Vertex ;

	TetGrid( const Vec& box, const Vec3i &res ) ;


	Index nNodes() const
	{ return (m_dim[0]+1) * (m_dim[1]+1) * (m_dim[2] + 1) ; }

	Index nCells() const
	{ return 6 * (m_dim[0]) * (m_dim[1]) * (m_dim[2]) ; }

	Index cellIndex( const Cell& cell ) const
	{
		return  (  (m_dim[2]) * (m_dim[1]) * cell[0]
				+  (m_dim[2]) * cell[1]
				+  cell[2] ) * 6 + cell[3] ;
	}

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	void locate( const Vec &x, Location& loc ) const ;

	using Base::interpolate ;
	void interpolate( const Location &loc , Interpolation& itp ) const ;

	void get_derivatives( const Location& loc, Derivatives& dc_dx ) const ;

	CellIterator cellBegin() const {
		return TetGridIterator( *this, Cell::Zero() ) ;
	}
	CellIterator cellEnd() const {
		return TetGridIterator( *this, Cell(m_dim[0],0,0,0) ) ;
	}

	void get_geo( const Cell &cell, CellGeo& geo ) const  ;

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

private:

	const Vec3i& dim() const { return m_dim ; }
	const Vec&    dx() const { return  m_dx ; }
	const Vec&   idx() const { return m_idx ; }



	Index nodeIndex( const Vertex& node ) const
	{
		return (m_dim[2]+1) * (m_dim[1]+1) * node[0]
			+  (m_dim[2]+1) * node[1]
			+  node[2] ;
	}

	void get_corner( const Vec3i &cell, Vec& corner ) const {
		corner = (cell.array().cast< Scalar >() * m_dx.array()).matrix() ;
	}

	Vec firstCorner( const Vec3i &cell ) const
	{ Vec corner ; get_corner( cell, corner ) ; return corner ; }

	Vec3i m_dim ;
	Vec   m_dx  ;
	Vec   m_idx  ;

	friend struct TetGridIterator ;
};


} //d6

#endif
