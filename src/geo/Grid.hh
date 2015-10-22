#ifndef D6_GRID_HH
#define D6_GRID_HH

#include "MeshBase.hh"
#include "Voxel.hh"

namespace d6 {

class Grid ;

struct GridIterator
{
	typedef Vec3i Cell ;
	typedef Voxel CellGeo ;

	const Grid& grid ;
	Cell cell ;

	GridIterator( const Grid& g, const Cell& c )
		: grid(g), cell(c)
	{}

	GridIterator& operator++() ;

	bool operator==( const GridIterator& o ) const
	{
		return cell == o.cell ;
	}
	bool operator!=( const GridIterator& o ) const
	{
		return cell != o.cell ;
	}

	Cell operator*() const {
		return cell ;
	}

	Index index() const ;

};

template < >
struct MeshTraits< Grid > {
	typedef GridIterator CellIterator ;

	typedef typename CellIterator::CellGeo CellGeo ;
	static constexpr Index NV = CellGeo::NV ;
	static constexpr Index NC = CellGeo::NC ;

};


class Grid : public MeshBase< Grid >
{
public:
	typedef MeshBase< Grid > Base ;

	typedef typename Base::Cell Cell ;
	typedef typename Base::CellGeo CellGeo ;
	typedef Vec3i 				Vertex ;

	Grid( const Vec& box, const Vec3i &res ) ;


	Index nNodes() const
	{ return (m_dim[0]+1) * (m_dim[1]+1) * (m_dim[2] + 1) ; }

	Index nCells() const
	{ return (m_dim[0]) * (m_dim[1]) * (m_dim[2]) ; }

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	void locate( const Vec &x, Location& loc ) const ;

	using Base::interpolate ;
	void interpolate( const Location &loc , Interpolation& itp ) const ;

	CellIterator cellBegin() const {
		return GridIterator( *this, Vec3i::Zero() ) ;
	}
	CellIterator cellEnd() const {
		return GridIterator( *this, Vec3i(m_dim[0],0,0) ) ;
	}

	void get_geo( const Cell &cell, CellGeo& geo ) const {
		get_corner( cell, geo.corner );
		geo.box = m_dx ;
	}

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) {
		ar & m_dim ;
		ar &  m_dx ;
		ar & m_idx ;
	}

	const Vec3i& dim() const { return m_dim ; }
	const Vec&    dx() const { return  m_dx ; }
	const Vec&   idx() const { return m_idx ; }

private:


	Index nodeIndex( const Vertex& node ) const
	{
		return (m_dim[2]+1) * (m_dim[1]+1) * node[0]
			+  (m_dim[2]+1) * node[1]
			+  node[2] ;
	}
	Index cellIndex( const Cell& cell ) const
	{
		return (m_dim[2]) * (m_dim[1]) * cell[0]
			+  (m_dim[2]) * cell[1]
			+  cell[2] ;
	}

	void get_corner( const Cell &cell, Vec& corner ) const {
		corner = (cell.array().cast< Scalar >() * m_dx.array()).matrix() ;
	}

	void clamp_cell( Cell& cell ) const ;

	Vec firstCorner( const Cell &cell ) const
	{ Vec corner ; get_corner( cell, corner ) ; return corner ; }

	Vec3i m_dim ;
	Vec   m_dx  ;
	Vec   m_idx  ;

	friend struct GridIterator ;
} ;

} //d6
#endif
