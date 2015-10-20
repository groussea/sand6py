#ifndef D6_GRID_HH
#define D6_GRID_HH

#include "MeshBase.hh"

namespace d6 {

class Grid ;

struct GridIterator
{
	typedef Vec3i Cell ;

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

	Index index() const ;

};

template < >
struct MeshTraits< Grid > {
	static constexpr Index NV = 8 ;

	typedef GridIterator CellIterator ;
};


class Grid : public MeshBase< Grid >
{
public:
	typedef MeshBase< Grid > Base ;

	typedef typename Base::Cell Cell ;
	typedef Vec3i 				Vertex ;

	Grid( const Vec& box, const Vec3i &res ) ;


	Index nNodes() const
	{ return (m_dim[0]+1) * (m_dim[1]+1) * (m_dim[2] + 1) ; }

	Index nCells() const
	{ return (m_dim[0]) * (m_dim[1]) * (m_dim[2]) ; }

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	void locate( const Vec &x, Location& loc ) const ;

	CellIterator cellBegin() const {
		return GridIterator( *this, Vec3i::Zero() ) ;
	}
	CellIterator cellEnd() const {
		return GridIterator( *this, Vec3i(m_dim[0],0,0) ) ;
	}

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
	Index cellIndex( const Cell& cell ) const
	{
		return (m_dim[2]) * (m_dim[1]) * cell[0]
			+  (m_dim[2]) * cell[1]
			+  cell[2] ;
	}

	void clamp_cell( Cell& cell ) const ;

	Vec firstCorner( const Cell &cell ) const
	{ return (cell.array().cast< Scalar >() * m_dx.array()).matrix() ; }

	Vec3i m_dim ;
	Vec   m_dx  ;
	Vec   m_idx  ;

	friend struct GridIterator ;
} ;

} //d6
#endif
