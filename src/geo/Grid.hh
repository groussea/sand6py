#ifndef D6_GRID_HH
#define D6_GRID_HH

#include "MeshBase.hh"

namespace d6 {

class Grid : public MeshBase< Grid >
{
public:
	typedef MeshBase< Grid > Base ;

	Grid( const Vec& box, const Vec3i &res ) ;


	size_t nNodes() const
	{ return (m_dim[0]+1) * (m_dim[1]+1) * (m_dim[2] + 1) ; }

	size_t nCells() const
	{ return (m_dim[0]) * (m_dim[1]) * (m_dim[2]) ; }

	Vec box() const
	{ return firstCorner( m_dim ) ; }

	void locate( const Vec &x, Location& loc ) ;

private:

	size_t nodeIndex( const Vec3i & node ) const
	{
		return (m_dim[2]+1) * (m_dim[1]+1) * node[0]
			+  (m_dim[2]+1) * node[1]
			+  node[2] ;
	}
	size_t cellIndex( const Vec3i & cell ) const
	{
		return (m_dim[2]) * (m_dim[1]) * cell[0]
			+  (m_dim[2]) * cell[1]
			+  cell[2] ;
	}

	void clamp_cell( Vec3i & cell ) const ;

	Vec firstCorner( const Vec3i &cell ) const
	{ return (cell.array().cast< Scalar >() * m_dx.array()).matrix() ; }

	Vec3i m_dim ;
	Vec   m_dx  ;
	Vec   m_idx  ;
} ;

} //d6
#endif
