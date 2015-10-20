#ifndef D6_GRID_HH
#define D6_GRID_HH

#include "MeshBase.hh"

namespace d6 {

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

	void clamp_cell( Cell& cell ) const ;

	Vec firstCorner( const Cell &cell ) const
	{ return (cell.array().cast< Scalar >() * m_dx.array()).matrix() ; }

	Vec3i m_dim ;
	Vec   m_dx  ;
	Vec   m_idx  ;
} ;

} //d6
#endif
