#include "Grid.hh"

#define MK_INDEX(i,j,k) (((i)<<2) + ((j)<<1) + (k))

namespace d6 {

Grid::Grid(const Vec &box, const Vec3i &res)
	: Base()
{
	m_dim = res ;
	m_dx.array() = box.array()/res.array().cast< Scalar >() ;
	m_idx.array() = 1./m_dx.array() ;
}

void Grid::clamp_cell( Vec3i & cell ) const
{
	cell = Vec3i(0,0,0).cwiseMax(cell).cwiseMin(m_dim-Vec3i(1,1,1)) ;
}

void Grid::locate(const Vec &x, Location &loc)
{
	Vec3i cell = (x.array()*m_idx.array()).matrix().cast< int >();
	clamp_cell(cell) ;

	loc.cidx = cellIndex( cell ) ;

	const Vec coo = Vec::Ones().array() - ( x - firstCorner( cell ) ).array() * m_idx.array() ;

	for( int i = 0 ; i < 2 ; ++i )
		for( int j = 0 ; j < 2 ; ++j )
			for( int k = 0 ; k < 2 ; ++k ) {
				const Vec3i corner (i,j,k) ;
				loc.nodes[ MK_INDEX(i,j,k) ] = nodeIndex( cell + corner ) ;
				// c_i(x) = i + (1 - 2*i )*x = [ x if i=0, 1 + -x if i = 1 ]
				Vec coeffs = corner.cast< Scalar >().array() + ( Vec::Ones() - 2*corner.cast< Scalar >() ).array() * coo.array() ;
				loc.coeffs[ MK_INDEX(i,j,k) ] = coeffs[0]*coeffs[1]*coeffs[2] ;
			}

}

} //d6
