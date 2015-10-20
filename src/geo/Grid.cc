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

void Grid::clamp_cell( Cell & cell ) const
{
	cell = Vec3i::Zero().cwiseMax(cell).cwiseMin(m_dim-Vec3i::Ones()) ;
}

void Grid::locate(const Vec &x, Location &loc) const
{
	loc.cell = (x.array()*m_idx.array()).matrix().cast< Index >();
	clamp_cell( loc.cell) ;

	const Vec coo = Vec::Ones().array() - ( x - firstCorner( loc.cell ) ).array() * m_idx.array() ;

	for( int i = 0 ; i < 2 ; ++i )
		for( int j = 0 ; j < 2 ; ++j )
			for( int k = 0 ; k < 2 ; ++k ) {
				const Cell corner (i,j,k) ;
				loc.nodes[ MK_INDEX(i,j,k) ] = nodeIndex( loc.cell + corner ) ;
				// c_i(x) = i + (1 - 2*i )*x = [ x if i=0, 1 + -x if i = 1 ]
				Vec coeffs = corner.cast< Scalar >().array() + ( Vec::Ones() - 2*corner.cast< Scalar >() ).array() * coo.array() ;
				loc.coeffs[ MK_INDEX(i,j,k) ] = coeffs[0]*coeffs[1]*coeffs[2] ;
			}

}

GridIterator& GridIterator::operator ++()
{
	++cell[2] ;
	if(cell[2] == grid.dim()[2]) {
		cell[2] = 0 ;
		++cell[1] ;
		if(cell[1] == grid.dim()[1]) {
			cell[1] = 0 ;
			++cell[0] ;
		}
	}

	return *this ;
}

Index GridIterator::index() const {
	return grid.cellIndex( cell ) ;
}

} //d6
