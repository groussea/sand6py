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
	loc.coords = x.array()*m_idx.array() ;
	loc.cell = loc.coords.cast< Index >();
	clamp_cell( loc.cell) ;

	loc.coords -= loc.cell.cast< Scalar >() ;
}

void Grid::interpolate(const Location &loc, Interpolation &itp) const
{
	for( int i = 0 ; i < 2 ; ++i )
		for( int j = 0 ; j < 2 ; ++j )
			for( int k = 0 ; k < 2 ; ++k ) {
				const Cell corner (i,j,k) ;
				const int idx = Voxel::cornerIndex( i, j, k ) ;
				itp.nodes[ idx ] = nodeIndex( loc.cell + corner ) ;
				itp.coeffs[ idx ] = Voxel::cornerCoeff( corner, loc.coords );
			}

}

void Grid::get_derivatives( const Location& loc, Derivatives& dc_dx ) const
{
	for( int i = 0 ; i < 2 ; ++i )
		for( int j = 0 ; j < 2 ; ++j )
			for( int k = 0 ; k < 2 ; ++k ) {
				const Cell corner (i,j,k) ;
				const int idx = Voxel::cornerIndex( i, j, k ) ;
				Voxel::getCornerDerivatives( corner, loc.coords, dc_dx.row( idx ) );
			}

	for (int k = 0 ; k < 3 ; ++k)
		dc_dx.col( k ) *= m_idx[k] ;
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
