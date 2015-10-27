#include "Grid.hh"
#include "BoundaryInfo.hh"

namespace d6 {

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

void Grid::make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const
{
	bc.resize( nNodes() );
	for( Index i = 0 ; i <= m_dim[1] ; ++i ) {
		for( Index j = 0 ; j <= m_dim[2] ; ++j ) {
			bc[ nodeIndex( Vertex(0       , i, j) ) ].set( mapper( "left"), Vec(-1,0,0) ) ;
			bc[ nodeIndex( Vertex(m_dim[0], i, j) ) ].set( mapper("right"), Vec( 1,0,0) ) ;
		}
	}
	for( Index i = 0 ; i <= m_dim[0] ; ++i ) {
		for( Index j = 0 ; j <= m_dim[2] ; ++j ) {
			bc[ nodeIndex( Vertex(i, 0       , j) ) ].set( mapper("front"), Vec(0,-1,0) ) ;
			bc[ nodeIndex( Vertex(i, m_dim[1], j) ) ].set( mapper( "back"), Vec(0, 1,0) ) ;
		}
	}
	for( Index i = 0 ; i <= m_dim[0] ; ++i ) {
		for( Index j = 0 ; j <= m_dim[1] ; ++j ) {
			bc[ nodeIndex( Vertex(i, j, 0       ) ) ].set( mapper("bottom"), Vec(0,0,-1) ) ;
			bc[ nodeIndex( Vertex(i, j, m_dim[2]) ) ].set( mapper(   "top"), Vec(0,0, 1) ) ;
		}
	}
}

} //d6
