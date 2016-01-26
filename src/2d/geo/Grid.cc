#include "Grid.hh"
#include "BoundaryInfo.hh"

namespace d6 {

GridIterator& GridIterator::operator ++()
{
	++cell[1] ;
	if(cell[1] == grid.dim()[1]) {
		cell[1] = 0 ;
		++cell[0] ;
	}

	return *this ;
}

Index GridIterator::index() const {
	return grid.cellIndex( cell ) ;
}


Grid::Grid(const Vec &box, const VecWi &res)
	: Base()
{
	m_dim = res ;
	set_box( box ) ;
}

void Grid::set_box( const Vec& box )
{
	m_dx.array() = box.array()/m_dim.array().cast< Scalar >() ;
	m_idx.array() = 1./m_dx.array() ;
}

void Grid::clamp_cell( Cell & cell ) const
{
	cell = Cell::Zero().max(cell).min(m_dim.array()-Cell::Ones()) ;
}

void Grid::locate(const Vec &x, Location &loc) const
{
	loc.coords = x.array()*m_idx.array() ;
	loc.cell = loc.coords.cast< Index >();
	clamp_cell( loc.cell) ;

	loc.coords -= loc.cell.cast< Scalar >().matrix() ;
}

void Grid::interpolate(const Location &loc, Interpolation &itp) const
{
	for( int i = 0 ; i < 2 ; ++i )
		for( int j = 0 ; j < 2 ; ++j ) {
			const Cell corner( i,j );
			const int idx = Voxel::cornerIndex( i, j ) ;
			itp.nodes[ idx ] = nodeIndex( loc.cell + corner ) ;
			itp.coeffs[ idx ] = Voxel::cornerCoeff( corner, loc.coords );
		}

}

void Grid::get_derivatives( const Location& loc, Derivatives& dc_dx ) const
{
	for( int i = 0 ; i < 2 ; ++i )
		for( int j = 0 ; j < 2 ; ++j ) {
				const Cell corner( i,j ) ;
				const int idx = Voxel::cornerIndex( i, j ) ;
				Voxel::getCornerDerivatives( corner, loc.coords, dc_dx.row( idx ) );
			}

	for (int k = 0 ; k < WD ; ++k)
		dc_dx.col( k ) *= m_idx[k] ;
}


void Grid::make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const
{
	for( Index i = 0 ; i <= m_dim[1] ; ++i ) {
		bc[ nodeIndex( Vertex(0       , i) ) ].combine( mapper( "left"), Vec(-1,0) ) ;
		bc[ nodeIndex( Vertex(m_dim[0], i) ) ].combine( mapper("right"), Vec( 1,0) ) ;
	}
	for( Index i = 0 ; i <= m_dim[0] ; ++i ) {
		bc[ nodeIndex( Vertex(i, 0       ) ) ].combine( mapper("bottom"), Vec(0,-1) ) ;
		bc[ nodeIndex( Vertex(i, m_dim[1]) ) ].combine( mapper(   "top"), Vec(0, 1) ) ;
	}
}

} //d6
