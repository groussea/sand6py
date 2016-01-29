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
	set_box( box ) ;
}

void Grid::set_box( const Vec& box )
{
	m_dx = box.array()/m_dim.array().cast< Scalar >() ;
}


void Grid::make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const
{
	for( Index i = 0 ; i <= m_dim[1] ; ++i ) {
		for( Index j = 0 ; j <= m_dim[2] ; ++j ) {
			bc[ nodeIndex( Vertex(0       , i, j) ) ].combine( mapper( "left"), Vec(-1,0,0) ) ;
			bc[ nodeIndex( Vertex(m_dim[0], i, j) ) ].combine( mapper("right"), Vec( 1,0,0) ) ;
		}
	}
	for( Index i = 0 ; i <= m_dim[0] ; ++i ) {
		for( Index j = 0 ; j <= m_dim[2] ; ++j ) {
			bc[ nodeIndex( Vertex(i, 0       , j) ) ].combine( mapper("front"), Vec(0,-1,0) ) ;
			bc[ nodeIndex( Vertex(i, m_dim[1], j) ) ].combine( mapper( "back"), Vec(0, 1,0) ) ;
		}
	}
	for( Index i = 0 ; i <= m_dim[0] ; ++i ) {
		for( Index j = 0 ; j <= m_dim[1] ; ++j ) {
			bc[ nodeIndex( Vertex(i, j, 0       ) ) ].combine( mapper("bottom"), Vec(0,0,-1) ) ;
			bc[ nodeIndex( Vertex(i, j, m_dim[2]) ) ].combine( mapper(   "top"), Vec(0,0, 1) ) ;
		}
	}
}

} //d6
