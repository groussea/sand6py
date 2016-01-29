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

Index GridIterator::index() const
{
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
