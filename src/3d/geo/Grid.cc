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

void Grid::boundaryInfo( const Location &loc, const BoundaryMapper& mapper, BoundaryInfo &info ) const
{
	constexpr Scalar eps = 1.e-6 ;

	const Vec &p = pos( loc ) ;
	const Vec &b = box() ;

	info.bc = BoundaryInfo::Interior ;

	if( p[0] < eps )
		info.combine( mapper( "left"), Vec(-1,0,0) ) ;
	if( p[0] > b[0] - eps )
		info.combine( mapper("right"), Vec( 1,0,0) ) ;

	if( p[1] < eps )
		info.combine( mapper("front"), Vec(0,-1,0) ) ;
	if( p[1] > b[1] - eps )
		info.combine( mapper( "back"), Vec(0, 1,0) ) ;

	if( p[2] < eps )
		info.combine( mapper("bottom"), Vec(0,0,-1) ) ;
	if( p[2] > b[2] - eps )
		info.combine( mapper(   "top"), Vec(0,0, 1) ) ;
}

} //d6
