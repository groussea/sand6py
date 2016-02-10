#include "TetGrid.hh"

#include "BoundaryInfo.hh"

namespace d6 {

TetGridIterator& TetGridIterator::operator ++()
{
	++cell[2] ;
	if(cell[2] == TetGrid::Nsub) {
		cell[2] = 0 ;
		++cell[1] ;
		if(cell[1] == grid.dim()[1]) {
			cell[1] = 0 ;
			++cell[0] ;
		}
	}

	return *this ;
}

Index TetGridIterator::index() const {
	return grid.cellIndex( cell ) ;
}

TetGrid::TetGrid(const Vec &box, const VecWi &res)
	: Base()
{
	m_dim = res ;
	m_dx = box.array()/res.array().cast< Scalar >() ;
}

Index TetGrid::nEdges() const {
	return ( m_dim*(m_dim+1) ).sum() + m_dim.prod() ;
}

void TetGrid::locate(const Vec &x, Location &loc) const
{
	Vec pos = x.array() / m_dx ;

	loc.cell.segment<WD>(0) =
			pos.cast< Index >().array().max( ArrWi::Zero() ).min(m_dim-ArrWi::Ones()) ;

	Tet geo ;
	for( loc.cell[WD] = 0 ; loc.cell[WD] < Nsub ; ++loc.cell[WD] )
	{
		get_geo(  loc.cell, geo ) ;
		geo.compute_coords( x, loc.coords );

		if( loc.coords.minCoeff() >= -1.e-12 )
			break ;
	}
}

void TetGrid::get_geo( const Cell &cell, CellGeo& geo ) const {
	const int color = (cell[0]%2) + (cell[1]%2) * 2 ;

	get_corner( cell.head<WD>(), geo.origin );
	geo.box = m_dx ;

	geo.update_geometry( color, cell[WD] ) ;

}

void TetGrid::boundaryInfo( const Location &loc, const BoundaryMapper& mapper, BoundaryInfo &info ) const
{
	constexpr Scalar eps = 1.e-6 ;

	const Vec &p = pos( loc ) ;
	const Vec &b = box() ;

	info.bc = BoundaryInfo::Interior ;

	if( p[0] < eps )
		info.combine( mapper( "left"), Vec(-1,0) ) ;
	if( p[0] > b[0] - eps )
		info.combine( mapper("right"), Vec( 1,0) ) ;

	if( p[1] < eps )
		info.combine( mapper("bottom"), Vec(0,-1) ) ;
	if( p[1] > b[1] - eps )
		info.combine( mapper(   "top"), Vec(0, 1) ) ;
}



}
