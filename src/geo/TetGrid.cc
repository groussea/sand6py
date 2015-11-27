#include "TetGrid.hh"

#include "BoundaryInfo.hh"

namespace d6 {

TetGridIterator& TetGridIterator::operator ++()
{
	++cell[3] ;
	if(cell[3] == 6) {
		cell[3] = 0 ;
		++cell[2] ;
		if(cell[2] == grid.dim()[2]) {
			cell[2] = 0 ;
			++cell[1] ;
			if(cell[1] == grid.dim()[1]) {
				cell[1] = 0 ;
				++cell[0] ;
			}
		}
	}

	return *this ;
}

Index TetGridIterator::index() const {
	return grid.cellIndex( cell ) ;
}

TetGrid::TetGrid(const Vec &box, const Vec3i &res)
	: Base()
{
	m_dim = res ;
	m_dx.array() = box.array()/res.array().cast< Scalar >() ;
	m_idx.array() = 1./m_dx.array() ;
}

void TetGrid::locate(const Vec &x, Location &loc) const
{
	Vec pos = x.array()*m_idx.array() ;

	loc.cell.segment<3>(0) =
			pos.cast< Index >().cwiseMax( Vec3i::Zero() ).cwiseMin(m_dim-Vec3i::Ones()) ;

//	pos -= loc.cell.segment<3>(0).cast< Scalar >() ;

	Tet geo ;
	for( loc.cell[3] = 0 ; loc.cell[3] < 6 ; ++loc.cell[3] )
	{
		get_geo(  loc.cell, geo ) ;
		geo.compute_coords( x, loc.coords );

		if( loc.coords.minCoeff() >= 0 )
			break ;
	}
}

void TetGrid::get_geo( const Cell &cell, CellGeo& geo ) const {
	const int color = (cell[0]%2) * 8 + (cell[1]%2) * 4 + (cell[2]%2) * 2 ;


	get_corner( cell.segment<3>(0), geo.origin );
	geo.box = m_dx ;

	geo.update_geometry( color, cell[3] ) ;

}

void TetGrid::interpolate(const Location &loc, Interpolation &itp) const
{

}

void TetGrid::get_derivatives( const Location& loc, Derivatives& dc_dx ) const
{
}


void TetGrid::make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const
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



}
