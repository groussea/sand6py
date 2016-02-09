#include "MeshShapeFunction.hh"

#include "Grid.hh"
#include "Voxel.hh"

#include "TetGrid.hh"
#include "Tet.hh"

namespace d6 {

template<>
void Linear<Grid>::interpolate( const Location& loc,
								typename Base::NodeList& nodes, typename Base::CoefList& coeffs ) const
{

	const Grid& g = mesh() ;

	for( int i = 0 ; i < 2 ; ++i ) {
		for( int j = 0 ; j < 2 ; ++j ) {
			const Grid::Cell corner( i,j );
			const int idx = Voxel::cornerIndex( i, j ) ;
			nodes[ idx ] = g.nodeIndex( loc.cell + corner ) ;
			coeffs[ idx ] = Voxel::cornerCoeff( corner, loc.coords );
		}
	}

}

template<>
void Linear<Grid>::get_derivatives( const Location& loc, typename Base::Derivatives& dc_dx ) const
{

	for( int i = 0 ; i < 2 ; ++i ) {
		for( int j = 0 ; j < 2 ; ++j ) {
				const Grid::Cell corner( i,j ) ;
				const int idx = Voxel::cornerIndex( i, j ) ;
				Voxel::getCornerDerivatives( corner, loc.coords, dc_dx.row( idx ) );
			}
	}

	for (int k = 0 ; k < WD ; ++k)
		dc_dx.col( k ) /= mesh().dx()[k] ;


}

template<>
void Linear<TetGrid>::interpolate( const Location& loc,
								typename Base::NodeList& nodes, typename Base::CoefList& coeffs ) const
{
	coeffs = loc.coords ;

	Tet geo ;
	mesh().get_geo(  loc.cell, geo ) ;

	for( Index k = 0 ; k < NI ; ++k ) {
		const VecWi v = ( geo.vertex(k).array() / mesh().dx().array()  + .5 ).cast<Index>() ;
		nodes[k] = mesh().nodeIndex( v ) ;
	}

}

template<>
void Linear<TetGrid>::get_derivatives( const Location& loc, typename Base::Derivatives& dc_dx ) const
{
	Tet geo ;
	mesh().get_geo(  loc.cell, geo ) ;

	geo.compute_derivatives( loc.coords, dc_dx ) ;
}

} // d6
