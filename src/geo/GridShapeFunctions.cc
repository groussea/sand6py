#include "MeshShapeFunction.hh"

#include "Grid.hh"
#include "Voxel.hh"

namespace d6 {

template<>
void Linear<Grid>::interpolate( const Location& loc, typename Base::Interpolation& itp ) const
{

	const Grid& g = mesh() ;

	for( int i = 0 ; i < 2 ; ++i ) {
		for( int j = 0 ; j < 2 ; ++j ) {
			const Grid::Cell corner( i,j );
			const int idx = Voxel::cornerIndex( i, j ) ;
			itp.nodes[ idx ] = g.nodeIndex( loc.cell + corner ) ;
			itp.coeffs[ idx ] = Voxel::cornerCoeff( corner, loc.coords );
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

} // d6
