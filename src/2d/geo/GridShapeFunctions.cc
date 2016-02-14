#include "MeshShapeFunction.hh"
#include "P2ShapeFunction.hh"

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


template<>
void P2<TetGrid>::list_nodes( const Location& loc, typename Base::NodeList& nodes ) const
{
	Tet geo ;
	mesh().get_geo(  loc.cell, geo ) ;

	TetGrid::Coords coords ;

	// Nodes correspond to those of double-res grid
	for( Index k = 0 ; k < NI ; ++k ) {
		dof_coords( geo, k, coords );
		const VecWi v = ( 2*geo.pos(coords).array() / mesh().dx().array()  + .5 ).cast<Index>() ;
		nodes[k] = (2*mesh().dim()[1]+1) * v[0] + v[1] ;
	}
}

template<>
void P2<TetGrid>::dof_coeffs( const typename MeshType::Coords& coords, typename Base::CoefList& coeffs) const
{
	// Order shoudl be consistent with locate_dof
	//   1
	//   |\  5
	//   | +
	//	 |  \  0
	// 3 *  /
	//   | +
	//   |/  4
	//   2

	const Scalar &c0 = coords[0] ;
	const Scalar &c1 = coords[1] ;
	const Scalar &c2 = coords[2] ;

	coeffs[0] = c0 * ( 2*c0 - 1 ) ;
	coeffs[1] = c1 * ( 2*c1 - 1 ) ;
	coeffs[2] = c2 * ( 2*c2 - 1 ) ;
	coeffs[3] = 4 * (c1*c2) ;
	coeffs[4] = 4 * (c0*c2) ;
	coeffs[5] = 4 * (c1*c0) ;
}

template<>
void P2<TetGrid>::dof_coeffs_tpz( const typename MeshType::Coords& coords, typename Base::CoefList& coeffs) const
{
	const Scalar &c0 = coords[0] ;
	const Scalar &c1 = coords[1] ;
	const Scalar &c2 = coords[2] ;

	coeffs[0] = std::max( 0., 2*c0 -1 ) ;
	coeffs[1] = std::max( 0., 2*c1 -1 ) ;
	coeffs[2] = std::max( 0., 2*c2 -1 ) ;

	coeffs[3] = std::max(0., 1 - std::max( 2*c0, c0 + std::fabs(c1-c2) )) ;
	coeffs[4] = std::max(0., 1 - std::max( 2*c1, c1 + std::fabs(c0-c2) )) ;
	coeffs[5] = std::max(0., 1 - std::max( 2*c2, c2 + std::fabs(c1-c0) )) ;
}

template<>
void P2<TetGrid>::get_derivatives( const Location& loc, typename Base::Derivatives& dc_dx ) const
{
	Tet geo ;
	mesh().get_geo(  loc.cell, geo ) ;

	typename Tet::Derivatives lin_dx ;
	geo.compute_derivatives( loc.coords, lin_dx ) ;

	const Scalar &c0 = loc.coords[0] ;
	const Scalar &c1 = loc.coords[1] ;
	const Scalar &c2 = loc.coords[2] ;

	dc_dx.row(0) = ( 4*c0 - 1 ) * lin_dx.row(0) ;
	dc_dx.row(1) = ( 4*c1 - 1 ) * lin_dx.row(1) ;
	dc_dx.row(2) = ( 4*c2 - 1 ) * lin_dx.row(2) ;
	dc_dx.row(3) = 4 * ( c2 * lin_dx.row(1) + c1 * lin_dx.row(2) ) ;
	dc_dx.row(4) = 4 * ( c0 * lin_dx.row(2) + c2 * lin_dx.row(0) ) ;
	dc_dx.row(5) = 4 * ( c1 * lin_dx.row(0) + c0 * lin_dx.row(1) ) ;

}

template<>
void P2<TetGrid>::build_visu_mesh( DynMatW& vertices, DynMati& indices ) const
{
	Grid tg ( Base::mesh().box(), 2 * Base::mesh().dim() ) ;
	Linear< Grid > lin ( tg ) ;

	lin.build_visu_mesh( vertices, indices );
}

} // d6
