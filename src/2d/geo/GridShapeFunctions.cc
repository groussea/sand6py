/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "MeshShapeFunction.hh"
#include "P2ShapeFunction.hh"

#include "TetGrid.hh"
#include "Tet.hh"

#include "Grid.hh"

namespace d6 {


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

template<>
void Edgewise<Grid>::list_nodes( const Location& loc, typename Base::NodeList& nodes ) const
{
	const Grid& grid = Base::mesh().derived() ;
	nodes[ 0 ] = loc.cell[0] * (grid.dim()[1]+1) + loc.cell[1]   ;
	nodes[ 1 ] = loc.cell[0] * (grid.dim()[1]+1) + loc.cell[1]+1 ;
	const Index off = grid.dim()[0] * (grid.dim()[1]+1) ;
	nodes[ 2 ] = loc.cell[1] * (grid.dim()[0]+1) + loc.cell[0]   + off;
	nodes[ 3 ] = loc.cell[1] * (grid.dim()[0]+1) + loc.cell[0]+1 + off;
}

template<>
void Edgewise<Grid>::dof_coeffs( const typename MeshType::Coords& coords, typename Base::CoefList& coeffs ) const
{
	coeffs[0] = 1.-coords[1] ;
	coeffs[1] =    coords[1] ;
	coeffs[2] = 1.-coords[0] ;
	coeffs[3] =    coords[0] ;
}

template<>
void Edgewise<Grid>::get_derivatives( const Location&, typename Base::Derivatives& dc_dx ) const
{
	dc_dx(0,0) =  0;
	dc_dx(0,1) = -1;
	dc_dx(1,0) =  0;
	dc_dx(1,1) =  1;

	dc_dx(2,0) = -1;
	dc_dx(2,1) =  0;
	dc_dx(3,0) =  1;
	dc_dx(3,1) =  0;

	for (int k = 0 ; k < WD ; ++k)
		dc_dx.col( k ) /= mesh().dx()[k] ;
}

} // d6
