/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
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

#include "Tet.hh"

#include "Tensor.hh"

namespace d6 {

void Tet::to_world( Vec& pos ) const {
	pos[0] = geometry.sign(0)*pos[0] ;
	pos[1] = geometry.sign(1)*pos[1] ;

}
void Tet::to_local( Vec& pos ) const {

	pos[0] = geometry.sign(0)*pos[0] ;
	pos[1] = geometry.sign(1)*pos[1] ;
}

void Tet::update_geometry( unsigned char symmetry, int num )
{
	geometry.sym = num == 1 ? 3 : 0 ;
	geometry.sym ^= symmetry ;

	Vec o (-.5,-.5)	 ;
	Vec o_rot = o ;
	to_world( o_rot );

	origin.array() += box*( o_rot - o ).array() ;

	Vec rot_box = box ;
	to_local( rot_box ) ;
	box = rot_box.array().abs() ;
}



void Tet::offset( int cornerIndex, Vec &v ) const
{

	v.setZero() ;

	switch (cornerIndex) {
		case 1:
			v[0] = box[0] ;
			break;
		case 2:
			v[1] = box[1] ;
			break ;
	}
}

void Tet::compute_vertices( Vertices& vertices ) const
{
	vertices.setZero() ;
	vertices(0,1) = box[0] ;
	vertices(1,2) = box[1] ;
}

void Tet::local_coords(const Vec &pos, Coords &coords) const
{
	const Vec scaled = pos.array()/box ;

	coords[0] = 1 - scaled[0] - scaled[1] ;
	coords[1] = scaled[0] ;
	coords[2] = scaled[1] ;
}

void Tet::compute_derivatives( const Coords & coords, Derivatives& dc_dx ) const
{
	(void) coords ;

	dc_dx.setZero() ;
	dc_dx(0, 0) = -1./box[ 0 ] ;
	dc_dx(0, 1) = -1./box[ 1 ] ;
	dc_dx(1, 0) =  1./box[ 0 ] ;
	dc_dx(2, 1) =  1./box[ 1 ] ;

	//to_world

	dc_dx.col(0) = geometry.sign(0)*dc_dx.col(0) ;
	dc_dx.col(1) = geometry.sign(1)*dc_dx.col(1) ;
}

Index Tet::sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const
{

	Index p = start ;

	if( N == 1 ) {

		const Scalar a = 1./6 ;
		const Scalar b = 2./3 ;

		constexpr Index NQ = QuadraturePoints<Tet,2>::NQ ;

		(void) N ;
		const Vec subBox = box.array() / std::pow( NQ*2, 1./WD) ; //Nsub.array().cast< Scalar >() ;

		VecS frame ;
		tensor_view( frame ).set_diag( Vec( .25 * subBox.array() * subBox.array() ) ) ;

		for( Index k = 0 ; k < NQ ; ++k ) {
			Coords coords = Coords::Constant(a) ;
			coords[k] = b ;
			points.col(p) = pos(coords) ;
			frames.col(p) = frame ;
			++p ;
		}
	} else {

		Scalar min = box.minCoeff() ;

		VecWi Nsub ;
		for( int k = 0 ; k < WD ; ++ k)
			Nsub[k] = N * std::round( box[k] / min ) ;

		const Vec subBox = box.array() / Nsub.array().cast< Scalar >() ;

		Vec rot_box = subBox ;
		to_world( rot_box ) ;

		VecS frame ;
		tensor_view( frame ).set_diag( Vec( .25 * rot_box.array() * rot_box.array() ) ) ;

		Vec local ;
		Coords coords ;

		for( int i = 0 ; i < Nsub[0] ; ++i )
			for( int j = 0 ; j < Nsub[1] ; ++j ) {

					local = (Vec(i+.5,j+.5).array() * subBox.array()).matrix() ;
					local_coords( local, coords );
					if( coords.minCoeff() < 0 )
						continue ;

					points.col(p) = pos( coords )  ;
					frames.col(p) = frame ;
					++p ;
				}
	}

	return p - start;
}

} //d6
