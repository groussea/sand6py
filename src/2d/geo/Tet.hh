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

#ifndef D6_TET_HH
#define D6_TET_HH

#include "utils/alg.hh"
#include "geo/geo.fwd.hh"

namespace d6 {

struct Tet {

	static constexpr Index NV = 3 ;
	static constexpr Index NC = 3 ;
	static constexpr Index NE = 3 ;
	static constexpr Index NF = 3 ;

	typedef Eigen::Matrix< Scalar, NC, 1  > Coords ;
	typedef Eigen::Matrix< Scalar, NC, WD > Derivatives ;

	typedef Eigen::Matrix< Scalar, WD, NV > Vertices ;
	typedef Eigen::Matrix< Scalar, WD, Eigen::Dynamic > Points ;
	typedef Eigen::Matrix< Scalar, SD, Eigen::Dynamic > Frames ;

	struct {
		unsigned char sym  : 2 ;

		inline char sign   ( int d ) const {
			return 1 - 2*( ( sym&(1<<d) ) >> d )  ;
		}

	} geometry ;

	Vec origin ;  //!< 2D coords of rectangle corner
	Arr box    ;  //!< Dimensions of cell

	int orientation ; //!< 1 bit part + 3 bits symmetry

	Tet()
	: geometry{0}
	{}

	Vec pos( const Coords& coords ) const {

		Vertices vtx  ;
		compute_vertices( vtx ) ;

		Vec p = vtx*coords ;

		to_world( p ) ;
		return p + origin ;
	}

	Vec center() const
	{
		//TODO optimize
		return pos( Coords::Constant(1./NC) ) ;
	}

	Vec vertex( int cornerIndex ) const
	{
		Vec v ;
		offset( cornerIndex, v ) ;

		to_world(v) ;
		return v + origin ;
	}
	void vertexCoords( int cornerIndex, Coords& coords ) const {
		coords.setZero() ;
		coords[cornerIndex] = 1. ;
	}
	void edgeCenterCoords( int edgeIndex, Coords& coords ) const {
		coords.setConstant( 1./2 ) ;
		coords[edgeIndex] = 0. ;
	}

	Scalar volume() const {
		return box.prod() / 2 ;
	}

	void compute_coords( const Vec& pos, Coords & coords ) const {
		Vec local = pos - origin ;
		to_local( local );
		local_coords( local, coords );
	}
	void compute_derivatives( const Coords & coords, Derivatives& dc_dx ) const ;

	Index sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const ;

	void update_geometry( unsigned char rotation, int num ) ;

private:

	void compute_vertices( Vertices& vertices ) const ;

	void offset( int cornerIndex, Vec &v ) const ;

	void to_world( Vec& pos ) const ;
	void to_local( Vec& pos ) const ;

	void local_coords( const Vec& pos, Coords& coords ) const ;

} ;

template<>
struct QuadraturePoints< Tet, 1 >
{
	static constexpr Index NQ = 1 ;
	typedef Eigen::Matrix< Scalar, Tet::NC, 1> QuadPoint ;

	static void get( const Tet&, Index, QuadPoint& qp ) {
		qp.setConstant(1./3) ;
	}

	static Scalar weight( const Tet& geo, Index ) {
		return geo.volume() ;
	}

private:
	typedef Eigen::Matrix< Scalar, Tet::NC, NQ> QuadPoints ;
} ;

template<>
struct QuadraturePoints< Tet, 2 >
{
	static constexpr Index NQ = 3 ;
	typedef Eigen::Matrix< Scalar, Tet::NC, 1> QuadPoint ;

	static void get( const Tet&, Index k, QuadPoint& qp ) {
		static constexpr Scalar a = 1./6 ;
		static constexpr Scalar b = 2./3 ;

		qp.setConstant(a) ;
		qp[k] = b ;
	}

	static Scalar weight( const Tet& geo, Index ) {
		return geo.volume()/NQ ;
	}

} ;

template<>
struct QuadraturePoints< Tet, 4 >
{
	static constexpr Index NQ = 6 ;
	typedef Eigen::Matrix< Scalar, Tet::NC, 1> QuadPoint ;

	static void get( const Tet&, Index k, QuadPoint& qp ) {
		const Scalar a = k<3
		        ? 0.091576213509771
		        : 0.445948490915965 ;
		const Scalar b = k<3
		        ? 0.816847572980459
		        : 0.108103018168070 ;

		qp.setConstant(a) ;
		qp[k<3?k:(k-3)] = b ;
	}

	static Scalar weight( const Tet& geo, Index k ) {
		return ( k<3
		         ? 0.109951743655322
		         : 0.223381589678011
		           ) * geo.volume() ;
	}

} ;


} //d6

#endif
