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

#ifndef D6_VOXEL_HH
#define D6_VOXEL_HH

#include "utils/alg.hh"
#include "geo/geo.fwd.hh"

namespace d6 {

struct Voxel {

	static constexpr Index NV = 4  ;
	static constexpr Index NC = WD ;
	static constexpr Index NE = 4 ;
	static constexpr Index NF = 4 ;

	typedef Eigen::Matrix< Scalar, NC, 1 > Coords ;
	typedef Eigen::Matrix< Scalar, WD, Eigen::Dynamic > Points ;
	typedef Eigen::Matrix< Scalar, SD, Eigen::Dynamic > Frames ;

	Vec origin ;  //!< 2D coords of first corner
	Vec box    ;  //!< Dimensions of cell

	template < typename Idx >
	static inline Idx cornerIndex( Idx i, Idx j ) {
		return (i << 0) + (j << 1) ;
	}
	static inline Index cornerIndex( const ArrWi& corner ) {
		return cornerIndex< Index >( corner[0], corner[1] ) ;
	}
	static inline Scalar cornerCoeff( const ArrWi& corner, const Coords &coords ) {
		// c_i(x) = i + (1 - 2*i )*( 1- x) = [ 1-x if i=0, 1 + (-1)(1-x) = x if i = 1 ]
		return ( corner.cast< Scalar >() + ( Coords::Ones().array() - 2*corner.cast< Scalar >() )
		        * ( Coords::Ones() - coords ).array() ).prod() ;
	}

	// \warning lacs scaling by 1./box
	template < typename Res >
	static void getCornerDerivatives( const ArrWi& corner, const Coords &coords, Res res ) {
		const Vec coeffs =  ( corner.cast< Scalar >() + ( Coords::Ones().array() - 2*corner.cast< Scalar >() )
		        * ( Coords::Ones() - coords ).array() ) ;
		Vec copy ;
		for( int k = 0 ; k < WD ; ++k ) {
			// d (c^k_i(x)) /dx _k = (2 * i - 1)
			copy = coeffs ; copy[k] = 2*corner[k] - 1 ;
			res[k] = copy.prod() ;
		}
	}

	Vec center() const {
		return origin + .5*box ;
	}

	Vec pos( const Coords& coords ) const {
		return origin.array() + coords.array()*box.array() ;
	}


	Vec vertex( int cornerIndex ) const {
		Coords coords ;
		vertexCoords( cornerIndex, coords );
		return origin + ( coords.array() * box.array() ).matrix() ;
	}

	static ArrWi corner( int cornerIndex ) {
		return ArrWi ( (cornerIndex&1)>>0, (cornerIndex&2)>>1 ) ;
	}

	void vertexCoords( int cornerIndex, Coords& coords ) const {
		coords = corner( cornerIndex ).cast< Scalar >() ;
	}

	Scalar volume() const { return box.prod() ; }

	Index sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const ;

} ;

template<>
struct QuadraturePoints< Voxel, 1 >
{
	static constexpr Index NQ = 1 ;
	typedef Eigen::Matrix< Scalar, Voxel::NC, 1> QuadPoint ;

	static void get( const Voxel&, Index, QuadPoint& qp ) {
		qp.setConstant(.5) ;
	}

	static Scalar weight( const Voxel& geo, Index ) {
		return geo.volume() ;
	}

private:
	typedef Eigen::Matrix< Scalar, Voxel::NC, NQ> QuadPoints ;
} ;

template<>
struct QuadraturePoints< Voxel, 2 >
{
	static constexpr Index NQ = 4 ;
	typedef Eigen::Matrix< Scalar, Voxel::NC, 1> QuadPoint ;

	static void get( const Voxel&, Index k, QuadPoint& qp ) {
		static const QuadPoints s_qps = Qps() ;

		qp = s_qps.col( k ) ;
	}

	static Scalar weight( const Voxel& geo, Index ) {
		return geo.volume()/NQ ;
	}

private:
	typedef Eigen::Matrix< Scalar, Voxel::NC, NQ> QuadPoints ;

	static QuadPoints Qps()  ;

} ;

} //d6

#endif


