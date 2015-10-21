#ifndef D6_VOXEL_HH
#define D6_VOXEL_HH

#include "utils/alg.hh"

namespace d6 {

struct Voxel {

	static constexpr Index NV = 8 ;
	static constexpr Index NC = 3 ;

	typedef Eigen::Matrix< Scalar, NC, 1> Coords ;

	Vec corner ;
	Vec box    ;

	template < typename Idx >
	static inline Idx cornerIndex( Idx i, Idx j, Idx k ) {
		return (i << 2) + (j << 1) + k ;
	}
	static inline Index cornerIndex( const Vec3i& corner ) {
		return cornerIndex< Index >( corner[0], corner[1], corner[2] ) ;
	}
	static inline Scalar cornerCoeff( const Vec3i& corner, const Coords &coords ) {
		// c_i(x) = i + (1 - 2*i )*( 1- x) = [ 1-x if i=0, 1 + (-1)(1-x) = x if i = 1 ]
		return ( corner.cast< Scalar >().array() + ( Coords::Ones() - 2*corner.cast< Scalar >() ).array()
				* ( Coords::Ones() - coords ).array() ).prod() ;
	}

} ;

} //d6

#endif


