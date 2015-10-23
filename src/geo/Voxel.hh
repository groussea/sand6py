#ifndef D6_VOXEL_HH
#define D6_VOXEL_HH

#include "utils/alg.hh"
#include <vector>

namespace d6 {

struct Voxel {

	static constexpr Index NV = 8 ;
	static constexpr Index NC = 3 ;

	typedef Eigen::Matrix< Scalar, NC, 1> Coords ;
	typedef Eigen::Matrix< Scalar, 3, Eigen::Dynamic > Points ;
	typedef Eigen::Matrix< Scalar, 6, Eigen::Dynamic > Frames ;

	Vec origin ;  //!< 3D coords of first corner
	Vec box    ;  //!< Dimensions of cell

	template < typename Idx >
	static inline Idx cornerIndex( Idx i, Idx j, Idx k ) {
		return (i << 0) + (j << 1) + (k << 2) ;
	}
	static inline Index cornerIndex( const Vec3i& corner ) {
		return cornerIndex< Index >( corner[0], corner[1], corner[2] ) ;
	}
	static inline Scalar cornerCoeff( const Vec3i& corner, const Coords &coords ) {
		// c_i(x) = i + (1 - 2*i )*( 1- x) = [ 1-x if i=0, 1 + (-1)(1-x) = x if i = 1 ]
		return ( corner.cast< Scalar >().array() + ( Coords::Ones() - 2*corner.cast< Scalar >() ).array()
				* ( Coords::Ones() - coords ).array() ).prod() ;
	}

	Vec center() const {
		return origin + .5*box ;
	}

	Vec vertex( int cornerIndex ) const {
		const Vec3i corner ( (cornerIndex&1)>>0, (cornerIndex&2)>>1, (cornerIndex&4)>>2 ) ;
		return origin + ( corner.cast< Scalar >().array() * box.array() ).matrix() ;
	}

	Scalar volume() const { return box.prod() ; }

	Index sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const ;

} ;

} //d6

#endif


