#ifndef D6_VOXEL_HH
#define D6_VOXEL_HH

#include "utils/alg.hh"

namespace d6 {

struct Voxel {

	static constexpr Index NV = 8 ;
	static constexpr Index NC = 3 ;
	static constexpr Index NQ = 8 ;

	typedef Eigen::Matrix< Scalar, NC, 1 > Coords ;
	typedef Eigen::Matrix< Scalar, NC, NQ> QuadPoints ;
	typedef Eigen::Matrix< Scalar,  1, NQ> QuadWeights ;
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

	// \warning lacs scaling by 1./box
	template < typename Res >
	static void getCornerDerivatives( const Vec3i& corner, const Coords &coords, Res res ) {
		const Vec coeffs =  ( corner.cast< Scalar >().array() + ( Coords::Ones() - 2*corner.cast< Scalar >() ).array()
				* ( Coords::Ones() - coords ).array() ) ;
		Vec copy ;
		for( int k = 0 ; k < 3 ; ++k ) {
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
		const Vec3i corner ( (cornerIndex&1)>>0, (cornerIndex&2)>>1, (cornerIndex&4)>>2 ) ;
		return origin + ( corner.cast< Scalar >().array() * box.array() ).matrix() ;
	}

	Scalar volume() const { return box.prod() ; }

	Index sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const ;

	static QuadPoints Qps()  ;

	void get_qp( QuadPoints& qp, QuadWeights& weights ) const {
		static const QuadPoints s_qps = Qps() ;
		qp = s_qps ;
		weights.setConstant( volume() / NQ ) ;
	}

} ;

} //d6

#endif


