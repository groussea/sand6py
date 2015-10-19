#ifndef D6_MESH_BASE_HH
#define D6_MESH_BASE_HH

#include "utils/alg.hh"

namespace d6 {

template <typename Derived >
struct MeshTraits {
	static constexpr unsigned NV = 8 ;
};

template <typename Derived>
class MeshBase {

public:

	typedef MeshTraits< Derived > Traits ;
	static constexpr unsigned NV = Traits::NV ;

	typedef Eigen::Matrix< size_t, NV, 1> NodeList ;
	typedef Eigen::Matrix< Scalar, NV, 1> CoefList ;

	struct Location {
		size_t   cidx   ; // Cell index
		NodeList nodes  ; // Adjacent node indices
		CoefList coeffs ; // Coordinates in cell
	};

	Derived& derived()
	{ return static_cast< Derived& >( *this ) ; }
	const Derived& derived() const
	{ return static_cast< const Derived& >( *this ) ; }

	size_t nNodes() const { return derived().nNodes() ; }
	size_t nCells() const { return derived().nCells() ; }

	Vec box() const { return derived().box() ; }
	Vec clamp_point( const Vec& p ) const {
		return Vec(0,0,0).cwiseMax(p).cwiseMin(box()) ;
	}

	void locate( const Vec &x, Location& loc ) {
		derived().locate( x, loc ) ;
	}

} ;

} //d6

#endif
