#ifndef D6_MESH_BASE_HH
#define D6_MESH_BASE_HH

#include "utils/alg.hh"

namespace d6 {

template <typename Derived >
struct MeshTraits {
	static constexpr Index NV = 8 ;
};

template <typename Derived>
class MeshBase {

public:

	typedef MeshTraits< Derived > Traits ;
	typedef typename Traits::CellIterator CellIterator ;
	typedef typename CellIterator::Cell Cell ;
	static constexpr Index NV = Traits::NV ;

	typedef Eigen::Matrix<  Index, NV, 1> NodeList ;
	typedef Eigen::Matrix< Scalar, NV, 1> CoefList ;

	struct Location {
		Cell      cell  ; // Cell
		NodeList nodes  ; // Adjacent node indices
		CoefList coeffs ; // Coordinates in cell
	};

	Derived& derived()
	{ return static_cast< Derived& >( *this ) ; }
	const Derived& derived() const
	{ return static_cast< const Derived& >( *this ) ; }

	Index nNodes() const { return derived().nNodes() ; }
	Index nCells() const { return derived().nCells() ; }

	Vec box() const { return derived().box() ; }
	Vec clamp_point( const Vec& p ) const {
		return Vec::Zero().cwiseMax(p).cwiseMin(box()) ;
	}

	void locate( const Vec &x, Location& loc ) const {
		derived().locate( x, loc ) ;
	}

	CellIterator cellBegin() const { return derived().cellBegin() ; }
	CellIterator cellEnd()   const { return derived().cellEnd() ; }

} ;

} //d6

#endif
