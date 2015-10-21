#ifndef D6_MESH_BASE_HH
#define D6_MESH_BASE_HH

#include "utils/alg.hh"

namespace d6 {

template <typename Derived >
struct MeshTraits {
};

template <typename Derived>
class MeshBase {

public:

	typedef MeshTraits< Derived > Traits ;
	typedef typename Traits::CellIterator CellIterator ;
	typedef typename CellIterator::Cell    Cell ;
	typedef typename Traits::CellGeo CellGeo ;

	static constexpr Index NV = Traits::NV ;
	static constexpr Index NC = Traits::NC ;

	typedef Eigen::Matrix< Scalar, NC, 1> Coords ;
	typedef Eigen::Matrix<  Index, NV, 1> NodeList ;
	typedef Eigen::Matrix< Scalar, NV, 1> CoefList ;

	struct Location {
		Cell      cell  ; // Cell
		Coords   coords ; // Coordinates in cell
	};
	struct Interpolation {
		NodeList nodes  ; // Adjacent node indices
		CoefList coeffs ; // Interpolation coeff for nodes
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
	void interpolate( const Location& loc, Interpolation& itp ) const {
		derived().interpolate( loc, itp ) ;
	}
	void interpolate( const Vec &x, Interpolation& itp ) const {
		Location loc ;
		derived().locate( x, loc ) ;
		derived().interpolate( loc, itp ) ;
	}

	CellIterator cellBegin() const { return derived().cellBegin() ; }
	CellIterator cellEnd()  const { return derived().cellEnd() ; }

	void get_geo( const Cell &cell, CellGeo& geo ) const {
		derived().get_geo( cell, geo ) ;
	}
} ;

} //d6

#endif
