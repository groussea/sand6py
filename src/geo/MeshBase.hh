#ifndef D6_MESH_BASE_HH
#define D6_MESH_BASE_HH

#include "geo.fwd.hh"
#include "MeshTraits.hh"

#include "utils/alg.hh"

namespace d6 {

struct BoundaryInfo ;
struct BoundaryMapper ;
typedef std::vector< BoundaryInfo > BoundaryConditions ;


template <typename Derived>
class MeshBase {

public:

	typedef MeshTraits< Derived > Traits ;
	typedef typename Traits::CellIterator CellIterator ;
	typedef typename Traits::Cell    Cell  ;
	typedef typename Traits::Cells   Cells ;
	typedef typename Traits::CellGeo CellGeo ;

	enum {
		NV = Traits::NV ,
		NC = Traits::NC ,
		NQ = Traits::NQ
	} ;

	typedef Eigen::Matrix< Scalar, NC, 1 > Coords ;
	typedef Eigen::Matrix<  Index, NV, 1 > NodeList ;
	typedef Eigen::Matrix< Scalar, NV, 1 > CoefList ;
	typedef Eigen::Matrix< Scalar, NV, 3 > Derivatives ;


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

	Index cellIndex( const Cell& cell ) const
	{ return derived().cellIndex( cell ) ; }

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
	void get_derivatives( const Location& loc, Derivatives& dc_dx ) const {
		derived().get_derivatives( loc, dc_dx ) ;
	}

	void interpolate( const Vec &x, Interpolation& itp ) const {
		Location loc ;
		derived().locate( x, loc ) ;
		derived().interpolate( loc, itp ) ;
	}
	void list_nodes( const Cell& cell, NodeList& list ) const {
		derived().list_nodes( cell, list ) ;
	}

	CellIterator cellBegin() const { return derived().cellBegin() ; }
	CellIterator cellEnd()  const { return derived().cellEnd() ; }

	void get_geo( const Cell &cell, CellGeo& geo ) const {
		derived().get_geo( cell, geo ) ;
	}

	void make_bc( const BoundaryMapper& mapper, BoundaryConditions &bc ) const {
		derived().make_bc( mapper, bc ) ;
	}

	Index nAdjacent( Index ) const { return NV ; }
} ;

} //d6

#endif
