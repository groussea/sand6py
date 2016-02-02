#ifndef D6_SHAPE_FUNCTION_BASE_HH
#define D6_SHAPE_FUNCTION_BASE_HH

#include "utils/alg.hh"
#include "geo.fwd.hh"

namespace d6 {

template <typename ShapeFunc>
struct ShapeFuncTraits
{
	static bool constexpr is_mesh_based = false ;
} ;


template< typename Derived >
struct ShapeFuncBase
{
	typedef ShapeFuncTraits< Derived > Traits ;

	typedef typename Traits::Location Location ;
	enum {
		NI = Traits::NI,
		NQ = Traits::NQ
	} ;
	static bool constexpr is_mesh_based = Traits::is_mesh_based ;
	typedef typename Traits::DOFDefinition DOFDefinition ;
	typedef typename Traits::QPIterator QPIterator ;

	typedef Eigen::Matrix<  Index, Traits::NI, 1  > NodeList ;
	typedef Eigen::Matrix< Scalar, Traits::NI, 1  > CoefList ;
	typedef Eigen::Matrix< Scalar, Traits::NI, WD > Derivatives ;

	struct Interpolation {
		NodeList nodes  ; // Adjacent node indices
		CoefList coeffs ; // Interpolation coeff for nodes
	};


	Index nDoF() const { return derived().nDoF() ; }

	const Derived& derived() const
	{ return static_cast< const Derived& >( *this ) ; }

	void interpolate( const Location& loc, Interpolation& itp ) const {
		derived().interpolate( loc, itp ) ;
	}
	void get_derivatives( const Location& loc, Derivatives& dc_dx ) const {
		derived().get_derivatives( loc, dc_dx ) ;
	}
	void list_nodes( const Location& loc, NodeList& list ) const {
		derived().list_nodes( loc, list ) ;
	}

	void compute_volumes( DynVec& volumes ) const {
		derived().compute_volumes( volumes ) ;
	}

	void all_dof_positions( DynMatW& vertices, std::vector<Index>& indexes, DynVec& totalVolumes  ) const
	{
		derived().all_dof_positions( vertices, indexes, totalVolumes ) ;
	}

	QPIterator qpBegin()
	{ return derived().qpBegin() ; }
	QPIterator qpEnd()
	{ return derived().qpEnd() ; }
};


} // d6

#endif
