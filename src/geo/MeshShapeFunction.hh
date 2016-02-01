#ifndef D6_MESH_SHAPE_FUNCTION_HH
#define D6_MESH_SHAPE_FUNCTION_HH

#include "ShapeFunctionBase.hh"

namespace d6 {

// MeshShapeFunction
template< template< typename  > class Interp, typename MeshT >
struct MeshShapeFunc ;

template< template< typename > class Interp, typename MeshT >
struct ShapeFuncTraits< MeshShapeFunc< Interp, MeshT>  >
{
	typedef MeshBase<MeshT> MeshType ;
	typedef typename MeshType::Location Location ;

	typedef MeshType DOFDefinition ;
	static bool constexpr is_mesh_based = true ;
} ;

template< template< typename  > class Interp, typename MeshT >
struct MeshShapeFunc : public ShapeFuncBase< Interp<MeshT> >
{
	typedef ShapeFuncBase< Interp<MeshT> > Base ;
	typedef MeshBase<MeshT> MeshType ;

	MeshShapeFunc ( const MeshType & mesh )
		: m_mesh( mesh )
	{}

	void compute_volumes( DynVec& volumes ) const
	{

		for( typename MeshType::CellIterator it = m_mesh.cellBegin() ; it != m_mesh.cellEnd() ; ++it ) {
			typename MeshType::CellGeo geo ;
			typename MeshType::NodeList nodes ;

			m_mesh.get_geo( *it, geo );
			const Scalar vol = geo.volume() / MeshType::NV ;

			Base::list_nodes( *it, nodes );
			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				volumes[nodes[k]] += vol ;
			}
		}
	}

	void list_nodes( const typename MeshType::Cell& cell, typename Base::NodeList& list ) const {
		typename Base::Location loc ;
		loc.cell = cell ;
		Base::derived().list_nodes( loc, list ) ;
	}

	void locate( const Vec&x, typename Base::Location & loc ) const {
		m_mesh.locate( x, loc ) ;
	}
	typename Base::Location locate( const Vec& x ) const {
		typename Base::Location loc ;
		m_mesh.locate( x, loc ) ;
		return loc ;
	}

	const MeshT& mesh()  const {
		return m_mesh.derived() ;
	}

protected:
	const MeshType& m_mesh ;

};

// Linear shape func

template< typename MeshT >
struct ShapeFuncTraits< Linear<MeshT>  > : public ShapeFuncTraits< MeshShapeFunc< Linear, MeshT >  >
{
	typedef ShapeFuncTraits< MeshShapeFunc< Linear, MeshT > > Base ;
	enum {
		NI = Base::MeshType::NV,
		NQ = Base::MeshType::CellGeo::NQ
	};

} ;

template <typename MeshT>
struct Linear : public MeshShapeFunc< Linear, MeshT >
{
	typedef MeshShapeFunc< Linear, MeshT > Base ;
	typedef MeshBase<MeshT> MeshType ;
	typedef typename Base::Location Location ;

	Linear ( const MeshType & mesh )
		: Base( mesh )
	{}

	Index nDoF() const { return Base::mesh().nNodes() ; }

	void interpolate( const Location& loc, typename Base::Interpolation& itp ) const ;
	void get_derivatives( const Location& loc, typename Base::Derivatives& dc_dx ) const ;

	using Base::list_nodes ;
	void list_nodes( const Location& loc, typename Base::NodeList& list ) const {
		typename Base::Interpolation itp ;
		interpolate( loc, itp );
		list = itp.nodes ;
	}

	void all_dof_positions( DynMatW& vertices  ) const
	{
		// FIXME implement
	}

} ;

// DG Linear shape func

template< typename MeshT >
struct ShapeFuncTraits< DGLinear<MeshT>  > : public ShapeFuncTraits< MeshShapeFunc< DGLinear, MeshT >  >
{
	typedef ShapeFuncTraits< MeshShapeFunc< DGLinear, MeshT > > Base ;
	enum {
		NI = Base::MeshType::NV,
		NQ = Base::MeshType::CellGeo::NQ
	};

} ;

template <typename MeshT>
struct DGLinear : public MeshShapeFunc< DGLinear, MeshT >
{
	typedef MeshShapeFunc< DGLinear, MeshT > Base ;
	typedef MeshBase<MeshT> MeshType ;
	typedef typename Base::Location Location ;

	DGLinear ( const MeshType & mesh )
		: Base( mesh )
	{}

	Index nDoF() const { return Base::mesh().nCells() * MeshType::NV ; }

	void interpolate( const Location& loc, typename Base::Interpolation& itp ) const {
		Linear<MeshT>( Base::mesh() ) .interpolate( loc, itp ) ;
		list_nodes( loc, itp.nodes ) ;
	}
	void get_derivatives( const Location& loc, typename Base::Derivatives& dc_dx ) const {
		Linear<MeshT>( Base::mesh() ) .interpolate( loc, dc_dx ) ;
	}

	// Orderin:g all nodes from cell 0, all nodes from cell, ... etc
	using Base::list_nodes ;
	void list_nodes( const Location& loc, typename Base::NodeList& list ) const {
		const Index cellIdx = Base::mesh().cellIndex( loc.cell ) ;
		for( Index k = 0 ; k < MeshType::NV ; ++ k ) {
			list[k] = cellIdx * MeshType::NV + k ;
		}
	}

	void all_dof_positions( DynMatW& vertices  ) const
	{
		//FIXME implement
	}

} ;


}

#endif
