#ifndef D6_MESH_SHAPE_FUNCTION_HH
#define D6_MESH_SHAPE_FUNCTION_HH

#include "ShapeFunctionBase.hh"

namespace d6 {

template< typename MeshT >

struct MeshQPIterator {

	typedef MeshBase< MeshT > MeshType ;

	explicit MeshQPIterator( const MeshType &mesh_, const typename MeshType::CellIterator& cIt )
		: mesh(mesh_), cellIt( cIt ), qp(0)
	{
		if( cellIt != mesh.cellEnd() ) update_cache();
	}

	bool operator==( const MeshQPIterator& o ) const
	{ return o.cellIt == cellIt && qp == o.qp ;}
	bool operator!=( const MeshQPIterator& o ) const
	{ return o.cellIt != cellIt || qp != o.qp ; }

	MeshQPIterator& operator++() {
		if( ++qp == MeshType::CellGeo::NQ ) {
			qp = 0 ;
			++cellIt ;
			if( cellIt != mesh.cellEnd() ) update_cache();
		}
	}

	Vec pos() const { return cached_qps.col(qp) ; }
	Scalar weight() const { return cached_qpw.col(qp) ; }
	typename MeshType::Location loc() const { return mesh.locate(pos()) ; }

private:
	const MeshType& mesh ;

	typename MeshType::CellIterator cellIt ;
	Index qp ;

	typename MeshType::Cell::QuadPoints  cached_qps ;
	typename MeshType::Cell::QuadWeights cached_qpw ;

	void update_cache() {
		typename MeshType::CellGeo geo ;
		mesh.get_geo( geo ) ;
		geo.get_qp( cached_qps, cached_qpw ) ;
	}
};


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
	typedef typename Base::QPIterator QPIterator ;

	MeshShapeFunc ( const MeshType & mesh )
		: m_mesh( mesh )
	{}

	void compute_volumes( DynVec& volumes ) const
	{
		volumes.resize( Base::nDoF() ) ;
		volumes.setZero() ;

		typename MeshType::CellGeo geo ;
		typename Base::NodeList nodes ;

		for( typename MeshType::CellIterator it = m_mesh.cellBegin() ; it != m_mesh.cellEnd() ; ++it ) {

			m_mesh.get_geo( *it, geo );
			const Scalar vol = geo.volume() / MeshType::NV ;

			list_nodes( *it, nodes );
			for( Index k = 0 ; k < MeshType::NV ; ++k ) {
				volumes[nodes[k]] += vol ;
			}
		}
	}

	void all_dof_positions( DynMatW& vertices, std::vector<Index>& indexes, DynVec& totalVolumes ) const
	{
		vertices.resize( WD, m_mesh.nNodes() );
		totalVolumes.resize( m_mesh.nNodes() );
		indexes.resize( Base::nDoF() );

		totalVolumes.setZero() ;

		typename MeshType::CellGeo cellGeo ;
		typename Base::NodeList meshNodes, shapeNodes ;

		for( typename MeshType::CellIterator it = m_mesh.cellBegin() ; it != m_mesh.cellEnd() ; ++it )
		{
			m_mesh.get_geo( *it, cellGeo ) ;

			list_nodes( *it, shapeNodes );
			Linear<MeshT>( m_mesh ).list_nodes( *it, meshNodes );

			const Scalar vol = cellGeo.volume() / meshNodes.rows() ;

			for( int k = 0 ; k < meshNodes.rows() ; ++k ) {
				indexes[shapeNodes[k]] = meshNodes[k] ;
				vertices.col( meshNodes[k] ) = cellGeo.vertex( k ) ;
				totalVolumes( meshNodes[k] ) += vol ;
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


	QPIterator qpBegin() const {
		return QPIterator( m_mesh, m_mesh.cellBegin() ) ;
	}

	QPIterator qpEnd() const {
		return QPIterator( m_mesh, m_mesh.cellEnd() ) ;
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
	typedef MeshQPIterator< MeshT > QPIterator ;

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
	typedef MeshQPIterator< MeshT > QPIterator ;

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

} ;


}

#endif
