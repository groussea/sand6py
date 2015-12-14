#ifndef D6_FORM_BUILDER_IMPL_HH
#define D6_FORM_BUILDER_IMPL_HH

#include "FormBuilder.hh"

#include "geo/MeshImpl.hh"
#include "geo/Particles.hh"

namespace d6 {

template < typename Func >
void FormBuilder::integrate_qp( const typename MeshType::Cells& cells, Func func ) const
{
	integrate_qp( cells.begin(), cells.end(), func ) ;
}

template < typename CellIterator, typename Func >
void FormBuilder::integrate_qp( const CellIterator& cellBegin, const CellIterator& cellEnd, Func func ) const
{
	typename MeshType::Location loc ;
	typename MeshType::Interpolation itp ;
	typename MeshType::Derivatives dc_dx ;

	Eigen::Matrix< Scalar, MeshType::NC, MeshType::NQ > qp ;
	Eigen::Matrix< Scalar,            1, MeshType::NQ > qp_weights ;

	typename MeshType::CellGeo cellGeo ;

	for( CellIterator it = cellBegin ; it !=  cellEnd ; ++it )
	{
		const typename MeshType::Cell& cell = *it ;
		loc.cell = cell ;
		m_mesh.get_geo( cell, cellGeo );

		cellGeo.get_qp( qp, qp_weights ) ;

		for ( int q = 0 ; q < MeshType::NQ ; ++q ) {
			loc.coords = qp.col(q) ;

			m_mesh.interpolate( loc, itp );
			m_mesh.get_derivatives( loc, dc_dx );

			func( qp_weights[q], loc, itp, dc_dx ) ;
		}
	}
}

template < typename Func >
void FormBuilder::integrate_node( const typename MeshType::Cells& cells, Func func ) const
{
	integrate_node( cells.begin(), cells.end(), func ) ;
}

template < typename CellIterator, typename Func >
void FormBuilder::integrate_node( const CellIterator& cellBegin, const CellIterator& cellEnd, Func func ) const
{
	typename MeshType::Location loc ;
	typename MeshType::Interpolation itp ;

	typename MeshType::CellGeo cellGeo ;

	itp.coeffs.setZero() ;

	for( CellIterator it = cellBegin ; it !=  cellEnd ; ++it )
	{
		const typename MeshType::Cell& cell = *it ;
		loc.cell = cell ;
		m_mesh.get_geo( cell, cellGeo );
		m_mesh.list_nodes( cell, itp.nodes );

		const Scalar w = cellGeo.volume() / MeshType::NV ;

		for ( int k = 0 ; k < MeshType::NV ; ++k ) {
			cellGeo.vertexCoords(k, loc.coords) ;
			itp.coeffs[k] = 1. ;

			func( w, loc, itp ) ;

			itp.coeffs[k] = 0. ;
		}
	}
}


template < typename Func >
void FormBuilder::integrate_particle( const Particles& particles, Func func ) const
{
	const size_t n = particles.count() ;

	typename MeshType::Location loc ;
	typename MeshType::Interpolation itp ;
	typename MeshType::Derivatives dc_dx ;

	for ( size_t i = 0 ; i < n ; ++i ) {
		m_mesh.locate( particles.centers().col(i), loc );
		m_mesh.interpolate( loc, itp );
		m_mesh.get_derivatives( loc, dc_dx );

		func( i, particles.volumes()[i], loc, itp, dc_dx ) ;
	}

}



} //d6

#endif
