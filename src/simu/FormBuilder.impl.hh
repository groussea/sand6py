#ifndef D6_FORM_BUILDER_IMPL_HH
#define D6_FORM_BUILDER_IMPL_HH

#include "FormBuilder.hh"

#include "geo/Grid.hh"
#include "geo/Particles.hh"

namespace d6 {

template < typename Func >
void FormBuilder::integrate_qp( const typename MeshType::Cells& cells, Func func ) const
{

	typename MeshType::Location loc ;
	typename MeshType::Interpolation itp ;
	typename MeshType::Derivatives dc_dx ;

	Eigen::Matrix< Scalar, MeshType::NC, MeshType::NQ > qp ;
	Eigen::Matrix< Scalar,            1, MeshType::NQ > qp_weights ;

	typename MeshType::CellGeo cellGeo ;

	for( const typename MeshType::Cell& cell : cells )
	{
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
void FormBuilder::integrate_particle( const Particles& particles, Func func ) const
{
	typename MeshType::Location loc ;
	typename MeshType::Interpolation itp ;
	typename MeshType::Derivatives dc_dx ;

	for ( unsigned i = 0 ; i < particles.count() ; ++i ) {
		m_mesh.locate( particles.centers().col(i), loc );
		m_mesh.interpolate( loc, itp );
		m_mesh.get_derivatives( loc, dc_dx );

		func( particles.volumes()[i], loc, itp, dc_dx ) ;
	}

}



} //d6

#endif
