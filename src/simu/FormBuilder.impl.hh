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
	typename Linear<MeshImpl>::Interpolation itp ;
	typename Linear<MeshImpl>::Derivatives dc_dx ;

	Linear<MeshImpl> shape( m_mesh ) ;

	for ( auto qpIt = shape.qpIterator(cellBegin) ; qpIt != shape.qpIterator(cellEnd) ; ++qpIt )
	{
		qpIt.locate( loc ) ;
		shape.interpolate( loc, itp );
		shape.get_derivatives( loc, dc_dx );

		func( qpIt.weight(), loc, itp, dc_dx ) ;
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
	typename Linear<MeshImpl>::Interpolation itp ;

	typename MeshType::CellGeo cellGeo ;

	itp.coeffs.setZero() ;

	for( CellIterator it = cellBegin ; it !=  cellEnd ; ++it )
	{
		const typename MeshType::Cell& cell = *it ;
		loc.cell = cell ;
		m_mesh.get_geo( cell, cellGeo );
		m_mesh.template shaped<Linear>().list_nodes( loc, itp.nodes );

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
	typename Linear<MeshImpl>::Interpolation itp ;
	typename Linear<MeshImpl>::Derivatives dc_dx ;

	for ( size_t i = 0 ; i < n ; ++i ) {
		m_mesh.locate( particles.centers().col(i), loc );
		m_mesh.template shaped<Linear>().interpolate( loc, itp );
		m_mesh.template shaped<Linear>().get_derivatives( loc, dc_dx );

		func( i, particles.volumes()[i], loc, itp, dc_dx ) ;
	}

}



} //d6

#endif
