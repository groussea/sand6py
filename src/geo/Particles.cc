#include "Particles.hh"

#include "Grid.hh"

#include "utils/Log.hh"

#include <bogus/Core/Utils/Timer.hpp>

namespace d6 {

const size_t Particles::s_MAX = 1.e6 ;

Particles::Particles()
	: m_count(0)
{
	resize(s_MAX) ;
}

void Particles::generate(const ScalarExpr &expr, const unsigned nSamples, const MeshType &mesh)
{
	bogus::Timer timer ;

	m_count = 0 ;

	// Uniform gen
	typename MeshType::CellGeo cellGeo ;

	for( typename MeshType::CellIterator it = mesh.cellBegin() ; it != mesh.cellEnd() ; ++it ) {
		mesh.get_geo( *it, cellGeo ) ;

		if( expr( cellGeo.center() ) == 0. )
			continue ;

		const Index n = cellGeo.sample_uniform( nSamples, m_count, m_centers, m_frames ) ;

		const Scalar volume = cellGeo.volume() / n ;
		m_volumes.segment( m_count, n ).setConstant( volume ) ;

		m_orient.block( 0, m_count, 1, n ).setConstant( 3 ) ; // Isotropic ori

		m_count += n ;
	}

	m_velocities.leftCols( count() ).setZero() ;

	Log::Verbose() << arg( "Generated %1 particles in %2 s ", m_count, timer.elapsed() ) << std::endl ;
}

void Particles::resize(size_t n)
{
	m_volumes.resize( n );

	m_centers.resize( 3, n);
	m_velocities.resize( 3, n);

	m_frames.resize( 6, n);
	m_orient.resize( 6, n);
}

} //d6
