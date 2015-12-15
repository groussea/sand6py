#include "Particles.hh"

#include "MeshImpl.hh"

#include "Tensor.hh"

#include "utils/Log.hh"

#include <bogus/Core/Utils/Timer.hpp>

namespace d6 {

const size_t Particles::s_MAX = 1.e7 ;

Particles::Particles()
	: m_count(0)
{
	resize(s_MAX) ;
}

void Particles::generate(const ScalarExpr &expr, const unsigned nSamples,
						 const MeshType &mesh, const bool alignOnCells,
						 const Vec& initialOri )
{
	bogus::Timer timer ;

	m_count = 0 ;

	// Uniform gen
	typename MeshType::CellGeo cellGeo ;

	for( typename MeshType::CellIterator it = mesh.cellBegin() ; it != mesh.cellEnd() ; ++it ) {
		mesh.get_geo( *it, cellGeo ) ;

		if( alignOnCells && expr( cellGeo.center() ) == 0. )
			continue ;

		Index n = cellGeo.sample_uniform( nSamples, m_count, m_centers, m_frames ) ;
		const Scalar volume = cellGeo.volume() / n ;

		for( size_t i = m_count ; i < m_count+n ; ) {
			if( !alignOnCells && expr( m_centers.col(i) ) == 0. ) {
				-- n ;
				m_centers.col(i) = m_centers.col(m_count+n) ;
				m_frames .col(i) = m_frames .col(m_count+n) ;
			} else ++i ;
		}

		m_volumes.segment( m_count, n ).setConstant( volume ) ;

		m_count += n ;
	}

	m_velocities.leftCols( count() ).setZero() ;

	Vec6 oriCoeffs ;
	{
		Mat oriTensor = Mat::Zero() ;
		oriTensor.diagonal() = initialOri ;
		tensor_view( oriCoeffs ).set( oriTensor) ;
	}

	m_orient.leftCols(m_count).colwise() = oriCoeffs ; // Isotropic ori

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

void Particles::EventLog::log(const Event& event)
{
	LockGuard lock( m_log_mutex ) ;
	m_log.back().push_back( event );
}

} //d6
