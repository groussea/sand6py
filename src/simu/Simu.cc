#include "Simu.hh"
#include "Phase.hh"

#include "Scenario.hh"
#include "RigidBody.hh"

#include "utils/Log.hh"
#include "utils/File.hh"
#include "utils/Config.hh"

#include "geo/MeshImpl.hh"
#include "geo/LevelSet.hh"

#include "geo/LevelSet.io.hh"
#include "geo/Particles.io.hh"

#include "utils/serialization.hh"

#include <boost/archive/binary_oarchive.hpp>

#include <bogus/Core/Utils/Timer.hpp>

namespace d6 {


Simu::Simu(const Config &config, const char *base_dir)
	: m_config(config), m_base_dir( base_dir ),
	  m_stats( m_base_dir ),
	  m_scenario( Scenario::parse( config ) ),
	  m_meshes{ std::unique_ptr<PrimalMesh>(new PrimalMesh( config.box, config.res )),
				#ifdef D6_UNSTRUCTURED_DUAL
				std::unique_ptr<  DualMesh>(new   DualMesh( m_particles.geo().centers() ))
				#else
				std::unique_ptr<  DualMesh>(new   DualMesh( config.box, config.res ))
				#endif
			   },
	  m_grains( new Phase( meshes() ) ),
	  m_solver( m_particles )
{

	m_particles.generate( config, meshes().primal(), *m_scenario );

	// Rigid bodies
	m_scenario->add_rigid_bodies( m_rigidBodies ) ;

	for( unsigned i = 0 ; i < m_rigidBodies.size() ; ++i ) {
		m_rbStresses.emplace_back( meshes().primal() );
		m_rbStresses.back().set_zero() ;
	}
}

Simu::~Simu()
{
}

void Simu::run()
{
	m_grains->fraction.set_zero();
	m_grains->stresses.set_zero();
	m_grains->velocity.set_zero();
	m_grains->sym_grad.set_zero();
	m_grains->spi_grad.set_zero();
	m_grains->geo_proj.set_zero();

	if( m_config.output ) {
		dump_particles( 0 ) ;
		dump_fields( 0 ) ;
	}
	m_particles.events().start();

	for( unsigned frame = 0 ; frame < m_config.nFrames ; ++ frame ) {
		bogus::Timer timer ;
		Log::Info() << "Starting frame " << (frame+1) << std::endl ;

		unsigned substeps = m_config.substeps ;
		if( substeps == 0 ) { //Adaptative timestepping
			substeps = std::ceil( std::max(1., m_stats.maxVelocity) / m_config.fps ) ;
		}

		for( unsigned s = 0 ; s < substeps ; ++ s ) {

			m_stats.frameId = frame ;
			m_stats.delta_t = 1./(m_config.fps * substeps) ;

			const Scalar t = m_config.time( frame ) + s * m_stats.delta_t ;
			Log::Verbose() << arg3( "Step %1/%2 \t t=%3 s",
									s+1, substeps,
									t * m_config.units().toSI( Units::Time ) ) << std::endl ;
			// Update external objects (moving boundaries,...)
			m_scenario->update( *this, t, m_stats.delta_t ) ;

			adapt_meshes();

			// Dump particles at last substep of each frame
			if( m_config.output && substeps == s+1 ) {
				dump_particles( frame+1 ) ;
				m_particles.events().clear();
			}
			m_particles.events().start();

			// Proper simulation step
			step( m_stats.delta_t ) ;

			// Log max particle velocity (useful for adaptative timestep )
			m_stats.maxVelocity = m_particles.geo().velocities().leftCols(m_particles.count()).lpNorm< Eigen::Infinity >() ;
			Log::Debug() << "Max particle vel: " << m_stats.maxVelocity << std::endl ;

			m_stats.dump();
		}

		Log::Info() << arg( "Frame done in %1 s", timer.elapsed() ) << std::endl ;

		if( m_config.output ) {
			dump_fields( frame+1 ) ;
		}
	}

	Log::Info() << "All done." << std::endl ;

}

void Simu::adapt_meshes()
{

#ifdef D6_UNSTRUCTURED_DUAL
	m_meshes.m_dual->resize( m_particles.count() ) ;
	m_meshes.m_dual->compute_weights_from_vertices( m_config ) ;

	// ; m_grains->stresses.set_zero() ;
	m_particles.events().replay( m_grains->stresses ) ;

	m_grains->stresses.fit_shape() ;
	m_grains->sym_grad.fit_shape() ;
	m_grains->spi_grad.fit_shape() ;

#endif


}

void Simu::step(const Scalar dt)
{
	bogus::Timer timer ;

	m_stats.nParticles = m_particles.count() ;
	m_stats.nNodes = meshes().primal().nNodes() ;

	//! Compute new grid velocities
	m_solver.step( m_config, dt, *m_grains, m_stats, m_rigidBodies, m_rbStresses ) ;
	const Scalar solverTime = timer.elapsed() ;

	//! Advance particles
	m_particles.update( m_config, dt, *m_grains ) ;

	for( RigidBody& rb: m_rigidBodies ) {
		rb.move( dt );
	}

	m_stats.totalTime = timer.elapsed() ;
	m_stats.advectionTime = timer.elapsed() - solverTime ;

	Log::Verbose() << arg( "Step done in %1 s", m_stats.totalTime ) << std::endl ;
}

void Simu::dump_particles( unsigned frame ) const
{
	FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame ) ) ) ;
	dir.makePath() ;
	if( ! dir.exists() )
		dir.makeDir() ;

	// Particles
	{
		std::ofstream ofs( dir.filePath("particles") );
		boost::archive::binary_oarchive oa(ofs);
		oa << m_particles.geo() ;
	}
	// Objects
	{
		std::ofstream ofs( dir.filePath("objects") );
		boost::archive::binary_oarchive oa(ofs);

		unsigned n = m_rigidBodies.size() ;
		oa << n ;
		LevelSet::register_derived(oa) ;
		for( unsigned i = 0 ; i < n ; ++i ) {
			const LevelSet* ptr = m_rigidBodies[i].levelSetPtr() ;
			oa << ptr ;
		}
	}

	// Log -- save at prevbious frame
	if( frame > 0 ) {
		FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame-1 ) ) );
		std::ofstream ofs( dir.filePath("log") );
		boost::archive::binary_oarchive oa(ofs);
		oa << m_particles.events() ;
	}
}


void Simu::dump_fields( unsigned frame ) const
{
	FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame ) ) ) ;
	dir.makePath() ;
	if( ! dir.exists() )
		dir.makeDir() ;

	// Grid
	{
		std::ofstream ofs( dir.filePath("meshes") );
		boost::archive::binary_oarchive oa(ofs);
		oa << meshes().primal() ;
		oa << meshes().  dual() ;
	}
	// Velocity, Stress, Phi
	{
		std::ofstream ofs( dir.filePath("fields") );
		boost::archive::binary_oarchive oa(ofs);
		oa << *m_grains ;
	}
}


} //d6
