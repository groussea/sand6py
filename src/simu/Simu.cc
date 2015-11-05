#include "Simu.hh"
#include "Phase.hh"

#include "Config.hh"
#include "RigidBody.hh"

#include "utils/Log.hh"
#include "utils/File.hh"

#include "geo/Grid.hh"
#include "geo/LevelSet.hh"

#include "geo/Particles.io.hh"
#include <boost/archive/binary_oarchive.hpp>

#include <bogus/Core/Utils/Timer.hpp>
#include <bogus/Core/Eigen/EigenSerialization.hpp>

namespace d6 {


Simu::Simu(const Config &config, const char *base_dir)
	: m_config(config), m_base_dir( base_dir ), m_solver( m_particles )
{
	m_mesh = new MeshImpl( config.box, config.res ) ;

	m_particles.generate( config, mesh() );

	m_grains = new Phase( mesh() ) ;

	// Rigid bodies
//	LevelSet::Ptr plane = LevelSet::make_plane() ;
//	plane->set_origin( Vec(0,0,1) ) ;

//	m_rigidBodies.emplace_back( plane );

	for( unsigned i = 0 ; i < m_rigidBodies.size() ; ++i ) {
		m_rbStresses.emplace_back( mesh() );
		m_rbStresses.back().set_zero() ;
	}
}

Simu::~Simu()
{
	delete m_mesh ;
	delete m_grains  ;
}

void Simu::run()
{
	m_grains->fraction.set_zero();
	m_grains->stresses.set_zero();
	m_grains->velocity.set_zero();
	m_grains->sym_grad.set_zero();
	m_grains->spi_grad.set_zero();

	if( m_config.output )
		dump( 0 ) ;

	for( unsigned frame = 0 ; frame < m_config.nFrames ; ++ frame ) {
		bogus::Timer timer ;
		Log::Info() << "Starting frame " << frame << std::endl ;

		for( unsigned s = 0 ; s < m_config.substeps ; ++ s ) {
			bogus::Timer timer ;
			Log::Verbose() << arg( "Step %1/%2", s+1, m_config.substeps ) << std::endl ;

			step() ;

			Log::Verbose() << arg( "Step done in %1 s", timer.elapsed() ) << std::endl ;
		}

		Log::Info() << arg( "Frame done in %1 s", timer.elapsed() ) << std::endl ;

		if( m_config.output )
			dump( frame+1 ) ;
		m_particles.geo().clear_log();
	}

	Log::Info() << "All done." << std::endl ;

}

void Simu::step()
{
	// TODO adapt mesh

	m_solver.step( m_config, *m_grains, m_rigidBodies, m_rbStresses ) ;

	m_particles.update( m_config, *m_grains ) ;
}

void Simu::dump( unsigned frame ) const
{
	// Dump frame data for viz
	FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame ) ) ) ;
	dir.makePath() ;
	if( ! dir.exists() )
		dir.makeDir() ;

	// Grid
	{
		std::ofstream ofs( dir.filePath("mesh") );
		boost::archive::binary_oarchive oa(ofs);
		oa << m_mesh->derived() ;
	}
	// Velocity, Stress, Phi
	{
		std::ofstream ofs( dir.filePath("fields") );
		boost::archive::binary_oarchive oa(ofs);
		oa << *m_grains ;
	}
	// Particles
	{
		std::ofstream ofs( dir.filePath("particles") );
		boost::archive::binary_oarchive oa(ofs);
		oa << m_particles.geo() ;
	}
}


} //d6
