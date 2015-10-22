#include "Simu.hh"
#include "Phase.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"
#include "utils/File.hh"

#include "geo/Grid.hh"

#include "geo/Particles.io.hh"
#include <boost/archive/binary_oarchive.hpp>

#include <bogus/Core/Utils/Timer.hpp>
#include <bogus/Core/Eigen/EigenSerialization.hpp>

namespace d6 {


Simu::Simu(const Config &config, const char *base_dir)
	: m_config(config), m_base_dir( base_dir )
{
	m_mesh = new MeshImpl( config.box, config.res ) ;

	m_particles.generate( config, mesh() );

	m_phase = new Phase( mesh() ) ;
}

Simu::~Simu()
{
	delete m_mesh ;
	delete m_phase  ;
}

void Simu::run()
{
	m_phase->fraction.set_zero();
	m_phase->stresses.set_zero();
	m_phase->velocity.set_zero();
	m_phase->sym_grad.set_zero();
	m_phase->spi_grad.set_zero();

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

		dump( frame+1 ) ;
	}

	Log::Info() << "All done." << std::endl ;

}

void Simu::step()
{
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
		oa << *m_phase ;
	}
	// Particles
	{
		std::ofstream ofs( dir.filePath("particles") );
		boost::archive::binary_oarchive oa(ofs);
		oa << m_particles.geo() ;
	}
}


} //d6
