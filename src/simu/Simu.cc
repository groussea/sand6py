#include "Simu.hh"
#include "Phase.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include "geo/Grid.hh"

#include <bogus/Core/Utils/Timer.hpp>

namespace d6 {


Simu::Simu(const Config &config)
	: m_config(config)
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

		dump() ;
	}

	Log::Info() << "All done." << std::endl ;

}

void Simu::step()
{

}

void Simu::dump() const
{
	// Dump frame data for viz

	// Grid
	// Velocity, Stress, Phi
	// Particles
}


} //d6
