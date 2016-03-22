#include "DiphasicSimu.hh"

#include "utils/Config.hh"

namespace d6 {


DiphasicSimu::DiphasicSimu(const Config &config, const char *base_dir)
	: Simu( config, base_dir ),
	  m_meshes{ std::unique_ptr<PrimalMesh>(new PrimalMesh( m_config.box, m_config.res, &m_particles.geo() )),
				std::unique_ptr<  DualMesh>(new   DualMesh( m_config.box, m_config.res, &m_particles.geo() ))
			   },
	  m_grains( new Phase( meshes() ) )
{
	m_particles.generate( config, meshes().primal(), *m_scenario );

	m_grains->serializeAllFields( m_config.exportAllFields ) ;

	m_grains->fraction.set_zero();
	m_grains->stresses.set_zero();
	m_grains->velocity.set_zero();
	m_grains->sym_grad.set_zero();
	m_grains->spi_grad.set_zero();
	m_grains->geo_proj.set_zero();
}

DiphasicSimu::~DiphasicSimu()
{

}

void DiphasicSimu::adapt_meshes()
{
	m_meshes.adapt( m_particles, m_grains );
}

void DiphasicSimu::update_fields(const Scalar dt)
{
	m_stats.nNodes = meshes().primal().nNodes() ;

	//! Compute new grid velocities
	m_solver.step( m_config, dt, *m_grains, m_stats, m_rigidBodies, m_rbStresses ) ;
}

void DiphasicSimu::move_particles(const Scalar dt)
{
	m_particles.update( m_config, dt, *m_grains ) ;
}

void DiphasicSimu::dump_fields( unsigned frame ) const
{
	FileInfo dir( FileInfo(m_base_dir).filePath( arg("frame-%1", frame ) ) ) ;
	dir.makePath() ;
	if( ! dir.exists() )
		dir.makeDir() ;

	// Grid
	{
		std::ofstream ofs( dir.filePath("meshes") );
		boost::archive::binary_oarchive oa(ofs);
		oa << meshes() ;
	}
	// Velocity, Stress, Phi
	{
		std::ofstream ofs( dir.filePath("fields") );
		boost::archive::binary_oarchive oa(ofs);
		oa << *m_grains ;
	}
}


} //d6
