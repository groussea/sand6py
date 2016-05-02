#include "DiphasicSolver.hh"

#include "DiphasicStepData.hh"

namespace d6 {

DiphasicSolver::DiphasicSolver(const DynParticles &particles)
	: m_particles(particles)
{

}

void DiphasicSolver::step(const Config &config, const Scalar dt, Phase &phase ) const
{

	DiphasicStepData stepData ;
	stepData.compute( m_particles, config, dt, phase );


}




} // d6
