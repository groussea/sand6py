#include "PhaseSolver.hh"

#include "Phase.hh"

namespace d6 {

PhaseSolver::PhaseSolver(const DynParticles &particles)
	: m_particles(particles)
{

}


void PhaseSolver::step( const Config &config, const MeshType &mesh, Phase &phase)
{

	phase.velocity.set_constant( Vec(.5,0.,-1) );

}

} //d6
