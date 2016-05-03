#include "DiphasicSolver.hh"

#include "DiphasicStepData.hh"

#include "mono/Phase.hh"

#include "simu/LinearSolver.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include <bogus/Core/Utils/Timer.hpp>
#include <bogus/Core/Block.impl.hpp>

namespace d6 {

DiphasicSolver::DiphasicSolver(const DynParticles &particles)
	: m_particles(particles)
{

}

void DiphasicSolver::step(const Config &config, const Scalar dt, Phase &phase ) const
{

	DiphasicStepData stepData ;
	stepData.compute( m_particles, config, dt, phase );

	solve( config, dt, stepData, phase) ;
}


void DiphasicSolver::solve(
	const Config& config, const Scalar dt, const DiphasicStepData& stepData ,
	Phase& phase ) const
{
	bogus::Timer timer ;

	(void) dt ;

	// Compute rhs of momentum conservation -- gravity + u(t)
	DynVec rhs ;
	{
		DynVec forces = stepData.forms.externalForces ;

		// Inertia
#ifndef  FULL_FEM
		forces += stepData.forms.linearMomentum  ;
#endif

		rhs = stepData.fullGridProj.vel * forces ;
	}

	// Solve unconstrained momentum equation
	DynVec u = stepData.forms.M_lumped_inv * rhs ;
	solveSDP( stepData.forms.A, stepData.forms.M_lumped_inv, rhs, u ) ;

	Log::Debug() << "Linear solve: " << timer.elapsed() << std::endl ;


	DynVec wvh ( WD * stepData.nPrimalNodes() ) ;
	wvh.setZero() ;

	// Output
	PrimalVectorField w( phase.velocity.shape() ) ;
	stepData.primalNodes.var2field( wvh, w ) ;

	// U_1
	phase.velocity.flatten() = u + w.flatten() ;

	// U_2
	PrimalScalarField ratio( phase.fraction.shape() ) ;
	ratio.flatten() = phase.fraction.flatten().array() / (1. - phase.fraction.flatten().array() ) ;
	w.multiply_by( ratio ) ;

	phase.geo_proj.flatten() = u - (config.alpha()+1)*w.flatten() ;
}


} // d6
