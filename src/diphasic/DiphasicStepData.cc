#include "DiphasicStepData.hh"

#include "mono/Phase.hh"
#include "mono/PhaseStepData.hh"

#include "simu/FormBuilder.impl.hh"
#include "simu/DynParticles.hh"

#include "geo/BoundaryInfo.hh"

#include "utils/Config.hh"
#include "utils/Log.hh"

#include <bogus/Core/Utils/Timer.hpp>

namespace d6 {

void DiphasicStepData::assembleMatrices(const Particles &particles,
		const Config &config, const Scalar dt, const DualShape &dShape,
										const PrimalScalarField &phiInt )
{

	const PrimalShape& pShape = phiInt.shape() ;

	PhaseStepData::computeProjectors( config, pShape, dShape, phiInt, primalNodes, dualNodes, 0, activeProj ) ;

}

void DiphasicStepData::compute(const DynParticles& particles,
		const Config &config, const Scalar dt, Phase &phase )
{
	const PrimalShape& pShape = phase.velocity.shape() ;
	const   DualShape& dShape = phase.stresses.shape() ;

	// Transfer particles quantities to grid
	PrimalScalarField intPhiPrimal  ( pShape ) ;
	PrimalVectorField intPhiVel     ( pShape ) ;

	// Not used, but required by particle interface
	DualScalarField intPhiDual    ( dShape ) ;
	DualScalarField intPhiInertia ( dShape ) ;
	DualScalarField intPhiCohesion( dShape ) ;
	DualTensorField intPhiOrient  ( dShape ) ;

	std::vector< bool > activeCells ;

#if defined(FULL_FEM)
	pShape.compute_tpz_mass   ( intPhiPrimal.flatten() ) ;
	dShape.compute_lumped_mass( intPhiDual  .flatten() ) ;
	intPhiVel.set_zero() ;
	intPhiInertia.set_zero() ;
	intPhiCohesion.set_zero() ;
	intPhiOrient.set_zero() ;
	activeCells.assign( pShape.mesh().nCells(), true ) ;
#else
	particles.integratePrimal( activeCells, intPhiPrimal, intPhiVel ) ;
	particles.integrateDual( intPhiDual, intPhiInertia, intPhiOrient, intPhiCohesion ) ;
#endif

	// Compute phi and grad_phi (for visualization purposes )
	PhaseStepData::computePhiAndGradPhi( intPhiPrimal, phase.fraction, phase.grad_phi ) ;

	// Active nodes
	PhaseStepData::computeActiveNodes( activeCells, pShape, dShape, primalNodes, dualNodes ) ;
	Log::Verbose() << "Active nodes: " << nPrimalNodes() << " / " << pShape.nDOF() << std::endl;
	Log::Verbose() << "  Dual nodes: " <<   nDualNodes() << " / " << dShape.nDOF() << std::endl;

	// Bilinear forms matrices
	assembleMatrices( particles.geo(), config, dt, dShape, intPhiPrimal );

	primalNodes.field2var( intPhiVel , forms.phiu ) ;
	dualNodes  .field2var( intPhiDual, forms.fraction ) ;


	// Volumes
	{
		DualScalarField volumes ( dShape ) ;
		dShape.compute_lumped_mass( volumes.flatten() );
		dualNodes.field2var( volumes, forms.volumes ) ;
	}

	// External forces
	PrimalVectorField gravity ( intPhiPrimal.shape() ) ;
	gravity.set_constant( config.gravity );
	gravity.multiply_by( intPhiPrimal ) ;
	gravity.flatten() *= config.volMass ;
	primalNodes.field2var( gravity, forms.externalForces ) ;
}

} // d6
