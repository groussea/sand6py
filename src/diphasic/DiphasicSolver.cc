#include "DiphasicSolver.hh"

#include "DiphasicStepData.hh"
#include "DiphasicFriction.hh"

#include "FluidPhase.hh"
#include "mono/Phase.hh"
#include "mono/PhaseSolver.hh"

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

void printEnergies( const Config &config, Phase &phase, FluidPhase &fluid  )
{

	 ////// STATS

	DynVec volumes ;
	phase.fraction.shape().compute_lumped_mass( volumes );

	PrimalScalarField mphi = phase.fraction ;
	mphi.flatten().array() *= volumes.array() ;

	PrimalVectorField phiu1 = phase.velocity ;
	phiu1.multiply_by( mphi ) ;

	const Scalar Ec1 = (config.alpha()+1)*phase.velocity.flatten().dot(phiu1.flatten()) ;

	mphi.flatten().array() = 1 - phase.fraction.flatten().array() ;
	mphi.flatten().array() *= volumes.array() ;

	PrimalVectorField phiu2 = fluid.velocity ;
	phiu2.multiply_by( mphi ) ;

	const Scalar Ec2 =                    fluid.velocity.flatten().dot(phiu2.flatten()) ;


	mphi.flatten().array() = 1 + config.alpha() * phase.fraction.flatten().array() ;
	mphi.flatten().array() *= volumes.array() ;

	PrimalVectorField betau = fluid.mavg_vel ;
	const Scalar Ecu = fluid.mavg_vel.flatten().dot(betau.flatten()) ;

	std::cout << "Ec1 \t " << 	Ec1 << "\n"
			  << "Ec2 \t " << 	Ec2 << "\n"
			  << "Ect \t " << 	Ec1 + Ec2 << "\n"
			  << "Ecu \t " << 	Ecu
			  << std::endl ;



}

void DiphasicSolver::step(const Config &config, const Scalar dt, Phase &phase, FluidPhase &fluid ) const
{

	DiphasicStepData stepData ;
	stepData.compute( m_particles, config, dt, fluid, phase );

//	printEnergies( config, phase, fluid);
	solve( config, dt, stepData, phase, fluid ) ;
//	printEnergies( config, phase, fluid);
}

void DiphasicSolver::addCohesionContrib (const Config&c, const DiphasicStepData &stepData,
									  DiphasicPrimalData& pbData, DynVec &l ) const
{
	//Cohesion : add \grad{ c phi } to rhs

	DynVec cohe_stress( pbData.H.rows() ) ;
	PhaseSolver::getCohesiveStress( c, stepData.cohesion,
					   stepData.forms.fraction/stepData.forms.volumes,
					   cohe_stress ) ;
	l.head(pbData.m()) -= pbData.G.transpose() * cohe_stress ;
	l.segment(pbData.m(),pbData.r()) -= pbData.H.transpose() * cohe_stress ;

}

void DiphasicSolver::solve(
	const Config& config, const Scalar dt, const DiphasicStepData& stepData ,
	Phase& phase, FluidPhase &fluid ) const
{
	typedef Eigen::SparseMatrix< Scalar > SM ;

	(void) dt ;

	// I Setup

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


	// FIXME avoid copies
	DiphasicPrimalData primal ;
	primal.A = stepData.forms.A ;
	primal.R = stepData.forms.R ;
	primal.M_lumped_inv = stepData.forms.M_lumped_inv ;
	primal.B = stepData.fullGridProj.pressure * ( stepData.forms.B * stepData.fullGridProj.vel ) ;
	primal.C = stepData.fullGridProj.pressure * ( stepData.forms.C * stepData.activeProj.vel ) ;

	const Index m  = primal.m() ;
	const Index r  = primal.r() ;
	const Index p  = primal.p() ;

	DynVec x ( primal.s() ) ;
	x.head(m) = fluid.mavg_vel.flatten() ; // = stepData.forms.M_lumped_inv * rhs ;
	auto w = x.segment( m, r ) ;
	x.segment(m+r,p) = fluid.pressure.flatten() ;

	stepData.primalNodes.field2var( fluid.velocity, w ) ;

	DynVec l ( primal.s() ) ;
	l.setZero() ;
	l.head(m) = rhs ;

	//  Add cohesion forces to rhs
	primal.H = stepData.forms.S.inv_sqrt *
			( stepData.activeProj.stress * ( stepData.forms.H * stepData.activeProj.vel ) ) ;
	primal.G = stepData.forms.S.inv_sqrt *
			( stepData.activeProj.stress * ( stepData.forms.G * stepData.fullGridProj.vel ) ) ;


	addCohesionContrib( config, stepData, primal, l );

	// II Solve unconstrained momentum equation

	SM M ;
	primal.makePenalizedEigenStokesMatrix( M, 1.e-8 );
	DiphasicPrimalData::MInvType M_fac ;
	DiphasicPrimalData::factorize( M, M_fac ) ;

	x = M_fac * l  ;

	// III Friction

	// M d(x,w,p) = (H G)' lambda
	// gamma = H(u+du) + G(w+dw) + q


	const Index n = stepData.dualNodes.count() ;

	// Inertia, mu(I) = \delta_mu * (1./ (1 + I0/I) ), I = dp * sqrt( rho ) * inertia, inertia = |D(U)|/sqrt(p)
	const Scalar I0bar = config.I0 / ( config.grainDiameter * std::sqrt( config.volMass )) ;
	primal.mu = DynArr::Constant( n, config.mu ) +
			config.delta_mu / ( 1. + I0bar / stepData.inertia.max(1.e-12) ) ;

	primal.k = primal.G * x.head(m) + primal.H * x.segment( m, r ) ;

	// Compressability beta(phi)
	{
		const DynArr intBeta = ( config.phiMax*stepData.forms.volumes
								 - stepData.forms.fraction );

		DynVec intBeta_s ( DynVec::Zero( n * SD ) ) ;
		component< SD >( intBeta_s, 0 ).array() = intBeta  * s_sqrt_2_d / dt  ;

		primal.k += ( stepData.forms.S.inv_sqrt * intBeta_s ).cwiseMax(0) ;
	}

	DynVec lambda ;
	stepData.dualNodes.field2var( phase.stresses, lambda ) ;

//	primal.dump( "primal.dd6" ) ;

	DiphasicFrictionSolver::Options options ;
	options.algorithm = DiphasicFrictionSolver::Options::PG_Fac_Red ;
	FrictionSolver::Stats stats ;
	DiphasicFrictionSolver( primal ).solve( options, M, M_fac, x, lambda, stats ) ;

	Log::Verbose() << arg3( "Friction: %1 iterations,\t err= %2,\t time= %3 ",
						   stats.nIterations(), stats.residual(), stats.time() ) << std::endl ;

	// IV  Output
	stepData.dualNodes.var2field( lambda, phase.stresses ) ;

	fluid.pressure.flatten() = x.segment(m+r, p) ;
	fluid.mavg_vel.flatten() = x.head(m) ;

	// U_1 = u + wh
	stepData.primalNodes.var2field( w, fluid.velocity ) ;
	phase.velocity.flatten() = x.head(m) + fluid.velocity.flatten() ;

	// U_2 = u - (alpha+1) phi/(1-phi) wh
	PrimalScalarField ratio( phase.fraction.shape() ) ;
	ratio.flatten() = phase.fraction.flatten().array() / (1. - phase.fraction.flatten().array().min(DiphasicStepData::s_maxPhi) ) ;
	fluid.velocity.multiply_by( ratio ) ;

	fluid.velocity.flatten() = x.head(m) - (config.alpha()+1)*fluid.velocity.flatten() ;

	std::cout << "U " << x.head(m).lpNorm< Eigen::Infinity >() << std::endl ;
	std::cout << "W " << ( fluid.velocity.flatten() - phase.velocity.flatten()).lpNorm< Eigen::Infinity >() << std::endl ;

/*

	const Scalar s = config.fluidVolMass / dt ;
	std::cout << "UMU: " <<   (stepData.forms.M_lumped * x.head(m)).dot( x.head(m) )/s << std::endl ;
	std::cout << "UAU: " <<   (stepData.forms.A * x.head(m)).dot( x.head(m) )/s << std::endl ;
//	std::cout << "U-grad-p: " <<   -( stepData.fullGridProj.vel * stepData.forms.B.transpose() * stepData.fullGridProj.pressure * x.tail(n) ).dot(x.head(m)) << std::endl ;
	std::cout << "U-grad-p: " <<   -( stepData.fullGridProj.vel * stepData.forms.B.transpose() * x.tail(n) ).dot(x.head(m))/s << std::endl ;
	std::cout << "U-dot-g: " <<  x.head(m).dot( stepData.fullGridProj.vel * stepData.forms.externalForces)/s  << std::endl ;
	std::cout << "U-dot-lm: " <<  x.head(m).dot(stepData.forms.linearMomentum)/s  << std::endl ;
*/
}


} // d6
