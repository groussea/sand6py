#include "DiphasicSolver.hh"

#include "DiphasicStepData.hh"
#include "DiphasicFriction.hh"

#include "FluidPhase.hh"
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

	// Solve unconstrained momentum equation
	const Index m  = stepData.forms.A.rows() ;
	const Index ma = stepData.forms.R.rows() ;
	const Index n  = stepData.forms.B.rows() ;

	DynVec x ( m+ma+n ) ;
	x.head(m) = fluid.mavg_vel.flatten() ; // = stepData.forms.M_lumped_inv * rhs ;
	x.tail(n) = fluid.pressure.flatten() ;

	const Scalar s = config.fluidVolMass / dt ;
	std::cout << "UMU: " <<   (stepData.forms.M_lumped * x.head(m)).dot( x.head(m) )/s << std::endl ;

	auto w = x.segment( m, ma ) ;
	stepData.primalNodes.field2var( fluid.velocity, w ) ;

	DynVec l ( m+ma+n ) ;
	l.head(m) = rhs ;
	l.segment( m, ma ).setZero() ;
	l.tail(n).setZero();

	// II Solve Linear

	// FIXME avoid copies
	DiphasicPrimalData primal ;
	primal.A = stepData.forms.A ;
	primal.R = stepData.forms.R ;
	primal.M_lumped_inv = stepData.forms.M_lumped_inv ;
	primal.B = stepData.fullGridProj.pressure * ( stepData.forms.B * stepData.fullGridProj.vel ) ;
	primal.C = stepData.fullGridProj.pressure * ( stepData.forms.C * stepData.activeProj.vel ) ;

	SM M ;
	primal.makePenalizedEigenStokesMatrix( M, 1.e-8 );
	DiphasicPrimalData::MInvType M_fac ;
	DiphasicPrimalData::factorize( M, M_fac ) ;

	x = M_fac * l  ;

	// III Friction

	// M d(x,w,p) = (H G)' lambda
	// gamma = H(u+du) + G(w+dw) + q


	const Index nc = stepData.dualNodes.count() ;

	primal.mu = DynVec::Constant( nc, config.mu ) ;

	primal.H = stepData.forms.S.inv_sqrt *
			( stepData.activeProj.stress * ( stepData.forms.H * stepData.activeProj.vel ) ) ;
	primal.G = stepData.forms.S.inv_sqrt *
			( stepData.activeProj.stress * ( stepData.forms.G * stepData.fullGridProj.vel ) ) ;

	primal.k = primal.G * x.head(m) + primal.H * x.segment( m, ma ) ;

	// Compressability beta(phi)
	{
		const DynArr intBeta = ( config.phiMax*stepData.forms.volumes
								 - stepData.forms.fraction );

		DynVec intBeta_s ( DynVec::Zero( nc * SD ) ) ;
		component< SD >( intBeta_s, 0 ).array() = intBeta  * s_sqrt_2_d / dt  ;

		primal.k += ( stepData.forms.S.inv_sqrt * intBeta_s ).cwiseMax(0) ;
	}

	DynVec lambda ;
	stepData.dualNodes.field2var( phase.stresses, lambda ) ;

	primal.dump( "primal.dd6" ) ;

	DiphasicFrictionSolver::Options options ;
	FrictionSolver::Stats stats ;
	DiphasicFrictionSolver( primal ).solve( options, M, M_fac, x, lambda, stats ) ;

	Log::Verbose() << arg3( "Friction: %1 iterations,\t err= %2,\t time= %3 ",
						   stats.nIterations(), stats.residual(), stats.time() ) << std::endl ;

	// IV  Output
	stepData.dualNodes.var2field( lambda, phase.stresses ) ;

	fluid.pressure.flatten() = x.tail(n) ;
	fluid.mavg_vel.flatten() = x.head(m) ;

	// U_1 = u + wh
	stepData.primalNodes.var2field( w, fluid.velocity ) ;
	phase.velocity.flatten() = x.head(m) + fluid.velocity.flatten() ;

	// U_2 = u - (alpha+1) phi/(1-phi) wh
	PrimalScalarField ratio( phase.fraction.shape() ) ;
	ratio.flatten() = phase.fraction.flatten().array() / (1. - phase.fraction.flatten().array().min(config.phiMax) ) ;
	fluid.velocity.multiply_by( ratio ) ;

	fluid.velocity.flatten() = x.head(m) - (config.alpha()+1)*fluid.velocity.flatten() ;

	std::cout << "U " << x.head(m).lpNorm< Eigen::Infinity >() << std::endl ;
	std::cout << "W " << ( fluid.velocity.flatten() - phase.velocity.flatten()).lpNorm< Eigen::Infinity >() << std::endl ;

/*
	std::cout << "UMU: " <<   (stepData.forms.M_lumped * x.head(m)).dot( x.head(m) )/s << std::endl ;
	std::cout << "UAU: " <<   (stepData.forms.A * x.head(m)).dot( x.head(m) )/s << std::endl ;
//	std::cout << "U-grad-p: " <<   -( stepData.fullGridProj.vel * stepData.forms.B.transpose() * stepData.fullGridProj.pressure * x.tail(n) ).dot(x.head(m)) << std::endl ;
	std::cout << "U-grad-p: " <<   -( stepData.fullGridProj.vel * stepData.forms.B.transpose() * x.tail(n) ).dot(x.head(m))/s << std::endl ;
	std::cout << "U-dot-g: " <<  x.head(m).dot( stepData.fullGridProj.vel * stepData.forms.externalForces)/s  << std::endl ;
	std::cout << "U-dot-lm: " <<  x.head(m).dot(stepData.forms.linearMomentum)/s  << std::endl ;
*/
}


} // d6
