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
#include <bogus/Core/Block.io.hpp>


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


	// D(U_g), W(u_g)

	{
		// Velocities gradient D(u) and W(u)

		DynVec int_phiDu = .5 * ( stepData.activeProj.stress * stepData.forms.G * phase.velocity.flatten() ).head( SD * stepData.nDualNodes() ) ;
		DynVec int_phiWu = .5 * stepData.forms.J *  phase.velocity.flatten() ;
		const DynVec int_phi = stepData.forms.fraction.max( 1.e-16 ) ;

		div_compwise<SD>( int_phiDu, int_phi ) ;
		div_compwise<RD>( int_phiWu, int_phi ) ;

		stepData.dualNodes.var2field( int_phiWu, phase.spi_grad ) ;
		stepData.dualNodes.var2field( int_phiDu, phase.sym_grad ) ;
	}

	// D(U_fluid)

	{
		PrimalScalarField volumes(phase.fraction.shape()) ;
		phase.fraction.shape().compute_lumped_mass( volumes.flatten() );

		fluid.sym_grad.flatten() = stepData.forms.D * fluid.velocity.flatten() ;
		fluid.sym_grad.divide_by_positive( volumes ) ;
	}
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
	// Step counter, only useful for dumping friction problem pbData
	static unsigned s_stepId = 0 ;

	typedef Eigen::SparseMatrix< Scalar > SM ;

	bogus::Timer timer ;

	// I Setup

	// Compute rhs of momentum conservation -- gravity + u(t)
	DynVec rhs ;
	{
		DynVec forces = stepData.forms.externalForces ;

		// Inertia
		forces += stepData.forms.linearMomentum  ;

		rhs = stepData.fullGridProj.vel * forces + stepData.forms.dirichletTerm ;
	}


	// FIXME avoid copies
	DiphasicPrimalData primal ;
	primal.A = stepData.forms.A ;
	primal.R = stepData.forms.R ;
	primal.R_visc = stepData.forms.R_visc ;
	primal.F = stepData.fullGridProj.vel * ( stepData.forms.F * stepData.activeProj.vel ) ;
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

	stepData.primalNodes.field2var( fluid.velocity, w ) ; //FIXME wrong warm start

	DynVec l ( primal.s() ) ;
	l.setZero() ;
	l.head(m) = rhs ;
	l.segment(m,r) = stepData.activeProj.vel * stepData.forms.fluctuMomentum ;
	l.segment(m,r) += stepData.activeProj.vel * stepData.forms.F.transpose() * stepData.dirichletVel ;
	l.segment(m+r,p) = stepData.fullGridProj.pressure *
	        stepData.forms.B * stepData.dirichletVel ;

	//  Add cohesion forces to rhs
	primal.H = stepData.forms.S.inv_sqrt *
	        ( stepData.activeProj.stress * ( stepData.forms.H * stepData.activeProj.vel ) ) ;
	primal.G = stepData.forms.S.inv_sqrt *
	        ( stepData.activeProj.stress * ( stepData.forms.G * stepData.fullGridProj.vel ) ) ;


	addCohesionContrib( config, stepData, primal, l );

	// II Solve unconstrained momentum equation

	timer.reset() ;

	SM M ;
	primal.makePenalizedEigenStokesMatrix( M, 1.e-8 );
	DiphasicPrimalData::MInvType M_fac ;
	Log::Debug() << "M assembly: " << timer.elapsed() << std::endl ;
	DiphasicPrimalData::factorize( M, M_fac ) ;
	Log::Debug() << "M factorization: " << timer.elapsed() << std::endl ;

	if( M_fac.block(0).factorization().info() != 0 ) {
		Log::Error() << "Stokes fac failed! "  << std::endl ;
		std::abort() ;
	}

	x = M_fac * l  ;

	Log::Debug() << "Stokes solver time: " << timer.elapsed() << std::endl ;


	// III Friction

	// M d(x,w,p) = (H G)' lambda
	// gamma = H(u+du) + G(w+dw) + q


	const Index n = stepData.dualNodes.count() ;

	// Inertia, mu(I) = \delta_mu * (1./ (1 + I0/I) ), I = dp * sqrt( rho ) * inertia, inertia = |D(U)|/sqrt(p)
	const Scalar I0bar = config.I0 / ( config.grainDiameter * std::sqrt( config.volMass )) ;
	primal.mu = DynArr::Constant( n, config.mu ) +
	        config.delta_mu / ( 1. + I0bar / stepData.inertia.max(1.e-12) ) ;

	primal.k = primal.G * x.head(m) + primal.H * x.segment( m, r )
	        + stepData.forms.S.inv_sqrt * stepData.activeProj.stress *
	          stepData.forms.G * stepData.dirichletVel ;

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

	// Dump problem data if requested
	if( config.dumpPrimalData > 0 && (++s_stepId % config.dumpPrimalData) == 0 ) {
		primal.dump( arg("primal-%1.dd6", s_stepId).c_str() ) ;
	}

	DiphasicFrictionSolver::Options options ;
	if(config.usePG)
		options.algorithm = DiphasicFrictionSolver::Options::PG ;
	else
		options.algorithm = DiphasicFrictionSolver::Options::GS ;
//	options.algorithm = DiphasicFrictionSolver::Options::ADMM ;

	options.tolerance = 1.e-8 ;
	options.useInfinityNorm = config.useInfNorm ;
	options.maxIterations = 1000 ;

	options.useCadoux = false ;

	FrictionSolver::Stats stats ;
	DiphasicFrictionSolver( primal ).solve( options, M_fac, x, lambda, stats ) ;

	Log::Verbose() << arg3( "Friction: %1 iterations,\t err= %2,\t time= %3 ",
	                       stats.nIterations(), stats.residual(), stats.time() ) << std::endl ;

	// IV  Output
	const Scalar sStk = std::sqrt( config.Stokes() ) ;

	stepData.dualNodes.var2field( lambda, phase.stresses ) ;

	fluid.pressure.flatten() = x.segment(m+r, p) ;
	fluid.mavg_vel.flatten() = x.head(m) ;

	// U_1 = u + sStk wh
	stepData.primalNodes.var2field( w*sStk, fluid.velocity ) ;
	phase.velocity.flatten() = x.head(m) + fluid.velocity.flatten() ;

	// U_2 = u - sStk (alpha+1) phi/(1-phi) wh
	PrimalScalarField ratio( phase.fraction.shape() ) ;
	ratio.flatten() = phase.fraction.flatten().array() / (1. - phase.fraction.flatten().array().min(config.phiMax) ) ;
	fluid.velocity.multiply_by( ratio ) ;

	fluid.velocity.flatten() = x.head(m) - (config.alpha()+1)*fluid.velocity.flatten() ;

	Log::Debug() << "U " << x.head(m).lpNorm< Eigen::Infinity >() << std::endl ;
	Log::Debug() << "W " << ( fluid.velocity.flatten() - phase.velocity.flatten()).lpNorm< Eigen::Infinity >() << std::endl ;
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
